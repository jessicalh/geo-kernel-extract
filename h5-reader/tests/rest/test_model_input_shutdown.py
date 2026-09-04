"""Graceful shutdown while Reader is exporting model inputs."""

from __future__ import annotations

import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import httpx
import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("H5READER_MODEL_INPUT_EXPORT_ENABLED") != "1",
    reason="model-input export uses an explicitly configured complete fixture",
)


def _export(base_url: str, output: Path, started: threading.Event) -> httpx.Response:
    with httpx.Client(base_url=base_url, timeout=300.0) as client:
        started.set()
        return client.post(
            "/api/model-input/export",
            json={"output_directory": str(output)},
        )


def test_shutdown_cancels_export_and_removes_partial_output(
    rest, tmp_path: Path
) -> None:
    before_state = rest.client.get("/ui/state").json()
    before_position = rest.client.post(
        "/positions", json={"atoms": [0], "frame": 0}
    ).json()["positions"][0]["position"]
    output = tmp_path / "export"
    output.mkdir()
    started = threading.Event()

    with ThreadPoolExecutor(max_workers=1) as executor:
        export_result = executor.submit(_export, rest.base_url, output, started)
        assert started.wait(timeout=5.0)

        deadline = time.monotonic() + 5.0
        while True:
            active = rest.client.post("/api/model-input/export", json={})
            if active.status_code == 409:
                break
            assert active.status_code == 400, active.text
            assert not export_result.done()
            assert time.monotonic() < deadline
            time.sleep(0.01)

        blocked_load = rest.client.post(
            "/api/run/load", json={"path": os.environ["H5READER_REST_FIXTURE"]}
        )
        assert blocked_load.status_code == 409, blocked_load.text
        assert blocked_load.json()["error"] == "another Reader operation is running"
        assert rest.client.get("/ui/state").json() == before_state
        after_position = rest.client.post(
            "/positions", json={"atoms": [0], "frame": 0}
        ).json()["positions"][0]["position"]
        assert after_position == before_position

        health = rest.client.get("/health")
        assert health.status_code == 200, health.text
        assert health.json()["ok"] is True

        with httpx.Client(base_url=rest.base_url, timeout=30.0) as client:
            shutdown = client.post("/shutdown")
        export = export_result.result(timeout=30.0)

    assert shutdown.status_code == 204, shutdown.text
    assert export.status_code == 409, export.text
    assert "cancelled" in export.json()["error"]
    assert not any(output.iterdir())

    rest.process.wait(timeout=30.0)
    assert rest.process.returncode == 0
