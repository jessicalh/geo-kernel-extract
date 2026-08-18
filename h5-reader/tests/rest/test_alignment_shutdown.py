"""Regression for graceful shutdown during a scientific alignment export."""

from __future__ import annotations

import os
import threading
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import httpx
import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("H5READER_ALIGNMENT_EXPORT_ENABLED") != "1",
    reason="alignment shutdown integration uses its dedicated target-free fixture",
)


def _export_request(
    base_url: str, output_root: Path, started: threading.Event
) -> httpx.Response:
    with httpx.Client(base_url=base_url, timeout=600.0) as client:
        started.set()
        return client.post(
            "/api/alignment/export",
            json={"output_root": str(output_root), "apply_display": False},
        )


def _wait_until_running(rest, future: Future[httpx.Response]) -> None:
    deadline = time.monotonic() + 30.0
    while time.monotonic() < deadline:
        response = rest.client.get("/api/alignment/export/status")
        assert response.status_code == 200, response.text
        if response.json()["running"]:
            return
        if future.done():
            pytest.fail("alignment export completed before its running state was observable")
        time.sleep(0.01)
    pytest.fail("alignment export did not enter its running state")


def test_shutdown_flushes_both_responses_before_teardown(rest, tmp_path: Path) -> None:
    output_root = tmp_path / "alignment"
    output_root.mkdir()
    started = threading.Event()

    with ThreadPoolExecutor(max_workers=1) as executor:
        export_future = executor.submit(
            _export_request, rest.base_url, output_root, started
        )
        assert started.wait(timeout=5.0)
        _wait_until_running(rest, export_future)

        shutdown = rest.client.post("/shutdown", timeout=60.0)
        export = export_future.result(timeout=120.0)

    assert shutdown.status_code == 204, shutdown.text
    assert export.status_code == 409, export.text
    assert export.json()["cancelled"] is True
    assert not any(output_root.iterdir())

    rest.process.wait(timeout=120.0)
    assert rest.process.returncode == 0
