"""End-to-end contract tests for current-scene video export."""

from __future__ import annotations

import shutil
import os
import subprocess
import time
from pathlib import Path

import pytest


def _wait_for_video(rest, *, timeout: float = 30.0) -> dict:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        response = rest.client.get("/api/video/export/status")
        assert response.status_code == 200, response.text
        status = response.json()
        if not status["running"]:
            return status
        time.sleep(0.02)
    raise AssertionError("scene video export did not finish")


def _assert_finalized_mp4(path: Path) -> None:
    assert path.is_file()
    data = path.read_bytes()
    assert len(data) > 1_024
    assert b"ftyp" in data[:32]
    assert b"mdat" in data
    assert b"moov" in data


def _decoded_frame_hashes(path: Path) -> list[str]:
    ffmpeg = shutil.which("ffmpeg")
    assert ffmpeg is not None
    result = subprocess.run(
        [
            ffmpeg,
            "-v",
            "error",
            "-i",
            str(path),
            "-map",
            "0:v:0",
            "-f",
            "framemd5",
            "-",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = [
        line for line in result.stdout.splitlines()
        if line and not line.startswith("#")
    ]
    return [row.rsplit(",", 1)[-1].strip() for row in rows]


def test_scene_video_export_reports_requested_frame_range(rest, tmp_path: Path) -> None:
    original_frame = rest.client.get("/frame/current").json()["frame"]
    output_path = tmp_path / "reader-scene.mp4"

    response = rest.client.post(
        "/api/video/export",
        json={
            "output_path": str(output_path),
            "start_frame": 0,
            "end_frame": 11,
            "frame_step": 1,
            "frames_per_second": 6,
        },
    )
    assert response.status_code == 202, response.text
    started = response.json()
    assert started["running"] is True
    assert started["frames_total"] == 12
    assert started["frames_per_second"] == 6

    health = rest.client.get("/health", timeout=2.0)
    assert health.status_code == 200, health.text
    assert health.json()["ok"] is True

    completed = _wait_for_video(rest)
    assert completed["ok"] is True
    assert completed["state"] == "completed"
    assert completed["frames_written"] == 12
    assert completed["last_frame"] == 11
    assert completed["frames_per_second"] == 6
    assert completed["file_size_bytes"] == output_path.stat().st_size
    assert completed["width"] >= 2 and completed["width"] % 2 == 0
    assert completed["height"] >= 2 and completed["height"] % 2 == 0
    assert rest.client.get("/frame/current").json()["frame"] == original_frame
    _assert_finalized_mp4(output_path)


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="ffmpeg is not installed")
def test_scene_video_export_decodes_exact_frame_range(rest, tmp_path: Path) -> None:
    output_path = tmp_path / "reader-scene-decoded.mp4"
    response = rest.client.post(
        "/api/video/export",
        json={
            "output_path": str(output_path),
            "start_frame": 0,
            "end_frame": 11,
            "frame_step": 1,
            "frames_per_second": 6,
        },
    )
    assert response.status_code == 202, response.text
    completed = _wait_for_video(rest)
    assert completed["state"] == "completed"

    hashes = _decoded_frame_hashes(output_path)
    assert len(hashes) == 12
    assert len(set(hashes)) > 1


def test_scene_video_stop_finalizes_partial_file(rest, tmp_path: Path) -> None:
    frame_count = rest.client.get("/frame/current").json()["count"]
    assert frame_count >= 12
    output_path = tmp_path / "reader-scene-stopped.mp4"

    response = rest.client.post(
        "/api/video/export",
        json={
            "output_path": str(output_path),
            "start_frame": 0,
            "end_frame": frame_count - 1,
            "frame_step": 1,
            "frames_per_second": 24,
        },
    )
    assert response.status_code == 202, response.text

    assert not rest.client.get("/ui/state").json()["context"]["editable"]
    blocked_context = rest.client.post("/ui/context", json={"action": "clear_selection"})
    assert blocked_context.status_code == 409, blocked_context.text

    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline:
        status = rest.client.get("/api/video/export/status").json()
        if status["running"] and status["frames_written"] > 0:
            break
        time.sleep(0.01)
    else:
        raise AssertionError("scene video export never produced a stoppable frame")

    stopped = rest.client.post("/api/video/export/stop", json={})
    assert stopped.status_code == 202, stopped.text
    assert stopped.json()["state"] == "stopping"

    completed = _wait_for_video(rest)
    assert completed["ok"] is True
    assert completed["state"] == "stopped"
    assert rest.client.get("/ui/state").json()["context"]["editable"]
    assert 0 < completed["frames_written"] < completed["frames_total"]
    assert completed["last_frame"] == completed["frames_written"] - 1
    assert completed["file_size_bytes"] == output_path.stat().st_size
    _assert_finalized_mp4(output_path)


def test_scene_video_export_recovers_from_frame_interference(
    rest, tmp_path: Path
) -> None:
    frame_count = rest.client.get("/frame/current").json()["count"]
    end_frame = min(frame_count - 1, 79)
    assert end_frame >= 11
    output_path = tmp_path / "reader-scene-interference.mp4"

    response = rest.client.post(
        "/api/video/export",
        json={
            "output_path": str(output_path),
            "start_frame": 0,
            "end_frame": end_frame,
            "frame_step": 1,
            "frames_per_second": 24,
        },
    )
    assert response.status_code == 202, response.text

    interference_requests = 0
    while interference_requests < 20:
        status = rest.client.get("/api/video/export/status").json()
        if not status["running"]:
            break
        frame = end_frame if interference_requests % 2 == 0 else 0
        moved = rest.client.post("/frame/set", json={"frame": frame})
        assert moved.status_code == 204, moved.text
        interference_requests += 1

    assert interference_requests > 0
    completed = _wait_for_video(rest, timeout=60.0)
    assert completed["state"] == "completed"
    assert completed["frames_written"] == end_frame + 1
    assert completed["last_frame"] == end_frame
    _assert_finalized_mp4(output_path)


@pytest.mark.skipif(os.environ.get("H5READER_EXPECT_DFT_AND_ML") != "1",
                    reason="requires the DFT fixture and working prediction runtime")
def test_removing_prediction_pin_during_export_loads_orca(rest, tmp_path):
    rest.client.post("/selection/pick", json={"atom": 0}).raise_for_status()
    response = rest.client.post("/dashboard/metric", json={
        "descriptor_id": "ml:experimental_shielding_t2", "anchor": {"atom": 0},
        "modes": ["static.tensor"],
    })
    response.raise_for_status()
    signal = response.json()["id"]
    removed = False
    try:
        output = tmp_path / "switch-to-orca.mp4"
        response = rest.client.post("/api/video/export", json={
            "output_path": str(output), "start_frame": 40, "end_frame": 43,
            "frames_per_second": 1,
        })
        assert response.status_code == 202, response.text
        status = rest.client.get("/api/video/export/status").json()
        assert status["running"] and status["frames_written"] == 0
        rest.client.post("/dashboard/metric/remove", json={"id": signal}).raise_for_status()
        removed = True
        completed = _wait_for_video(rest)
        assert completed["state"] == "completed", completed
        assert completed["frames_written"] == 4
        _assert_finalized_mp4(output)
    finally:
        if rest.client.get("/api/video/export/status").json()["running"]:
            rest.client.post("/api/video/export/stop", json={}).raise_for_status()
            _wait_for_video(rest)
        if not removed:
            rest.client.post("/dashboard/metric/remove", json={"id": signal}).raise_for_status()
