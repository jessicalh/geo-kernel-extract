"""Visible tensors must survive video capture and changes of selection."""

import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from PIL import Image
import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("H5READER_EXPECT_AUTO_ML") != "1",
    reason="requires the packaged no-ORCA prediction fixture",
)


def post(rest, route, body):
    response = rest.client.post(route, json=body)
    response.raise_for_status()
    return response


def ml_state(rest):
    return rest.client.get("/ui/state").json()["experimentalShieldingMl"]


def wait_tensor(rest, atom, frame):
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        state = ml_state(rest)
        tensor = state["tensorDisplay"]
        if tensor["atom"] == atom and tensor["frame"] == frame and tensor["displayed"]:
            return tensor
        time.sleep(0.05)
    pytest.fail(f"tensor did not appear: {state}")


def wait_video(rest):
    deadline = time.monotonic() + 180
    while time.monotonic() < deadline:
        status = rest.client.get("/api/video/export/status").json()
        if not status["running"]:
            return status
        time.sleep(0.05)
    pytest.fail(f"video did not finish: {status}")


def export(rest, path, first, last):
    post(rest, "/api/video/export", {
        "output_path": str(path), "start_frame": first, "end_frame": last,
        "frames_per_second": 1,
    })


def decode_frame(video, index, output):
    subprocess.run([
        shutil.which("ffmpeg"), "-v", "error", "-n", "-ss", str(index),
        "-i", str(video), "-frames:v", "1", str(output),
    ], check=True, capture_output=True)
    with Image.open(output) as image:
        return np.asarray(image.convert("RGB"), dtype=np.int16)


def shielding_group(rest):
    tree = rest.client.get("/inspector/tree").json()
    return next(item for item in tree if item["field"] == "Shielding tensor (Predicted)")


def test_pinned_tensor_keeps_its_identity_and_checkbox_without_focus(rest):
    post(rest, "/selection/pick", {"atom": 16})
    wait_tensor(rest, 16, 0)
    signal = post(rest, "/dashboard/metric", {
        "descriptor_id": "ml:experimental_shielding_t2", "anchor": {"atom": 16},
        "modes": ["static.tensor"],
    }).json()["id"]
    first = shielding_group(rest)
    identity = next(item["value"] for item in first["children"] if item["field"] == "Atom")
    assert identity
    try:
        for atom in (17, None):
            if atom is None:
                post(rest, "/selection/clear", {})
            else:
                post(rest, "/selection/pick", {"atom": atom})
            wait_tensor(rest, 16, 0)
            group = shielding_group(rest)
            assert group["checked"]
            assert next(item["value"] for item in group["children"] if item["field"] == "Atom") == identity
            for visible in (False, True):
                post(rest, "/overlay", {"name": "shielding", "visible": visible})
                assert shielding_group(rest)["checked"] is visible
                assert ml_state(rest)["tensorDisplay"]["displayed"] is visible
    finally:
        post(rest, "/dashboard/metric/remove", {"id": signal})
    assert not ml_state(rest)["tensorDisplay"]["active"]
    assert all(item["field"] != "Shielding tensor (Predicted)"
               for item in rest.client.get("/inspector/tree").json())


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="ffmpeg is required to inspect video pixels")
def test_predicted_video_matches_independently_prepared_frames(rest, tmp_path):
    post(rest, "/selection/pick", {"atom": 16})
    wait_tensor(rest, 16, 0)
    post(rest, "/filter", {"residues": [0]})
    post(rest, "/camera/inspect_atom", {"atom": 16, "distance": 13})
    try:
        video = tmp_path / "predicted.mp4"
        export(rest, video, 0, 3)
        completed = wait_video(rest)
        assert completed["state"] == "completed", completed
        assert completed["frames_written"] == 4
        for frame in range(4):
            post(rest, "/frame/set", {"frame": frame})
            wait_tensor(rest, 16, frame)
            reference_video = tmp_path / f"reference-{frame}.mp4"
            export(rest, reference_video, frame, frame)
            assert wait_video(rest)["state"] == "completed"
            actual = decode_frame(video, frame, tmp_path / f"actual-{frame}.png")
            expected = decode_frame(reference_video, 0, tmp_path / f"expected-{frame}.png")
            # Gold/pink arrow interiors, excluding the molecule's yellow sulfur,
            # red oxygen and blue nitrogen. Tolerate H.264 quantisation at edges.
            red, green, blue = expected.transpose(2, 0, 1)
            gold = (red > 140) & (green > 80) & (green < 0.85 * red) & (blue < 80)
            pink = (red > 140) & (blue > 90) & (green < 0.7 * np.minimum(red, blue))
            mask = gold | pink
            assert mask.sum() > 100, "reference must contain visible tensor arrows"
            difference = np.abs(actual - expected).mean(axis=2)[mask]
            assert np.median(difference) < 10, (frame, np.median(difference))
            assert np.percentile(difference, 95) < 45, (frame, np.percentile(difference, 95))
    finally:
        post(rest, "/filter", {"residues": []})


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="ffmpeg is required to inspect video pixels")
def test_video_waits_for_the_new_trajectory_envelope(rest, tmp_path):
    post(rest, "/selection/pick", {"atom": 16})
    post(rest, "/filter", {"residues": [0]})
    post(rest, "/camera/inspect_atom", {"atom": 16, "distance": 13})
    try:
        post(rest, "/overlay", {"name": "trajectory", "visible": True})
        video = tmp_path / "envelope-immediate.mp4"
        export(rest, video, 0, 0)
        result = wait_video(rest)
        assert result["state"] == "completed", result
        assert result["frames_written"] == 1
        reference = tmp_path / "envelope-settled.mp4"
        export(rest, reference, 0, 0)
        assert wait_video(rest)["state"] == "completed"
        actual = decode_frame(video, 0, tmp_path / "envelope-immediate.png")
        expected = decode_frame(reference, 0, tmp_path / "envelope-settled.png")
        assert np.percentile(np.abs(actual - expected), 99) < 15
    finally:
        post(rest, "/overlay", {"name": "trajectory", "visible": False})
        post(rest, "/filter", {"residues": []})


def test_hidden_tensor_export_does_not_queue_predictions(rest, tmp_path):
    post(rest, "/selection/pick", {"atom": 16})
    wait_tensor(rest, 16, 0)
    post(rest, "/overlay", {"name": "shielding", "visible": False})
    try:
        export(rest, tmp_path / "hidden.mp4", 4, 5)
        assert wait_video(rest)["state"] == "completed"
        assert not ml_state(rest)["inferenceRunning"]
        assert ml_state(rest)["tensorDisplay"]["resident"]
    finally:
        post(rest, "/overlay", {"name": "shielding", "visible": True})


def test_stop_while_prediction_is_pending(rest, tmp_path):
    post(rest, "/selection/pick", {"atom": 16})
    wait_tensor(rest, 16, 0)
    export(rest, tmp_path / "stopped.mp4", 10, 13)
    assert ml_state(rest)["inferenceRunning"]
    before = rest.client.get("/api/video/export/status").json()
    assert before["running"] and before["frames_written"] == 0
    assert rest.client.get("/health", timeout=2).status_code == 200
    post(rest, "/api/video/export/stop", {})
    stopped = wait_video(rest)
    assert stopped["frames_written"] == 0
    assert stopped["state"] == "failed"
    assert stopped["error"] == "video export stopped before the first frame was written"
    wait_tensor(rest, 16, 0)
    assert rest.client.get("/api/video/export/status").json() == stopped


@pytest.mark.skipif(os.environ.get("H5READER_EXPECT_ML_PROCESS_FAILURE") != "1",
                    reason="requires an inference helper configured to fail")
def test_failed_prediction_and_cached_failure_end_export(rest, tmp_path):
    post(rest, "/selection/pick", {"atom": 16})
    deadline = time.monotonic() + 30
    while ml_state(rest)["inferenceRunning"]:
        assert time.monotonic() < deadline
        time.sleep(0.05)
    assert not ml_state(rest)["tensorDisplay"]["resident"]
    for attempt in range(2):
        export(rest, tmp_path / f"failed-{attempt}.mp4", 1, 1)
        result = wait_video(rest)
        assert result["state"] == "failed", result
        assert "Shielding prediction could not be displayed" in result["error"]
        assert result["frames_written"] == 0
        assert not Path(result["output_path"]).exists()


@pytest.mark.skipif(os.environ.get("H5READER_TEST_VIDEO_SHUTDOWN") != "1",
                    reason="run separately because it closes the test Reader")
def test_shutdown_while_prediction_is_pending(rest, tmp_path):
    post(rest, "/selection/pick", {"atom": 16})
    wait_tensor(rest, 16, 0)
    export(rest, tmp_path / "shutdown.mp4", 20, 23)
    assert ml_state(rest)["inferenceRunning"]
    post(rest, "/shutdown", {})
    assert rest.process.wait(timeout=30) == 0
