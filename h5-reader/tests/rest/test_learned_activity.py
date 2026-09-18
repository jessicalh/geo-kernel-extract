"""Captured hidden tensors use the molecule's normal display and lifecycle."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path

import numpy as np
import pytest


def fixture_capture(rest):
    path = Path(os.environ["H5READER_REST_FIXTURE"])
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest["kind"] == "trajectory":
        root = path.parent / manifest["trajectory"]["extraction_dir"]
        directory = sorted((root / "npys").glob("frame_*"))[0]
        original = int(directory.name.removeprefix("frame_"))
        raw = np.load(directory / "pos.npy").astype(float)
    else:
        directory = path.parent / manifest["single_pose"]["pose_dir"]
        original = 0
        raw = np.load(directory / "pos.npy").astype(float)
    assert rest.client.get("/protein/atoms").json()["count"] == len(raw)
    atoms = [0, 1, 2]
    tensor = np.diag([-1.0, -1.0, 2.0])
    six = tensor[[0, 0, 0, 1, 1, 2], [0, 1, 2, 1, 2, 2]]
    capture = {
        "schema_version": 1,
        "model": "synthetic test only",
        "coordinate_frame": "extraction",
        "atom_count": len(raw),
        "atom_indices": atoms,
        "frames": [
            {
                "frame_index": original,
                "positions": raw[atoms].tolist(),
                "channels": [
                    {
                        "name": "test / 2e 0",
                        "tensors": [six.tolist(), (0.5 * six).tolist(), [0] * 6],
                    }
                ],
            }
        ],
    }
    return capture, raw, tensor


def load_capture(rest, tmp_path, capture):
    path = tmp_path / "activity.json"
    path.write_text(json.dumps(capture), encoding="utf-8")
    return rest.client.post("/api/learned-activity/load", json={"path": str(path)})


def matrix(six):
    xx, xy, xz, yy, yz, zz = six
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def test_activity_scale_rotation_missing_frame_and_reload(rest, tmp_path):
    capture, raw, tensor = fixture_capture(rest)
    try:
        response = load_capture(rest, tmp_path, capture)
        assert response.status_code == 200, response.text
        state = response.json()
        assert state["reference_magnitude"] == pytest.approx(2.0)
        assert [s["peak_radius"] for s in state["samples"]] == pytest.approx(
            [1.5, 0.75, 0.0]
        )
        for kind in ["all_atom_fit", "backbone_fit"]:
            assert (
                rest.client.post("/transform", json={"kind": kind}).status_code == 204
            )
            state = rest.client.get("/api/learned-activity").json()
            response = rest.client.post(
                "/positions", json={"frame": 0, "atoms": list(range(len(raw)))}
            )
            displayed = np.array(
                [item["position"] for item in response.json()["positions"]]
            )
            u, _, vt = np.linalg.svd(
                (raw - raw.mean(0)).T @ (displayed - displayed.mean(0))
            )
            rotation = vt.T @ u.T
            assert np.linalg.det(rotation) == pytest.approx(1.0)
            np.testing.assert_allclose(
                matrix(state["samples"][0]["tensor"]),
                rotation @ tensor @ rotation.T,
                atol=1e-7,
            )
            np.testing.assert_allclose(
                state["samples"][0]["center"], displayed[0], atol=1e-7
            )
        assert (
            rest.client.post("/api/learned-activity", json={"radius": 2.0}).status_code
            == 200
        )
        assert rest.client.get("/api/learned-activity").json()["samples"][1][
            "peak_radius"
        ] == pytest.approx(1.0)
        if rest.client.get("/frame/current").json()["count"] > 1:
            rest.client.post("/frame/set", json={"frame": 1}).raise_for_status()
            assert rest.client.get("/api/learned-activity").json()["samples"] == []
            rest.client.post("/frame/set", json={"frame": 0}).raise_for_status()
            assert len(rest.client.get("/api/learned-activity").json()["samples"]) == 3
        rest.client.post(
            "/api/learned-activity", json={"visible": False}
        ).raise_for_status()
        assert rest.client.get("/api/learned-activity").json()["samples"] == []
        rest.client.post(
            "/api/run/load",
            json={"path": os.environ["H5READER_REST_FIXTURE"]},
            timeout=180,
        ).raise_for_status()
        assert not rest.client.get("/api/learned-activity").json()["loaded"]
    finally:
        rest.client.post("/api/learned-activity/clear").raise_for_status()


@pytest.mark.parametrize(
    "defect",
    [
        "wrong_position",
        "wrong_frame",
        "wrong_atoms",
        "nontraceless",
        "duplicate_atom",
        "wrong_basis",
    ],
)
def test_rejected_activity_preserves_previous_capture(rest, tmp_path, defect):
    good, _, _ = fixture_capture(rest)
    assert load_capture(rest, tmp_path, good).status_code == 200
    before = rest.client.get("/api/learned-activity").json()
    bad = copy.deepcopy(good)
    if defect == "wrong_position":
        bad["frames"][0]["positions"][0][0] += 1.0
    elif defect == "wrong_frame":
        bad["frames"][0]["frame_index"] = 10_000_000
    elif defect == "wrong_atoms":
        bad["atom_count"] += 1
    elif defect == "nontraceless":
        bad["frames"][0]["channels"][0]["tensors"][0][0] += 1
    elif defect == "duplicate_atom":
        bad["atom_indices"][1] = bad["atom_indices"][0]
    elif defect == "wrong_basis":
        bad["coordinate_frame"] = "model"
    try:
        response = load_capture(rest, tmp_path, bad)
        assert response.status_code == 400, response.text
        assert rest.client.get("/api/learned-activity").json() == before
        assert rest.client.get("/health").json()["ok"]
        assert (
            rest.client.post("/api/learned-activity", json={"radius": 0}).status_code
            == 400
        )
    finally:
        rest.client.post("/api/learned-activity/clear").raise_for_status()
