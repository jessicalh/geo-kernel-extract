"""REST loading uses the same run path as the Reader UI."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest


OUTPUT_FILES = {
    "positions.npy",
    "scalars.npy",
    "scalar_valid.npy",
    "applicability.npy",
    "l1.npy",
    "rotations.npy",
    "translations.npy",
    "fit_status.npy",
    "fit_rmsd_A.npy",
    "atoms.csv",
    "bonds.csv",
    "frames.csv",
}


def _load_path_from_environment(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        pytest.skip(f"{name} is not configured")
    path = Path(value).resolve()
    if not path.exists():
        pytest.fail(f"{name} does not exist: {path}")
    return path


def _lgs_file(path: Path) -> Path:
    if path.is_file():
        return path
    candidates = sorted(
        child for child in path.iterdir()
        if child.is_file() and child.suffix.lower() == ".lgs"
    )
    assert len(candidates) == 1, f"expected one .LGS in {path}"
    return candidates[0]


def _position(rest, atom: int = 0, frame: int = 0):
    response = rest.client.post("/positions", json={"atoms": [atom], "frame": frame})
    assert response.status_code == 200, response.text
    return response.json()["positions"][0]["position"]


def test_run_load_switches_real_runs_and_exports_the_new_run(rest, tmp_path):
    original = _load_path_from_environment("H5READER_REST_FIXTURE")
    replacement = _load_path_from_environment("H5READER_REST_RELOAD_FIXTURE")
    before = rest.client.get("/ui/state").json()
    before_atoms = rest.client.get("/protein/atoms").json()["count"]

    try:
        response = rest.client.post(
            "/api/run/load", json={"path": str(replacement)}, timeout=900.0
        )
        assert response.status_code == 200, response.text
        assert response.json() == {"ok": True}

        loaded = rest.client.get("/ui/state").json()
        loaded_atoms = rest.client.get("/protein/atoms").json()["count"]
        assert loaded["loaded"] is True
        if _lgs_file(replacement).resolve() != _lgs_file(original).resolve():
            assert (loaded["protein"], loaded["frames"], loaded_atoms) != (
                before["protein"],
                before["frames"],
                before_atoms,
            )

        output = tmp_path / "replacement-export"
        output.mkdir()
        response = rest.client.post(
            "/api/model-input/export",
            json={"output_directory": str(output)},
            timeout=900.0,
        )
        assert response.status_code == 200, response.text
        assert response.json()["frames"] == loaded["frames"]
        assert response.json()["atoms"] == loaded_atoms
        assert {path.name for path in output.iterdir()} == OUTPUT_FILES
        positions = np.load(output / "positions.npy")
        assert positions.shape == (
            loaded["frames"],
            loaded_atoms,
            3,
        )
        np.testing.assert_allclose(
            positions[0, 0] - positions[0, 1],
            np.asarray(_position(rest, 0)) - np.asarray(_position(rest, 1)),
            rtol=0.0,
            atol=2.0e-5,
        )
    finally:
        restore_path = original.parent if original.is_file() else original
        if len(list(restore_path.glob("*.LGS"))) != 1:
            restore_path = original
        response = rest.client.post("/api/run/load", json={"path": str(restore_path)})
        assert response.status_code == 200, response.text

    restored = rest.client.get("/ui/state").json()
    assert restored["protein"] == before["protein"]
    assert restored["frames"] == before["frames"]
    assert rest.client.get("/protein/atoms").json()["count"] == before_atoms
    assert rest.client.get("/health").status_code == 200


def test_run_load_rejects_bad_requests_without_replacing_the_run(rest, tmp_path):
    before = rest.client.get("/ui/state").json()
    before_atoms = rest.client.get("/protein/atoms").json()["count"]
    before_position = _position(rest)

    response = rest.client.post("/api/run/load", json={"path": "", "dataset": "none"})
    assert response.status_code == 400

    missing = tmp_path / "missing.LGS"
    response = rest.client.post("/api/run/load", json={"path": str(missing)})
    assert response.status_code == 409
    assert response.json()["error"]

    malformed = tmp_path / "malformed.LGS"
    malformed.write_text("{", encoding="utf-8")
    response = rest.client.post("/api/run/load", json={"path": str(malformed)})
    assert response.status_code == 409
    assert response.json()["error"]

    ambiguous = tmp_path / "ambiguous"
    ambiguous.mkdir()
    (ambiguous / "first.LGS").write_text("{}", encoding="utf-8")
    (ambiguous / "second.LGS").write_text("{}", encoding="utf-8")
    response = rest.client.post("/api/run/load", json={"path": str(ambiguous)})
    assert response.status_code == 409
    assert response.json()["error"]

    source_lgs = _lgs_file(_load_path_from_environment("H5READER_REST_FIXTURE"))
    source_manifest = json.loads(source_lgs.read_text(encoding="utf-8"))
    if source_manifest.get("kind") == "trajectory":
        trajectory = source_manifest["trajectory"]
        for key in ("md_dir", "topology_top", "extraction_dir", "extraction_manifest"):
            value = Path(trajectory[key])
            if not value.is_absolute():
                value = (source_lgs.parent / value).resolve()
            trajectory[key] = str(value)
        corrupt_h5 = tmp_path / "corrupt.h5"
        corrupt_h5.write_bytes(b"not an HDF5 file")
        trajectory["trajectory_h5"] = str(corrupt_h5)
        source_manifest.pop("dft", None)
        late_failure = tmp_path / "late-failure.LGS"
        late_failure.write_text(json.dumps(source_manifest), encoding="utf-8")
        response = rest.client.post(
            "/api/run/load", json={"path": str(late_failure)}
        )
        assert response.status_code == 409
        assert response.json()["error"]

    after = rest.client.get("/ui/state").json()
    assert after["loaded"] is True
    assert after["protein"] == before["protein"]
    assert after["frames"] == before["frames"]
    assert rest.client.get("/protein/atoms").json()["count"] == before_atoms
    assert _position(rest) == before_position
