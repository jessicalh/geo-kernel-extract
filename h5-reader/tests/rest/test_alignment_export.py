"""End-to-end contract test for the target-free scientific alignment export."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import threading
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import h5py
import httpx
import numpy as np
import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("H5READER_ALIGNMENT_EXPORT_ENABLED") != "1",
    reason="alignment export integration uses its dedicated target-free fixture",
)


NUMERIC_CONTRACTS = {
    "aligned_positions.npy": (np.dtype("<f8"), 3),
    "rotations.npy": (np.dtype("<f8"), 3),
    "translations.npy": (np.dtype("<f8"), 2),
    "original_frame_index.npy": (np.dtype("<u8"), 1),
    "time_ps.npy": (np.dtype("<f8"), 1),
    "fit_atom_index.npy": (np.dtype("<u8"), 1),
    "fit_reference_positions.npy": (np.dtype("<f8"), 2),
    "fit_rmsd_A.npy": (np.dtype("<f8"), 1),
    "fit_singular_values.npy": (np.dtype("<f8"), 2),
    "fit_status.npy": (np.dtype("u1"), 1),
    "ca_rotations.npy": (np.dtype("<f8"), 3),
    "ca_translations.npy": (np.dtype("<f8"), 2),
    "ca_fit_atom_index.npy": (np.dtype("<u8"), 1),
    "ca_fit_reference_positions.npy": (np.dtype("<f8"), 2),
    "ca_fit_rmsd_A.npy": (np.dtype("<f8"), 1),
    "ca_fit_singular_values.npy": (np.dtype("<f8"), 2),
    "ca_fit_status.npy": (np.dtype("u1"), 1),
}

BYTE_COPIES = (
    "atoms_category_info.npy",
    "residues.npy",
    "bonds.npy",
    "rings.npy",
    "ring_membership.npy",
    "extraction_manifest.json",
    "native_mopac_complete.json",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while block := source.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _load_arrays(member_directory: Path) -> dict[str, np.ndarray]:
    arrays = {name: np.load(member_directory / name) for name in NUMERIC_CONTRACTS}
    for name, (dtype, rank) in NUMERIC_CONTRACTS.items():
        assert arrays[name].dtype == dtype, name
        assert arrays[name].ndim == rank, name
        assert np.isfinite(arrays[name]).all(), name
    return arrays


def _assert_reader_positions_match_export(
    rest, aligned: np.ndarray, frame: int, atoms: list[int]
) -> None:
    response = rest.client.post("/positions", json={"frame": frame, "atoms": atoms})
    assert response.status_code == 200, response.text
    displayed = np.asarray(
        [item["position"] for item in response.json()["positions"]],
        dtype=np.float64,
    )
    np.testing.assert_allclose(displayed, aligned[frame, atoms], atol=1e-12, rtol=0.0)


def _post_alignment_export(
    base_url: str,
    output_root: Path,
    apply_display: bool,
    started: threading.Event,
) -> httpx.Response:
    with httpx.Client(base_url=base_url, timeout=600.0) as client:
        started.set()
        return client.post(
            "/api/alignment/export",
            json={
                "output_root": str(output_root),
                "apply_display": apply_display,
            },
        )


def _wait_until_export_is_running(rest, future: Future[httpx.Response]) -> None:
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline:
        status = rest.client.get("/api/alignment/export/status")
        assert status.status_code == 200, status.text
        if status.json()["running"]:
            return
        if future.done():
            pytest.fail("alignment export completed before its running state was observable")
        time.sleep(0.01)
    pytest.fail("alignment export did not enter its running state")


def test_alignment_export_contract_and_resume(rest, tmp_path: Path) -> None:
    output_root = tmp_path / "alignment"
    output_root.mkdir()
    interrupted = output_root / ".bmr4095.staging.interrupted"
    interrupted.mkdir()
    (interrupted / "partial.npy").write_bytes(b"interrupted")

    started = threading.Event()
    with ThreadPoolExecutor(max_workers=1) as executor:
        export_future = executor.submit(
            _post_alignment_export,
            rest.base_url,
            output_root,
            True,
            started,
        )
        assert started.wait(timeout=2.0)
        _wait_until_export_is_running(rest, export_future)

        probe_started = time.monotonic()
        health = rest.client.get("/health", timeout=2.0)
        frame_state = rest.client.get("/frame/current", timeout=2.0)
        probe_elapsed = time.monotonic() - probe_started
        assert health.status_code == 200, health.text
        assert health.json()["ok"] is True
        assert frame_state.status_code == 200, frame_state.text
        assert probe_elapsed < 2.0

        response = export_future.result(timeout=600.0)
    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["ok"] is True
    assert payload["already_complete"] is False
    assert payload["display_applied"] is True
    member_directory = Path(payload["output_directory"])
    assert member_directory.is_dir()
    assert interrupted.is_dir()
    assert not (output_root / "COMPLETE.json").exists()

    manifest_path = member_directory / "export.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["complete"] is True
    assert manifest["contains_targets"] is False
    assert manifest["protected_targets_opened"] == 0
    assert manifest["transform_convention"] == (
        "aligned[t,a] = rotations[t] @ raw[t,a] + translations[t]"
    )
    assert set(manifest["files"]) == set(NUMERIC_CONTRACTS) | set(BYTE_COPIES)

    arrays = _load_arrays(member_directory)
    frame_count = manifest["dimensions"]["frames"]
    atom_count = manifest["dimensions"]["atoms"]
    primary_count = manifest["dimensions"]["primary_fit_atoms"]
    ca_count = manifest["dimensions"]["ca_fit_atoms"]
    assert frame_count == 100
    assert manifest["dimensions"]["residues"] > 0
    assert manifest["dimensions"]["bonds"] >= 0
    assert manifest["dimensions"]["rings"] == (
        manifest["dimensions"]["aromatic_rings"]
        + manifest["dimensions"]["saturated_rings"]
    )
    assert manifest["dimensions"]["ring_memberships"] >= 0
    assert manifest["source"]["provenance_authority"] == (
        "READER_ALIGNMENT_EXPORT_CONTRACT.md"
    )
    assert arrays["aligned_positions.npy"].shape == (frame_count, atom_count, 3)
    assert arrays["rotations.npy"].shape == (frame_count, 3, 3)
    assert arrays["translations.npy"].shape == (frame_count, 3)
    assert arrays["fit_atom_index.npy"].shape == (primary_count,)
    assert arrays["ca_fit_atom_index.npy"].shape == (ca_count,)
    assert np.all(arrays["fit_status.npy"] == 0)
    assert np.all(arrays["ca_fit_status.npy"] == 0)

    source = manifest["source"]
    extraction_directory = Path(source["extraction_directory"])
    h5_path = Path(source["trajectory_h5"])
    with h5py.File(h5_path, "r") as trajectory:
        raw_positions = np.asarray(
            trajectory["/trajectory/positions/xyz"], dtype=np.float64
        ).transpose(1, 0, 2)
        h5_frames = np.asarray(
            trajectory["/trajectory/frames/original_index"], dtype=np.uint64
        )
        h5_times = np.asarray(
            trajectory["/trajectory/frames/time_ps"], dtype=np.float64
        )
        h5_fit_indices = np.asarray(
            trajectory["/trajectory/rmsd_tracking/atom_indices"], dtype=np.uint64
        )

    np.testing.assert_array_equal(arrays["original_frame_index.npy"], h5_frames)
    np.testing.assert_array_equal(arrays["time_ps.npy"], h5_times)
    np.testing.assert_array_equal(arrays["fit_atom_index.npy"], h5_fit_indices)
    np.testing.assert_array_equal(h5_frames, np.arange(15, 1501, 15, dtype=np.uint64))
    np.testing.assert_allclose(
        h5_times, np.arange(150.0, 15001.0, 150.0), atol=1e-9, rtol=0.0
    )
    for row, original_frame in enumerate(h5_frames):
        source_positions = np.load(
            extraction_directory
            / "npys"
            / f"frame_{int(original_frame):06d}"
            / "pos.npy"
        )
        np.testing.assert_allclose(
            raw_positions[row], source_positions, atol=1e-10, rtol=1e-12
        )

    rotations = arrays["rotations.npy"]
    translations = arrays["translations.npy"]
    aligned = arrays["aligned_positions.npy"]
    reconstructed = np.einsum("tij,taj->tai", rotations, raw_positions)
    reconstructed += translations[:, None, :]
    np.testing.assert_allclose(aligned, reconstructed, atol=1e-12, rtol=0.0)
    inverse = np.einsum("tji,taj->tai", rotations, aligned - translations[:, None, :])
    np.testing.assert_allclose(inverse, raw_positions, atol=1e-10, rtol=0.0)

    identity = np.eye(3)[None, :, :]
    orthogonality = np.matmul(rotations, np.transpose(rotations, (0, 2, 1)))
    assert np.max(np.linalg.norm(orthogonality - identity, axis=(1, 2))) <= 1e-8
    assert np.max(np.abs(np.linalg.det(rotations) - 1.0)) <= 1e-8
    ca_rotations = arrays["ca_rotations.npy"]
    ca_orthogonality = np.matmul(ca_rotations, np.transpose(ca_rotations, (0, 2, 1)))
    assert np.max(np.linalg.norm(ca_orthogonality - identity, axis=(1, 2))) <= 1e-8
    assert np.max(np.abs(np.linalg.det(ca_rotations) - 1.0)) <= 1e-8
    ca_translations = arrays["ca_translations.npy"]
    ca_aligned = np.einsum("tij,taj->tai", ca_rotations, raw_positions)
    ca_aligned += ca_translations[:, None, :]
    ca_round_trip = np.einsum(
        "tji,taj->tai",
        ca_rotations,
        ca_aligned - ca_translations[:, None, :],
    )
    np.testing.assert_allclose(ca_round_trip, raw_positions, atol=1e-10, rtol=0.0)

    validation = manifest["validation"]
    assert validation["ca_max_rotation_orthogonality_frobenius_error"] <= 1e-8
    assert abs(validation["ca_determinant_min"] - 1.0) <= 1e-8
    assert abs(validation["ca_determinant_max"] - 1.0) <= 1e-8
    assert validation["ca_max_round_trip_reconstruction_error_A"] <= 1e-10

    for name in BYTE_COPIES:
        source_path = extraction_directory / name
        if name == "extraction_manifest.json":
            source_path = Path(
                next(
                    item["path"]
                    for item in source["input_files"]
                    if item["role"] == "copied_source"
                    and Path(item["path"]).name == name
                )
            )
        assert _sha256(member_directory / name) == _sha256(source_path)

    table = (output_root / "members.tsv").read_text(encoding="utf-8")
    assert payload["member_id"] in table
    assert ".staging." not in table

    state = rest.client.get("/transform").json()
    assert state["kind"] == "scientific_alignment"
    assert state["window"] == 0
    assert state["subset_atoms"] == arrays["fit_atom_index.npy"].tolist()
    display_atoms = sorted({0, atom_count // 2, atom_count - 1})
    for frame in (0, frame_count // 2, frame_count - 1):
        _assert_reader_positions_match_export(
            rest, arrays["aligned_positions.npy"], frame, display_atoms
        )

    smoothing = rest.client.post("/transform/smoothing", json={"window": 1})
    assert smoothing.status_code == 409
    assert rest.client.get("/transform").json()["kind"] == "scientific_alignment"

    ordinary_display = rest.client.post(
        "/transform", json={"kind": "all_atom_fit", "reference_frame": 0}
    )
    assert ordinary_display.status_code == 204, ordinary_display.text
    assert rest.client.get("/transform").json()["kind"] == "all_atom_fit"

    for _ in range(3):
        repeated = rest.client.post(
            "/api/alignment/export",
            json={"output_root": str(output_root), "apply_display": True},
            timeout=600.0,
        )
        assert repeated.status_code == 200, repeated.text
        repeated_payload = repeated.json()
        assert repeated_payload["ok"] is True
        assert repeated_payload["already_complete"] is True
        assert repeated_payload["display_applied"] is True
        assert Path(repeated_payload["output_directory"]) == member_directory
    assert rest.client.get("/transform").json()["kind"] == "scientific_alignment"
    _assert_reader_positions_match_export(
        rest, arrays["aligned_positions.npy"], frame_count - 1, display_atoms
    )

    overlap = rest.client.post(
        "/api/alignment/export",
        json={"output_root": source["calcset_root"], "apply_display": False},
        timeout=600.0,
    )
    assert overlap.status_code == 409
    assert overlap.json()["ok"] is False

    invalid_root = tmp_path / "invalid-final"
    invalid_member = invalid_root / payload["member_id"]
    invalid_member.mkdir(parents=True)
    (invalid_member / "not-an-export.txt").write_text(
        "do not replace", encoding="ascii"
    )
    invalid = rest.client.post(
        "/api/alignment/export",
        json={"output_root": str(invalid_root), "apply_display": False},
        timeout=600.0,
    )
    assert invalid.status_code == 409
    assert (invalid_member / "not-an-export.txt").read_text(
        encoding="ascii"
    ) == "do not replace"

    shutil.rmtree(interrupted)


def test_alignment_export_cancellation_cleans_staging(rest, tmp_path: Path) -> None:
    output_root = tmp_path / "cancelled-alignment"
    output_root.mkdir()
    started = threading.Event()

    with ThreadPoolExecutor(max_workers=1) as executor:
        export_future = executor.submit(
            _post_alignment_export,
            rest.base_url,
            output_root,
            False,
            started,
        )
        assert started.wait(timeout=2.0)
        _wait_until_export_is_running(rest, export_future)

        cancelled = rest.client.post("/api/alignment/export/cancel")
        assert cancelled.status_code == 202, cancelled.text
        assert cancelled.json()["cancel_requested"] is True

        response = export_future.result(timeout=30.0)

    assert response.status_code == 409, response.text
    payload = response.json()
    assert payload["ok"] is False
    assert payload["cancelled"] is True
    assert not any(path.name.startswith(".") for path in output_root.iterdir())
    assert not any(path.is_dir() for path in output_root.iterdir())
    status = rest.client.get("/api/alignment/export/status")
    assert status.status_code == 200
    assert status.json()["running"] is False
