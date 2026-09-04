"""Reader-side contract test for the structural model-input export."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("H5READER_MODEL_INPUT_EXPORT_ENABLED") != "1",
    reason="model-input export uses an explicitly configured complete fixture",
)


FILES = {
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

ATOM_COLUMNS = [
    "atom_index",
    "chain_id",
    "residue_number",
    "insertion_code",
    "residue_index",
    "previous_residue_index",
    "next_residue_index",
    "amber_atom_name",
    "bmrb_atom_name",
    "parent_atom_index",
    "element",
    "iupac_atom",
    "iupac_residue",
    "amber_residue",
    "terminal_state",
    "residue_variant",
    "formal_charge",
    "polar_h_kind",
    "atom_role",
    "hybridisation",
    "backbone_role",
    "planar_group",
    "locant",
    "ring_position",
    "aromatic",
    "exchangeable",
]


def _csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        return list(reader.fieldnames or []), list(reader)


def _decoded(values: np.ndarray) -> list[str]:
    return [
        value.decode("ascii").rstrip("\x00 ")
        if isinstance(value, bytes)
        else str(value).rstrip("\x00 ")
        for value in values
    ]


def _assert_atoms_match_sources(
    atoms: list[dict[str, str]],
    source_atoms: np.ndarray,
    source_residues: np.ndarray,
    frame_directory: Path,
) -> None:
    text_columns = {
        "chain_id": "chain_id",
        "insertion_code": "insertion_code",
        "amber_atom_name": "amber_atom_name",
        "bmrb_atom_name": "bmrb_atom_name",
        "iupac_atom": "iupac_atom_name",
        "iupac_residue": "iupac_residue_3letter",
        "amber_residue": "amber_residue_3letter",
    }
    integer_columns = {
        "atom_index": "atom_index",
        "residue_number": "residue_number",
        "residue_index": "residue_index",
        "parent_atom_index": "parent_atom_index",
        "element": "element",
        "terminal_state": "terminal_state",
        "residue_variant": "residue_variant_index",
        "formal_charge": "formal_charge",
        "polar_h_kind": "polar_h_kind",
        "backbone_role": "backbone_role",
        "planar_group": "planar_group",
        "locant": "locant",
        "ring_position": "ring_position_primary",
        "aromatic": "aromatic",
        "exchangeable": "is_exchangeable",
    }

    for output_name, source_name in text_columns.items():
        assert [row[output_name] for row in atoms] == _decoded(
            source_atoms[source_name]
        )
    for output_name, source_name in integer_columns.items():
        np.testing.assert_array_equal(
            [int(row[output_name]) for row in atoms], source_atoms[source_name]
        )

    residues_by_index = {
        int(residue["residue_index"]): residue for residue in source_residues
    }
    for output_name, source_name in (
        ("previous_residue_index", "prev_residue_index"),
        ("next_residue_index", "next_residue_index"),
    ):
        expected = [
            int(residues_by_index[int(atom["residue_index"])][source_name])
            for atom in source_atoms
        ]
        assert [int(row[output_name]) for row in atoms] == expected

    np.testing.assert_array_equal(
        [int(row["atom_role"]) for row in atoms],
        np.load(frame_directory / "enrichment_role.npy"),
    )
    np.testing.assert_array_equal(
        [int(row["hybridisation"]) for row in atoms],
        np.load(frame_directory / "enrichment_hybridisation_class.npy"),
    )


def _source_directories(frame_id: str) -> tuple[Path, Path]:
    lgs_path = Path(os.environ["H5READER_REST_FIXTURE"])
    lgs = json.loads(lgs_path.read_text(encoding="utf-8"))
    if lgs["kind"] == "trajectory":
        topology = (lgs_path.parent / lgs["trajectory"]["extraction_dir"]).resolve()
        return topology, topology / "npys" / frame_id
    topology = (lgs_path.parent / lgs["single_pose"]["pose_dir"]).resolve()
    return topology, topology


def _expected_scalar_arrays(
    frame_directory: Path,
    residue_index: np.ndarray,
    formal_charge: int,
) -> tuple[np.ndarray, np.ndarray]:
    atom_count = len(residue_index)
    values = np.zeros((atom_count, 35), dtype=np.float32)
    valid = np.zeros((atom_count, 35), dtype=np.uint8)
    channel = 0

    def masked(
        source: np.ndarray, mask: np.ndarray, *, residue_axis: bool = False
    ) -> None:
        nonlocal channel
        atom_mask = np.asarray(mask[residue_index] if residue_axis else mask)
        atom_values = np.asarray(source[residue_index] if residue_axis else source)
        valid[:, channel] = atom_mask.astype(np.uint8)
        values[:, channel] = np.where(atom_mask == 1, atom_values, 0.0).astype(
            np.float32
        )
        channel += 1

    for stem in (
        "tau_N_CA_C",
        "angle_N_CA_CB",
        "angle_CB_CA_C",
        "angle_Cprev_N_CA",
        "angle_CA_C_Nnext",
    ):
        masked(
            np.load(frame_directory / f"{stem}.npy"),
            np.load(frame_directory / f"{stem}_valid.npy"),
            residue_axis=True,
        )

    torsion_valid = np.load(frame_directory / "dssp_torsion_valid.npy")
    for stem in ("dssp_torsion_cos", "dssp_torsion_sin"):
        torsion = np.load(frame_directory / f"{stem}.npy")
        for component in range(6):
            masked(
                torsion[:, component],
                torsion_valid[:, component],
                residue_axis=True,
            )

    omega_valid = np.load(frame_directory / "omega_valid.npy")
    masked(
        np.load(frame_directory / "omega_cos.npy"),
        omega_valid,
        residue_axis=True,
    )
    masked(
        np.load(frame_directory / "omega_sin.npy"),
        omega_valid,
        residue_axis=True,
    )

    values[:, channel] = np.load(frame_directory / "omega_is_xpro.npy")[
        residue_index
    ].astype(np.float32)
    valid[:, channel] = 1
    channel += 1

    masked(
        np.load(frame_directory / "pyramidalization.npy"),
        np.load(frame_directory / "pyramidalization_valid.npy"),
    )

    secondary_structure = np.load(frame_directory / "dssp_ss8.npy")
    values[:, channel : channel + 8] = secondary_structure.astype(np.float32)
    valid[:, channel : channel + 8] = 1
    channel += 8

    dssp_observed = np.load(frame_directory / "dssp_observed.npy")
    hbond_energy = np.load(frame_directory / "dssp_hbond_energy.npy")
    for component in range(4):
        masked(hbond_energy[:, component], dssp_observed)
    masked(np.load(frame_directory / "dssp_ppii.npy"), dssp_observed)

    values[:, channel] = np.float32(formal_charge)
    valid[:, channel] = 1
    channel += 1
    assert channel == 35
    return values, valid


def _expected_cb_residuals(
    source_positions: np.ndarray,
    rotation: np.ndarray,
    atom_rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    atom_count = len(atom_rows)
    vectors = np.zeros((atom_count, 1, 3), dtype=np.float32)
    applicable = np.zeros(atom_count, dtype=np.float32)
    residue_indices = atom_rows["residue_index"]
    atom_names = np.char.decode(atom_rows["iupac_atom_name"], "ascii")

    for residue in np.unique(residue_indices):
        residue_atoms = np.flatnonzero(residue_indices == residue)
        named = {atom_names[index].rstrip("\x00 "): index for index in residue_atoms}
        if not {"N", "CA", "C", "CB"}.issubset(named):
            continue
        n = source_positions[named["N"]] - source_positions[named["CA"]]
        c = source_positions[named["C"]] - source_positions[named["CA"]]
        ideal_cb = (
            source_positions[named["CA"]]
            - 0.58273431 * np.cross(c, n)
            - 0.56802827 * n
            - 0.54067466 * c
        )
        residual = rotation @ (source_positions[named["CB"]] - ideal_cb)
        vectors[residue_atoms, 0, :] = residual.astype(np.float32)
        applicable[residue_atoms] = 1.0
    return vectors, applicable


def test_model_input_export_matches_reader_state(rest, tmp_path: Path) -> None:
    extra_option = rest.client.post(
        "/api/model-input/export",
        json={"output_directory": str(tmp_path), "extra": True},
    )
    assert extra_option.status_code == 400

    occupied = tmp_path / "occupied"
    occupied.mkdir()
    marker = occupied / "keep.txt"
    marker.write_text("unchanged", encoding="utf-8")
    refused = rest.client.post(
        "/api/model-input/export",
        json={"output_directory": str(occupied)},
    )
    assert refused.status_code == 409
    assert set(occupied.iterdir()) == {marker}
    assert marker.read_text(encoding="utf-8") == "unchanged"

    output = tmp_path / "export"
    output.mkdir()
    response = rest.client.post(
        "/api/model-input/export",
        json={"output_directory": str(output)},
        timeout=900.0,
    )
    assert response.status_code == 200, response.text
    counts = response.json()
    assert set(counts) == {"frames", "atoms", "bonds"}
    frame_count = counts["frames"]
    atom_count = counts["atoms"]
    bond_count = counts["bonds"]
    assert frame_count > 0
    assert atom_count > 0
    assert bond_count >= 0
    assert {path.name for path in output.iterdir()} == FILES

    positions = np.load(output / "positions.npy")
    scalars = np.load(output / "scalars.npy")
    scalar_valid = np.load(output / "scalar_valid.npy")
    applicability = np.load(output / "applicability.npy")
    l1 = np.load(output / "l1.npy")
    rotations = np.load(output / "rotations.npy")
    translations = np.load(output / "translations.npy")
    fit_status = np.load(output / "fit_status.npy")
    fit_rmsd = np.load(output / "fit_rmsd_A.npy")

    assert positions.dtype == np.dtype("<f4")
    assert scalars.dtype == np.dtype("<f4")
    assert scalar_valid.dtype == np.dtype("u1")
    assert applicability.dtype == np.dtype("<f4")
    assert l1.dtype == np.dtype("<f4")
    assert rotations.dtype == np.dtype("<f8")
    assert translations.dtype == np.dtype("<f8")
    assert fit_status.dtype == np.dtype("u1")
    assert fit_rmsd.dtype == np.dtype("<f8")

    assert positions.shape == (frame_count, atom_count, 3)
    assert scalars.shape == (frame_count, atom_count, 35)
    assert scalar_valid.shape == scalars.shape
    assert applicability.shape == (frame_count, atom_count, 8)
    assert l1.shape == (frame_count, atom_count, 1, 3)
    assert rotations.shape == (frame_count, 3, 3)
    assert translations.shape == (frame_count, 3)
    assert fit_status.shape == (frame_count,)
    assert fit_rmsd.shape == (frame_count,)

    assert np.isfinite(positions).all()
    assert np.isfinite(scalars).all()
    assert np.isfinite(applicability).all()
    assert np.isfinite(l1).all()
    assert np.isfinite(rotations).all()
    assert np.isfinite(translations).all()
    assert np.isfinite(fit_rmsd).all()
    assert set(np.unique(scalar_valid)).issubset({0, 1})
    assert set(np.unique(applicability)).issubset({0.0, 1.0})
    assert np.all(scalars[scalar_valid == 0] == 0.0)
    assert np.all(l1[applicability[..., 0] == 0] == 0.0)
    np.testing.assert_array_equal(applicability[..., 1:7], scalar_valid[..., 5:11])
    np.testing.assert_array_equal(applicability[..., 7], scalar_valid[..., 20])
    assert np.all(scalar_valid[..., 19] == 1)
    assert np.all(scalar_valid[..., 21:29] == 1)
    assert np.all(scalar_valid[..., 34] == 1)
    assert np.all(fit_status == 0)
    assert np.all(fit_rmsd >= 0.0)
    np.testing.assert_allclose(
        positions.mean(axis=1, dtype=np.float64), 0.0, atol=2e-5, rtol=0.0
    )
    for rotation in rotations:
        np.testing.assert_allclose(
            rotation @ rotation.T, np.eye(3), atol=1e-8, rtol=0.0
        )
        assert np.linalg.det(rotation) > 0.0

    atom_header, atoms = _csv_rows(output / "atoms.csv")
    assert atom_header == ATOM_COLUMNS
    assert len(atoms) == atom_count
    assert [int(row["atom_index"]) for row in atoms] == list(range(atom_count))

    bond_header, bonds = _csv_rows(output / "bonds.csv")
    assert bond_header == ["bond_index", "atom_a", "atom_b", "order", "category"]
    assert len(bonds) == bond_count
    assert [int(row["bond_index"]) for row in bonds] == list(range(bond_count))
    assert all(0 <= int(row["atom_a"]) < atom_count for row in bonds)
    assert all(0 <= int(row["atom_b"]) < atom_count for row in bonds)

    frame_header, frames = _csv_rows(output / "frames.csv")
    assert frame_header == [
        "frame_row",
        "frame_id",
        "original_frame_index",
        "time_ps",
    ]
    assert len(frames) == frame_count
    assert [int(row["frame_row"]) for row in frames] == list(range(frame_count))

    trajectory = frames[0]["frame_id"] != "static"
    if trajectory:
        for row in frames:
            source_frame = int(row["original_frame_index"])
            assert row["frame_id"] == f"frame_{source_frame:06d}"
            assert np.isfinite(float(row["time_ps"]))
        transform = rest.client.get("/transform")
        assert transform.status_code == 200
        assert transform.json()["kind"] == "scientific_alignment"
    else:
        assert frame_count == 1
        assert frames[0]["original_frame_index"] == ""
        assert frames[0]["time_ps"] == ""
        np.testing.assert_array_equal(rotations[0], np.eye(3))
        np.testing.assert_array_equal(translations[0], np.zeros(3))
        assert fit_rmsd[0] == 0.0

    all_atoms = list(range(atom_count))
    sample_frames = sorted({0, frame_count // 2, frame_count - 1})
    for frame in sample_frames:
        topology_directory, frame_directory = _source_directories(
            frames[frame]["frame_id"]
        )
        source_atoms = np.load(topology_directory / "atoms_category_info.npy")
        source_positions = np.load(frame_directory / "pos.npy").astype(np.float64)
        residue_index = source_atoms["residue_index"].astype(np.int64)
        formal_charge = int(source_atoms["formal_charge"].sum(dtype=np.int64))

        if frame == sample_frames[0]:
            source_residues = np.load(topology_directory / "residues.npy")
            _assert_atoms_match_sources(
                atoms, source_atoms, source_residues, frame_directory
            )
            source_bonds = np.load(topology_directory / "bonds.npy")
            np.testing.assert_array_equal(
                [int(row["atom_a"]) for row in bonds],
                source_bonds["atom_index_a"],
            )
            np.testing.assert_array_equal(
                [int(row["atom_b"]) for row in bonds],
                source_bonds["atom_index_b"],
            )
            np.testing.assert_array_equal(
                [int(row["order"]) for row in bonds],
                source_bonds["bond_order"],
            )
            np.testing.assert_array_equal(
                [int(row["category"]) for row in bonds],
                source_bonds["bond_category"],
            )

        expected_scalars, expected_valid = _expected_scalar_arrays(
            frame_directory, residue_index, formal_charge
        )
        np.testing.assert_array_equal(scalars[frame], expected_scalars)
        np.testing.assert_array_equal(scalar_valid[frame], expected_valid)

        aligned = source_positions @ rotations[frame].T + translations[frame]
        aligned -= aligned.mean(axis=0)
        np.testing.assert_allclose(
            positions[frame], aligned.astype(np.float32), atol=2e-5, rtol=0.0
        )

        expected_l1, expected_l1_valid = _expected_cb_residuals(
            source_positions, rotations[frame], source_atoms
        )
        np.testing.assert_allclose(l1[frame], expected_l1, atol=2e-5, rtol=0.0)
        np.testing.assert_array_equal(applicability[frame, :, 0], expected_l1_valid)

        live = rest.client.post("/positions", json={"frame": frame, "atoms": all_atoms})
        assert live.status_code == 200, live.text
        displayed = np.asarray(
            [entry["position"] for entry in live.json()["positions"]],
            dtype=np.float64,
        )
        displayed -= displayed.mean(axis=0)
        np.testing.assert_allclose(positions[frame], displayed, atol=2e-5, rtol=0.0)
