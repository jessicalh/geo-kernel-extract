#!/usr/bin/env python3
"""NPY boundary checks for the local-geometry producer additions.

Invoked by the C++ production-path forcing test.  Every array is opened with
allow_pickle=False; coordinate oracles below use only emitted coordinates and
the invariant topology sidecars.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


LOCAL_SCALARS = (
    "tau_N_CA_C",
    "angle_N_CA_CB",
    "angle_CB_CA_C",
    "angle_Cprev_N_CA",
    "angle_CA_C_Nnext",
    "cb_deviation",
)


def load(path: Path) -> np.ndarray:
    return np.load(path, allow_pickle=False)


def backbone_separation(residues: np.ndarray, start: int, goal: int) -> int:
    if start == goal:
        return 0
    distance = [-1] * len(residues)
    distance[start] = 0
    queue = [start]
    for current in queue:
        for field in ("prev_residue_index", "next_residue_index"):
            neighbour = int(residues[field][current])
            if neighbour < 0 or distance[neighbour] >= 0:
                continue
            distance[neighbour] = distance[current] + 1
            if neighbour == goal:
                return distance[neighbour]
            queue.append(neighbour)
    return -1


def check_output(output_dir: Path) -> None:
    sidecar_dir = output_dir if (output_dir / "bonds.npy").is_file() else output_dir.parent

    pos = load(output_dir / "pos.npy")
    residue_index = load(output_dir / "residue_index.npy")
    bonds = load(sidecar_dir / "bonds.npy")
    residues = load(sidecar_dir / "residues.npy")
    n_atoms = pos.shape[0]
    n_residues = residues.shape[0]
    n_bonds = bonds.shape[0]

    assert pos.shape == (n_atoms, 3) and pos.dtype == np.float64

    bond_length = load(output_dir / "bond_length.npy")
    bond_direction = load(output_dir / "bond_direction.npy")
    bond_valid = load(output_dir / "bond_geometry_valid.npy")
    assert bond_length.shape == (n_bonds,) and bond_length.dtype == np.float64
    assert bond_direction.shape == (n_bonds, 3) and bond_direction.dtype == np.float64
    assert bond_valid.shape == (n_bonds,) and bond_valid.dtype == np.uint8
    np.testing.assert_array_equal(bonds["bond_index"], np.arange(n_bonds, dtype=np.int32))
    delta = pos[bonds["atom_index_b"]] - pos[bonds["atom_index_a"]]
    direct_length = np.linalg.norm(delta, axis=1)
    direct_valid = np.isfinite(delta).all(axis=1) & (direct_length > 1.0e-15)
    direct_direction = np.zeros_like(delta)
    direct_direction[direct_valid] = delta[direct_valid] / direct_length[direct_valid, None]
    np.testing.assert_allclose(bond_length, direct_length, rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(bond_direction, direct_direction, rtol=0.0, atol=1.0e-12)
    np.testing.assert_array_equal(bond_valid, direct_valid.astype(np.uint8))

    for stem in LOCAL_SCALARS:
        value = load(output_dir / f"{stem}.npy")
        valid = load(output_dir / f"{stem}_valid.npy")
        assert value.shape == (n_residues,) and value.dtype == np.float64
        assert valid.shape == (n_residues,) and valid.dtype == np.uint8
        np.testing.assert_array_equal(valid, np.isfinite(value).astype(np.uint8))
    cb_vector = load(output_dir / "cb_residual_vector.npy")
    cb_vector_valid = load(output_dir / "cb_residual_vector_valid.npy")
    assert cb_vector.shape == (n_residues, 3) and cb_vector.dtype == np.float64
    assert cb_vector_valid.shape == (n_residues,) and cb_vector_valid.dtype == np.uint8
    np.testing.assert_array_equal(
        cb_vector_valid, np.isfinite(cb_vector).all(axis=1).astype(np.uint8)
    )
    cb_deviation = load(output_dir / "cb_deviation.npy")
    valid_cb = cb_vector_valid.astype(bool)
    np.testing.assert_allclose(
        np.linalg.norm(cb_vector[valid_cb], axis=1), cb_deviation[valid_cb],
        rtol=0.0, atol=1.0e-12,
    )

    torsion_angle = load(output_dir / "dssp_torsion_angle.npy")
    torsion_sin = load(output_dir / "dssp_torsion_sin.npy")
    torsion_cos = load(output_dir / "dssp_torsion_cos.npy")
    torsion_valid = load(output_dir / "dssp_torsion_valid.npy")
    for array in (torsion_angle, torsion_sin, torsion_cos):
        assert array.shape == (n_residues, 6) and array.dtype == np.float64
    assert torsion_valid.shape == (n_residues, 6) and torsion_valid.dtype == np.uint8
    valid_torsion = torsion_valid.astype(bool)
    np.testing.assert_array_equal(valid_torsion, np.isfinite(torsion_angle))
    np.testing.assert_array_equal(valid_torsion, np.isfinite(torsion_sin))
    np.testing.assert_array_equal(valid_torsion, np.isfinite(torsion_cos))
    np.testing.assert_allclose(
        torsion_sin[valid_torsion], np.sin(torsion_angle[valid_torsion]),
        rtol=0.0, atol=1.0e-12,
    )
    np.testing.assert_allclose(
        torsion_cos[valid_torsion], np.cos(torsion_angle[valid_torsion]),
        rtol=0.0, atol=1.0e-12,
    )

    omega = load(output_dir / "omega_actual.npy")
    omega_sin = load(output_dir / "omega_sin.npy")
    omega_cos = load(output_dir / "omega_cos.npy")
    omega_valid = load(output_dir / "omega_valid.npy")
    assert omega.shape == omega_sin.shape == omega_cos.shape == omega_valid.shape == (n_residues,)
    assert omega.dtype == omega_sin.dtype == omega_cos.dtype == np.float64
    assert omega_valid.dtype == np.uint8
    valid_omega = omega_valid.astype(bool)
    np.testing.assert_array_equal(valid_omega, np.isfinite(omega))
    np.testing.assert_array_equal(valid_omega, np.isfinite(omega_sin))
    np.testing.assert_array_equal(valid_omega, np.isfinite(omega_cos))
    np.testing.assert_allclose(omega_sin[valid_omega], np.sin(omega[valid_omega]), atol=1.0e-12)
    np.testing.assert_allclose(omega_cos[valid_omega], np.cos(omega[valid_omega]), atol=1.0e-12)

    partner = load(output_dir / "dssp_hbond_partner_residue_index.npy")
    energy = load(output_dir / "dssp_hbond_energy.npy")
    observed = load(output_dir / "dssp_observed.npy").astype(bool)
    assert partner.shape == energy.shape == (n_atoms, 4)
    assert partner.dtype == np.int32 and energy.dtype == np.float64
    assert np.all((partner == -1) | ((partner >= 0) & (partner < n_residues)))
    missing_partner = observed[:, None] & (partner == -1)
    assert np.all(energy[missing_partner] == 0.0)
    for ri in range(n_residues):
        atom_rows = np.flatnonzero(residue_index == ri)
        if atom_rows.size:
            assert np.all(partner[atom_rows] == partner[atom_rows[0]])

    pair_index = load(output_dir / "hbond_pairs_index.npy")
    pair_geometry = load(output_dir / "hbond_pairs_geometry.npy")
    pair_angle_valid = load(output_dir / "hbond_pairs_angle_valid.npy")
    assert pair_index.ndim == 2 and pair_index.shape[1] == 6 and pair_index.dtype == np.int32
    assert pair_geometry.shape == (pair_index.shape[0], 5) and pair_geometry.dtype == np.float64
    assert pair_angle_valid.shape == (pair_index.shape[0],) and pair_angle_valid.dtype == np.uint8
    assert pair_index.shape[0] > 0
    assert len({(int(row[1]), int(row[4])) for row in pair_index}) == pair_index.shape[0]

    for row, geometry, angle_is_valid in zip(pair_index, pair_geometry, pair_angle_valid):
        donor_residue, donor_n, donor_h, acceptor_residue, acceptor_o, separation = map(int, row)
        assert int(residue_index[donor_n]) == donor_residue
        assert int(residue_index[donor_h]) == donor_residue
        assert int(residue_index[acceptor_o]) == acceptor_residue
        assert separation == backbone_separation(residues, donor_residue, acceptor_residue)

        h_to_o = pos[acceptor_o] - pos[donor_h]
        distance = np.linalg.norm(h_to_o)
        np.testing.assert_allclose(geometry[0], distance, rtol=0.0, atol=1.0e-12)
        np.testing.assert_allclose(geometry[2:5], h_to_o / distance, rtol=0.0, atol=1.0e-12)

        h_to_n = pos[donor_n] - pos[donor_h]
        if np.isfinite(h_to_n).all() and np.linalg.norm(h_to_n) > 1.0e-12:
            expected_angle = np.arccos(np.clip(
                np.dot(h_to_n, h_to_o) / (np.linalg.norm(h_to_n) * distance),
                -1.0, 1.0,
            ))
            assert angle_is_valid == 1
            np.testing.assert_allclose(geometry[1], expected_angle, rtol=0.0, atol=1.0e-12)
        else:
            assert angle_is_valid == 0 and np.isnan(geometry[1])


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: npy_local_geometry_acceptance.py STATIC_DIR FRAME_DIR")
    check_output(Path(sys.argv[1]))
    check_output(Path(sys.argv[2]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
