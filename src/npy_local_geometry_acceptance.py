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

ROTATION_ANGLE = 0.731
ROTATION_AXIS = np.array([0.37, -0.61, 0.704], dtype=np.float64)
ROTATION_AXIS /= np.linalg.norm(ROTATION_AXIS)
TRANSLATION = np.array([3.25, -1.75, 2.5], dtype=np.float64)


def acceptance_rotation() -> np.ndarray:
    """Rodrigues matrix matching DirectionalAcceptanceRotation() in C++."""
    x, y, z = ROTATION_AXIS
    cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    c = np.cos(ROTATION_ANGLE)
    s = np.sin(ROTATION_ANGLE)
    return c * np.eye(3) + s * cross + (1.0 - c) * np.outer(
        ROTATION_AXIS, ROTATION_AXIS
    )


def load(path: Path) -> np.ndarray:
    return np.load(path, allow_pickle=False)


def reconstruct_full9(values: np.ndarray) -> np.ndarray:
    """Project-native [T0,T1_xyz,T2_-2..+2] -> Cartesian 3x3."""
    flat = np.asarray(values, dtype=np.float64).reshape(-1, 9)
    out = np.empty((flat.shape[0], 3, 3), dtype=np.float64)
    t0 = flat[:, 0]
    t1 = flat[:, 1:4]
    t2 = flat[:, 4:9]
    sxy = t2[:, 0] / np.sqrt(2.0)
    syz = t2[:, 1] / np.sqrt(2.0)
    szz = t2[:, 2] / np.sqrt(1.5)
    sxz = t2[:, 3] / np.sqrt(2.0)
    dxy = t2[:, 4] * np.sqrt(2.0)
    sxx = (-szz + dxy) / 2.0
    syy = (-szz - dxy) / 2.0
    out[:, 0, 0] = sxx + t0
    out[:, 1, 1] = syy + t0
    out[:, 2, 2] = szz + t0
    out[:, 0, 1] = sxy + t1[:, 2]
    out[:, 1, 0] = sxy - t1[:, 2]
    out[:, 0, 2] = sxz - t1[:, 1]
    out[:, 2, 0] = sxz + t1[:, 1]
    out[:, 1, 2] = syz + t1[:, 0]
    out[:, 2, 1] = syz - t1[:, 0]
    return out


def reconstruct_t2(values: np.ndarray) -> np.ndarray:
    flat = np.asarray(values, dtype=np.float64).reshape(-1, 5)
    full = np.zeros((flat.shape[0], 9), dtype=np.float64)
    full[:, 4:] = flat
    return reconstruct_full9(full)


def rotate_rank2(values: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    return np.einsum("ai,nij,bj->nab", rotation, values, rotation)


def assert_same_finite_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left_finite = np.isfinite(left).all(axis=tuple(range(1, left.ndim)))
    right_finite = np.isfinite(right).all(axis=tuple(range(1, right.ndim)))
    np.testing.assert_array_equal(left_finite, right_finite)
    return left_finite


def assert_vector_covariance(
    original: np.ndarray,
    transformed: np.ndarray,
    rotation: np.ndarray,
    *,
    atol: float,
) -> None:
    original = np.asarray(original).reshape(-1, 3)
    transformed = np.asarray(transformed).reshape(-1, 3)
    valid = assert_same_finite_rows(original, transformed)
    np.testing.assert_allclose(
        transformed[valid], original[valid] @ rotation.T,
        rtol=2.0e-10, atol=atol,
    )
    assert np.isnan(transformed[~valid]).all()


def assert_full9_covariance(
    original: np.ndarray,
    transformed: np.ndarray,
    rotation: np.ndarray,
    *,
    atol: float,
) -> None:
    original = np.asarray(original).reshape(-1, 9)
    transformed = np.asarray(transformed).reshape(-1, 9)
    valid = assert_same_finite_rows(original, transformed)
    lhs = reconstruct_full9(transformed[valid])
    rhs = rotate_rank2(reconstruct_full9(original[valid]), rotation)
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-9, atol=atol)
    assert np.isnan(transformed[~valid]).all()


def assert_t2_covariance(
    original: np.ndarray,
    transformed: np.ndarray,
    rotation: np.ndarray,
    *,
    atol: float,
) -> None:
    original = np.asarray(original).reshape(-1, 5)
    transformed = np.asarray(transformed).reshape(-1, 5)
    valid = assert_same_finite_rows(original, transformed)
    lhs = reconstruct_t2(transformed[valid])
    rhs = rotate_rank2(reconstruct_t2(original[valid]), rotation)
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-9, atol=atol)
    assert np.isnan(transformed[~valid]).all()


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
    # Exercise the real final reader on every serialized producer array. This
    # catches accidental object/pickle payloads even for invariant sidecars;
    # the directional arrays receive stronger name-level checks below.
    serialized = {path.stem: load(path) for path in sorted(output_dir.glob("*.npy"))}
    assert serialized, output_dir

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


def check_directional_covariance(original_dir: Path, transformed_dir: Path) -> None:
    """Compare two complete calculator reruns at the serialized NPY boundary."""
    q = acceptance_rotation()
    np.testing.assert_allclose(q @ q.T, np.eye(3), rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(np.linalg.det(q), 1.0, rtol=0.0, atol=2.0e-15)

    original_pos = load(original_dir / "pos.npy")
    transformed_pos = load(transformed_dir / "pos.npy")
    np.testing.assert_allclose(
        transformed_pos, original_pos @ q.T + TRANSLATION,
        rtol=0.0, atol=3.0e-12,
    )

    def pair(stem: str) -> tuple[np.ndarray, np.ndarray]:
        left_path = original_dir / f"{stem}.npy"
        right_path = transformed_dir / f"{stem}.npy"
        assert left_path.is_file() == right_path.is_file(), stem
        return load(left_path), load(right_path)

    def vector(stem: str, atol: float = 2.0e-10) -> None:
        left, right = pair(stem)
        assert_vector_covariance(left, right, q, atol=atol)

    for stem in (
        "bond_direction",
        "ring_direction_to_center",
        "bs_total_B",
        "bs_ring_B_field",
        "hm_ring_B_field",
        "hbond_nearest_dir",
        "cb_residual_vector",
        "eeq_coulomb_E",
        "eeq_coulomb_E_backbone",
        "eeq_coulomb_E_sidechain",
        "eeq_coulomb_E_aromatic",
    ):
        vector(stem)

    # The fixed-sign cross(c,n) ideal-L-CB construction has an SO(3), not an
    # O(3), contract.  This production rerun uses a proper transform, under
    # which its scalar norm is invariant and its global residual is polar.
    cb_deviation, cb_deviation_q = pair("cb_deviation")
    np.testing.assert_allclose(
        cb_deviation_q, cb_deviation, rtol=2.0e-10, atol=2.0e-10,
        equal_nan=True,
    )

    mc_direction, mc_direction_q = pair("mc_nearest_co_dir")
    assert_vector_covariance(mc_direction, mc_direction_q, q, atol=3.0e-10)
    mc_midpoint, mc_midpoint_q = pair("mc_nearest_co_midpoint")
    midpoint_valid = assert_same_finite_rows(mc_midpoint, mc_midpoint_q)
    np.testing.assert_allclose(
        mc_midpoint_q[midpoint_valid],
        mc_midpoint[midpoint_valid] @ q.T + TRANSLATION,
        rtol=0.0, atol=3.0e-10,
    )
    assert np.isnan(mc_midpoint_q[~midpoint_valid]).all()

    # Sparse rows retain identity/order while their vector blocks rotate.
    spatial, spatial_q = pair("spatial_neighbors")
    np.testing.assert_array_equal(spatial_q[:, :2], spatial[:, :2])
    assert_vector_covariance(spatial[:, 2:5], spatial_q[:, 2:5], q, atol=2.0e-12)
    np.testing.assert_allclose(spatial_q[:, 5], spatial[:, 5], rtol=0.0, atol=2.0e-12)

    hbond_index, hbond_index_q = pair("hbond_pairs_index")
    np.testing.assert_array_equal(hbond_index_q, hbond_index)
    hbond, hbond_q = pair("hbond_pairs_geometry")
    np.testing.assert_allclose(hbond_q[:, :2], hbond[:, :2], rtol=0.0, atol=2.0e-11)
    assert_vector_covariance(hbond[:, 2:5], hbond_q[:, 2:5], q, atol=2.0e-11)

    ring_geometry, ring_geometry_q = pair("ring_geometry")
    np.testing.assert_array_equal(ring_geometry_q[:, :3], ring_geometry[:, :3])
    np.testing.assert_allclose(
        ring_geometry_q[:, 3:6], ring_geometry[:, 3:6] @ q.T + TRANSLATION,
        rtol=0.0, atol=3.0e-12,
    )
    assert_vector_covariance(
        ring_geometry[:, 6:9], ring_geometry_q[:, 6:9], q, atol=3.0e-12
    )
    np.testing.assert_allclose(
        ring_geometry_q[:, 9], ring_geometry[:, 9], rtol=0.0, atol=2.0e-12
    )

    ring_pairs, ring_pairs_q = pair("ring_pair_geometry")
    np.testing.assert_allclose(ring_pairs_q, ring_pairs, rtol=2.0e-12, atol=2.0e-12)

    bs_cyl, bs_cyl_q = pair("bs_ring_B_cylindrical")
    np.testing.assert_allclose(bs_cyl_q, bs_cyl, rtol=3.0e-10, atol=3.0e-16)
    for stem in (
        "bs_per_type_T0", "hm_per_type_T0", "pq_per_type_T0",
        "piquad_axial_scalar_per_type_T0", "piquad_quad_scalar",
        "ringchi_scalar", "ringchi_per_type_T0",
    ):
        left, right = pair(stem)
        np.testing.assert_allclose(right, left, rtol=3.0e-10, atol=3.0e-12,
                                   equal_nan=True)

    # Full project-native tensors: reconstruct before applying Q T Q^T.
    mc_full9 = (
        "mc_peptide_co_fixed", "mc_peptide_co_bo", "mc_peptide_co_rhombic",
        "mc_peptide_cn_fixed", "mc_peptide_cn_bo",
        "mc_backbone_other_fixed", "mc_backbone_other_bo",
        "mc_sidechain_co_fixed", "mc_sidechain_co_bo",
        "mc_sidechain_other_fixed", "mc_sidechain_other_bo",
        "mc_disulfide_fixed", "mc_disulfide_bo",
        "mc_aromatic_fixed", "mc_aromatic_bo",
        "mc_backbone_xh_fixed", "mc_backbone_xh_bo",
        "mc_sidechain_xh_fixed", "mc_sidechain_xh_bo",
        "mc_s_h_fixed", "mc_s_h_bo",
        "mc_nearest_co_T2", "mc_nearest_cn_T2",
        "sidechain_co_fixed_T2", "sidechain_co_bo_T2",
    )
    for stem in ("bs_shielding", "hm_shielding", "eeq_coulomb_efg", *mc_full9):
        left, right = pair(stem)
        assert_full9_covariance(left, right, q, atol=3.0e-9)

    # Per-ring-type T1 is eight Cartesian Levi-Civita-dual blocks. For this
    # proper Q it follows Q; the independent improper test pins det(Q)Q.
    for stem in ("bs_per_type_T1", "hm_per_type_T1"):
        left, right = pair(stem)
        assert_vector_covariance(left.reshape(-1, 3), right.reshape(-1, 3), q,
                                 atol=3.0e-9)

    native_t2 = (
        "bs_per_type_T2", "hm_per_type_T2",
        "eeq_coulomb_efg_backbone", "eeq_coulomb_efg_sidechain",
        "eeq_coulomb_efg_aromatic",
    )
    for stem in native_t2:
        left, right = pair(stem)
        assert_t2_covariance(left.reshape(-1, 5), right.reshape(-1, 5), q,
                             atol=3.0e-9)

    # EFG full9 structural zeros are part of the serialized contract.
    for directory in (original_dir, transformed_dir):
        full = load(directory / "eeq_coulomb_efg.npy")
        np.testing.assert_allclose(full[:, :4], 0.0, rtol=0.0, atol=2.0e-12)

    # Shared sparse ring row identity, invariant geometry, and three tensor
    # blocks (BS G, symmetric-traceless HM H, HM G).
    ring, ring_q = pair("ring_contributions")
    np.testing.assert_array_equal(ring_q[:, :3], ring[:, :3])
    np.testing.assert_allclose(ring_q[:, 3:9], ring[:, 3:9], rtol=3.0e-10,
                               atol=3.0e-10)
    for start in (9, 18, 27):
        assert_full9_covariance(ring[:, start:start + 9],
                                ring_q[:, start:start + 9], q, atol=4.0e-9)
    np.testing.assert_allclose(ring_q[:, 36:40], ring[:, 36:40],
                               rtol=3.0e-10, atol=3.0e-10)
    np.testing.assert_allclose(ring[:, 18:22], 0.0, rtol=0.0, atol=2.0e-12)
    np.testing.assert_allclose(ring_q[:, 18:22], 0.0, rtol=0.0, atol=2.0e-12)

    # Ring-local PiQuad coefficients are invariant under a global proper
    # rotation; its serialized frame axes rotate with the conformation.
    for stem in ("piquad_local_tensor", "piquad_local_T2", "piquad_local_geometry"):
        left, right = pair(stem)
        np.testing.assert_allclose(right, left, rtol=3.0e-10, atol=3.0e-10,
                                   equal_nan=True)
    local_frame, local_frame_q = pair("piquad_local_frame")
    assert_vector_covariance(local_frame.reshape(-1, 3),
                             local_frame_q.reshape(-1, 3), q, atol=3.0e-10)

    sidechain_frame, sidechain_frame_q = pair("sidechain_co_frame")
    np.testing.assert_allclose(
        sidechain_frame_q[:, :3], sidechain_frame[:, :3] @ q.T + TRANSLATION,
        rtol=0.0, atol=3.0e-12,
    )
    assert_vector_covariance(sidechain_frame[:, 3:].reshape(-1, 3),
                             sidechain_frame_q[:, 3:].reshape(-1, 3), q,
                             atol=3.0e-10)
    for stem in ("sidechain_co_source_bonds", "sidechain_co_frame_quality"):
        left, right = pair(stem)
        np.testing.assert_allclose(right, left, rtol=3.0e-10, atol=3.0e-10,
                                   equal_nan=True)

    # Live SASA sums a lab-fixed 92-point Fibonacci lattice.  The continuum
    # normal is polar, but the finite lattice is not O(3)-closed.  On this
    # committed full-protein fixture, atoms with at least 10 A^2 exposed area
    # (several exposed sample cells, rather than a cancelling one-cell fringe)
    # have cos(Qn,n') >= 0.90.  The area bound is 6 A^2: at most four cells at
    # the largest expanded atomic radius.  Up to 5% of zero/nonzero statuses
    # may switch at the sampling boundary; the missing normal-valid mask is
    # recorded as a source-contract issue rather than hidden as exact O(3).
    sasa, sasa_q = pair("sasa_normal")
    expected_sasa = sasa @ q.T
    sasa_area, sasa_area_q = pair("atom_sasa")
    np.testing.assert_allclose(sasa_area_q, sasa_area, rtol=0.0, atol=6.0)
    sasa_fraction, sasa_fraction_q = pair("atom_sasa_fraction")
    element, element_q = pair("element")
    np.testing.assert_array_equal(element_q, element)
    # Independent serialized-value reconstruction of the producer's Bondi +
    # 1.4 A probe normalization. Unknown elements use the live 1.70 A default.
    bondi_radius = np.full(element.shape, 1.70, dtype=np.float64)
    for atomic_number, radius in {
        1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 16: 1.80,
    }.items():
        bondi_radius[element == atomic_number] = radius
    sasa_denominator = 4.0 * np.pi * (bondi_radius + 1.4) ** 2
    np.testing.assert_allclose(
        sasa_fraction, sasa_area / sasa_denominator,
        rtol=2.0e-15, atol=2.0e-15,
    )
    np.testing.assert_allclose(
        sasa_fraction_q, sasa_area_q / sasa_denominator,
        rtol=2.0e-15, atol=2.0e-15,
    )
    assert np.all(
        np.abs(sasa_fraction_q - sasa_fraction)
        <= 6.0 / sasa_denominator + 2.0e-15
    )
    original_zero = np.linalg.norm(sasa, axis=1) < 1.0e-14
    rotated_zero = np.linalg.norm(sasa_q, axis=1) < 1.0e-14
    populated = ~original_zero & ~rotated_zero
    stable = populated & (np.minimum(sasa_area, sasa_area_q) > 10.0)
    assert stable.any()
    cosine = np.sum(expected_sasa[stable] * sasa_q[stable], axis=1)
    assert np.min(cosine) >= 0.90, np.min(cosine)
    assert np.mean(original_zero != rotated_zero) <= 0.05

    # Translation-independent signed geometry is invariant for proper Q.
    for stem in (
        "omega_actual", "omega_deviation", "omega_sin", "omega_cos",
        "aromatic_chi2", "pucker_Q", "pucker_theta",
    ):
        left, right = pair(stem)
        np.testing.assert_allclose(right, left, rtol=0.0, atol=5.0e-11,
                                   equal_nan=True)
    # DSSP is rerun through its production text boundary, whose coordinate
    # formatting is lower precision than the in-memory calculators.
    for stem in ("dssp_torsion_angle", "dssp_torsion_sin", "dssp_torsion_cos"):
        left, right = pair(stem)
        np.testing.assert_allclose(right, left, rtol=0.0, atol=5.0e-2,
                                   equal_nan=True)
    # PPII is a chirality-conditioned libdssp class, not a 0e scalar.  It is
    # stable under this proper rerun; an improper transform reverses the
    # signed phi/psi region and has no homogeneous class transformation law.
    for stem in ("dssp_observed", "dssp_ss8", "dssp_ppii"):
        left, right = pair(stem)
        np.testing.assert_array_equal(right, left)


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: npy_local_geometry_acceptance.py "
            "STATIC_DIR TRANSFORMED_STATIC_DIR FRAME_DIR"
        )
    static_dir = Path(sys.argv[1])
    transformed_static_dir = Path(sys.argv[2])
    frame_dir = Path(sys.argv[3])
    check_output(static_dir)
    check_output(transformed_static_dir)
    check_output(frame_dir)
    check_directional_covariance(static_dir, transformed_static_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
