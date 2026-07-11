"""Producer-SDK read-back for typed side-chain carbonyl anisotropy data."""

import numpy as np

from nmr_extract import (
    CATALOG,
    ShieldingTensor,
    SidechainCarbonylAnisotropyGroup,
    load,
)
from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


N_ATOMS = 5
N_SOURCES = 3


def test_sidechain_carbonyl_catalog_contract():
    expected = {
        "sidechain_co_source_bonds": (8, "sidechain_co_source"),
        "sidechain_co_frame": (12, "sidechain_co_source"),
        "sidechain_co_frame_quality": (4, "sidechain_co_source"),
        "sidechain_co_fixed_T2": (9, "atom"),
        "sidechain_co_bo_T2": (9, "atom"),
        "sidechain_co_scalar_audit": (4, "atom"),
    }
    for stem, (cols, axis) in expected.items():
        spec = CATALOG[stem]
        assert spec.required is True
        assert spec.group == "sidechain_carbonyl_anisotropy"
        assert spec.cols == cols
        assert spec.native_axis == axis
        assert spec.mechanism == "bond_anisotropy"

    for stem in ("sidechain_co_fixed_T2", "sidechain_co_bo_T2"):
        spec = CATALOG[stem]
        assert spec.wrapper is ShieldingTensor
        assert spec.tensor_basis == \
            "project_native_full9_spherical_tensor_v1"
        assert spec.tensor_component_order == \
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        assert spec.tensor_frame == "conformation_cartesian_xyz"


def test_sidechain_carbonyl_group_reads_all_producer_arrays(tmp_path):
    write_required_sdk_npys(
        tmp_path,
        N_ATOMS,
        n_residues=2,
        n_sidechain_co_sources=N_SOURCES,
    )
    write_minimal_topology_sidecar(tmp_path, N_ATOMS, n_residues=2)

    source_bonds = np.asarray([
        [0, 0, 1, 0, 4, 2, 1, 1],
        [2, 2, 3, 1, 4, 7, 2, 1],
        [3, 2, 4, 1, 4, 7, 2, 1],
    ], dtype=np.int32)
    frame = np.arange(N_SOURCES * 12, dtype=np.float64).reshape(
        N_SOURCES, 12) / 10.0
    quality = np.asarray([
        [1.23, 0.0, 1.1, 1.0],
        [1.25, 1.0e-15, 1.2, 1.0],
        [1.26, 2.0e-15, 1.3, 1.0],
    ], dtype=np.float64)
    fixed = np.arange(N_ATOMS * 9, dtype=np.float64).reshape(
        N_ATOMS, 9) / 7.0
    bo = np.full((N_ATOMS, 9), np.nan, dtype=np.float64)
    audit = np.asarray([
        [1.0, np.nan, 1.0, 2.1],
        [2.0, np.nan, 2.0, 2.2],
        [3.0, np.nan, 3.0, 2.3],
        [4.0, np.nan, 0.0, np.nan],
        [5.0, np.nan, 1.0, 2.5],
    ], dtype=np.float64)

    for stem, data in {
        "sidechain_co_source_bonds": source_bonds,
        "sidechain_co_frame": frame,
        "sidechain_co_frame_quality": quality,
        "sidechain_co_fixed_T2": fixed,
        "sidechain_co_bo_T2": bo,
        "sidechain_co_scalar_audit": audit,
    }.items():
        np.save(tmp_path / f"{stem}.npy", data)

    protein = load(tmp_path)
    group = protein.sidechain_carbonyl_anisotropy
    assert isinstance(group, SidechainCarbonylAnisotropyGroup)
    assert group.source_bonds.dtype == np.int32
    np.testing.assert_array_equal(group.source_bonds, source_bonds)
    np.testing.assert_allclose(group.frame, frame)
    np.testing.assert_allclose(group.frame_quality, quality)
    assert isinstance(group.fixed_T2, ShieldingTensor)
    assert isinstance(group.bo_T2, ShieldingTensor)
    np.testing.assert_allclose(group.fixed_T2.data, fixed)
    np.testing.assert_array_equal(np.isnan(group.bo_T2.data), np.ones_like(
        bo, dtype=bool))
    np.testing.assert_allclose(group.scalar_audit, audit, equal_nan=True)
