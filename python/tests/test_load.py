"""Tests for nmr_extract against real extraction data.

Three extraction directories:
  GEO_ONLY — geometry-only extraction (1Q8K fleet test, no MOPAC/APBS)
  BASELINE — single-protein extraction (P84477), includes MOPAC/APBS/Orca
  MUTANT   — mutant comparison (A0A7C5FAR6), includes delta arrays

GEO_ONLY is the primary test target (fast, no external tools).
BASELINE and MUTANT are secondary (may be unavailable on some machines).
"""

import numpy as np
import os
import pytest
import tempfile
import torch
from e3nn.o3 import Irreps
from pathlib import Path

from nmr_extract import (
    load,
    Protein,
    RingType,
    ShieldingTensor,
    EFGTensor,
    E3nnTensor,
    PositionField,
    VectorField,
    MagneticVectorField,
    PerRingTypeT0,
    PerRingTypeT1,
    PerRingTypeT2,
    McConnellNearFieldCounts,
    RingContributions,
    RingGeometry,
    CATALOG,
    project_t2_to_e3nn,
    project_full9_to_e3nn,
)
from nmr_extract._protein import PiQuadrupoleGroup
from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


REPO_ROOT = Path(__file__).resolve().parents[2]

# Primary: geometry-only fleet extraction (always available when fixtures are present)
GEO_ONLY = Path(os.environ.get(
    "NMR_SDK_GEO_ONLY",
    REPO_ROOT / "tests/data/sdk_geo_only/1Q8K/1Q8K_10023_01_boltzmann_minimum_8.86e-02_frame001",
))

# Secondary: legacy extractions (may not exist)
BASELINE = Path(os.environ.get(
    "NMR_BASELINE_FEATURES",
    REPO_ROOT / "baseline_features/P84477",
))
MUTANT = Path(os.environ.get(
    "NMR_SDK_MUTANT_TEST",
    Path(tempfile.gettempdir()) / "sdk_mutant_test",
))

MCCONNELL_CATEGORIES = (
    "peptide_co",
    "peptide_cn",
    "backbone_other",
    "sidechain_co",
    "sidechain_other",
    "disulfide",
    "aromatic_zeroed",
    "backbone_xh",
    "sidechain_xh",
    "s_h",
)
MCCONNELL_CHANNELS = tuple(
    f"{category}_{channel}"
    for category in MCCONNELL_CATEGORIES
    for channel in ("fixed", "bo")
)

RETIRED_NPY_STEMS = {
    "disp_per_ring_spherical",
    "disp_per_ring_tensor",
    "disp_per_type_T2",
    "disp_shielding",
    "hbond_nearest_spherical",
    "hbond_nearest_tensor",
    "hbond_shielding",
    "pq_per_type_T2",
    "pq_shielding",
    "ringchi_shielding",
}


def test_hbond_nearest_direction_catalog_names_explicit_donor_h_source():
    description = CATALOG["hbond_nearest_dir"].description
    assert "explicit donor H source point" in description
    assert "target atom" in description
    assert "midpoint" not in description


def test_pi_quadrupole_dense_aliases_are_axis_correct_and_backward_compatible(
        tmp_path):
    """A18's two names are N x 8 arrays, not sparse P-row properties."""
    legacy = PerRingTypeT0(np.arange(8, dtype=np.float64)[None, :])
    preferred = PerRingTypeT0(
        np.arange(8, dtype=np.float64)[None, :] + 100.0)

    group = PiQuadrupoleGroup(
        per_type_T0=legacy,
        axial_scalar_per_type_T0=preferred,
    )
    assert group.pq_per_type_T0 is legacy
    assert group.piquad_axial_scalar_per_type_T0 is preferred

    pre_alias = PiQuadrupoleGroup(
        per_type_T0=legacy,
        axial_scalar_per_type_T0=None,
    )
    assert pre_alias.piquad_axial_scalar_per_type_T0 is legacy

    # Exercise the real top-level required-file gate, not only direct group
    # construction: a pre-alias extraction has the legacy required payload
    # but no additive preferred-name duplicate and must still load.
    write_required_sdk_npys(tmp_path, n_atoms=1, n_residues=1)
    write_minimal_topology_sidecar(tmp_path, n_atoms=1, n_residues=1)
    preferred_path = tmp_path / "piquad_axial_scalar_per_type_T0.npy"
    assert not preferred_path.exists()
    loaded = load(tmp_path)
    assert loaded.pi_quadrupole.axial_scalar_per_type_T0 is None
    np.testing.assert_array_equal(
        loaded.pi_quadrupole.piquad_axial_scalar_per_type_T0.data,
        loaded.pi_quadrupole.pq_per_type_T0.data,
    )


def test_mcconnell_all_ten_category_views_load_from_current_contract(
        tmp_path):
    write_required_sdk_npys(tmp_path, n_atoms=1, n_residues=1)
    write_minimal_topology_sidecar(tmp_path, n_atoms=1, n_residues=1)
    mc = load(tmp_path).mcconnell
    assert set(mc.fixed) == set(MCCONNELL_CATEGORIES)
    assert set(mc.bo) == set(MCCONNELL_CATEGORIES)
    assert len(mc.fixed) == len(MCCONNELL_CATEGORIES) == 10
    assert len(mc.bo) == len(MCCONNELL_CATEGORIES) == 10
    for category in MCCONNELL_CATEGORIES:
        assert mc.fixed[category] is getattr(mc, f"{category}_fixed")
        assert mc.bo[category] is getattr(mc, f"{category}_bo")


# ── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def geo():
    if not Path(GEO_ONLY).exists():
        pytest.skip("GEO_ONLY extraction not available")
    try:
        return load(str(GEO_ONLY))
    except (ValueError, KeyError, FileNotFoundError) as e:
        pytest.skip(f"GEO_ONLY extraction stale (re-emit needed): {e}")


@pytest.fixture(scope="module")
def baseline():
    if not Path(BASELINE).exists():
        pytest.skip("baseline extraction not available")
    try:
        return load(str(BASELINE))
    except (ValueError, KeyError, FileNotFoundError) as e:
        pytest.skip(f"baseline extraction stale: {e}")


@pytest.fixture(scope="module")
def mutant():
    if not Path(MUTANT).exists():
        pytest.skip("mutant extraction not available")
    try:
        return load(str(MUTANT))
    except (ValueError, KeyError, FileNotFoundError) as e:
        pytest.skip(f"mutant extraction stale: {e}")


# ── Identity ────────────────────────────────────────────────────────


class TestIdentity:

    def test_protein_id(self, geo):
        assert "1Q8K" in geo.protein_id

    def test_n_atoms_consistent(self, geo):
        N = geo.n_atoms
        assert geo.pos.data.shape == (N, 3)
        assert geo.element.shape == (N,)
        assert geo.residue_type.shape == (N,)
        assert geo.residue_index.shape == (N,)

    def test_element_values(self, geo):
        valid = {1, 6, 7, 8, 16}
        assert set(geo.element.tolist()).issubset(valid)

    def test_pos_is_affine_positionfield(self, geo):
        assert isinstance(geo.pos, PositionField)
        assert not hasattr(geo.pos, "irreps")

    def test_large_protein(self, geo):
        """1Q8K has 4876 atoms — this is a real fleet protein."""
        assert geo.n_atoms == 4876


# ── Irreps on every tensor type ─────────────────────────────────────


class TestIrreps:

    def test_spherical_tensor(self, geo):
        st = geo.biot_savart.shielding
        assert isinstance(st, ShieldingTensor)
        assert not hasattr(st, "irreps")
        assert st.tensor_basis == "project_native_full9_spherical_tensor_v1"
        assert st.component_order == \
            "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        e3 = st.to_e3nn()
        assert isinstance(e3, E3nnTensor)
        assert e3.irreps == Irreps("1x0e + 1x1e + 1x2e")
        assert e3.irreps.dim == 9
        assert st.to_e3nn_T2().irreps == Irreps("1x2e")

    def test_efg_tensor(self, geo):
        efg = geo.coulomb.efg_backbone
        assert isinstance(efg, EFGTensor)
        assert not hasattr(efg, "irreps")
        assert efg.tensor_basis == "project_native_t2_isometric_real_tesseral_v1"
        assert efg.component_order == "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        assert efg.to_e3nn().irreps == Irreps("1x2e")
        assert efg.to_e3nn_T2().irreps == Irreps("1x2e")

    def test_vector_field(self, geo):
        v = geo.coulomb.E
        assert isinstance(v, VectorField)
        assert v.irreps == Irreps("1x1o")

    def test_per_ring_type_T0(self, geo):
        t0 = geo.biot_savart.per_type_T0
        assert isinstance(t0, PerRingTypeT0)
        assert t0.irreps == Irreps("8x0e")

    def test_per_ring_type_T2(self, geo):
        t2 = geo.biot_savart.per_type_T2
        assert isinstance(t2, PerRingTypeT2)
        assert not hasattr(t2, "irreps")
        assert t2.tensor_basis == "project_native_t2_isometric_real_tesseral_v1"
        assert t2.component_order == "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
        assert t2.to_e3nn().irreps == Irreps("8x2e")
        assert t2.to_e3nn().irreps.dim == 40

    def test_per_ring_type_T1(self, geo):
        t1 = geo.biot_savart.per_type_T1
        assert isinstance(t1, PerRingTypeT1)
        assert t1.irreps == Irreps("8x1e")
        assert t1.irreps.dim == 24

    def test_magnetic_vector_field(self, geo):
        b = geo.biot_savart.total_B
        assert isinstance(b, MagneticVectorField)
        assert b.irreps == Irreps("1x1e")

    def test_mcconnell_channels(self, geo):
        for name in MCCONNELL_CHANNELS:
            channel = getattr(geo.mcconnell, name)
            assert isinstance(channel, ShieldingTensor)
            assert not hasattr(channel, "irreps")
            assert channel.to_e3nn().irreps == Irreps("1x0e + 1x1e + 1x2e")
            assert channel.to_e3nn().irreps.dim == 9


# ── Torch conversion ────────────────────────────────────────────────


class TestTorch:

    def test_shielding_to_torch(self, geo):
        t = geo.biot_savart.shielding.torch()
        assert isinstance(t, torch.Tensor)
        assert t.dtype == torch.float64
        assert t.shape == (geo.n_atoms, 9)

    def test_vectorfield_to_torch(self, geo):
        t = geo.biot_savart.total_B.torch()
        assert isinstance(t, torch.Tensor)
        assert t.shape == (geo.n_atoms, 3)

    def test_torch_shares_memory(self, geo):
        t = geo.biot_savart.shielding.torch()
        d = geo.biot_savart.shielding.data
        assert t.numpy().base is d or np.shares_memory(t.numpy(), d)


# ── SphericalTensor components ──────────────────────────────────────


class TestSphericalComponents:

    def test_T0_shape(self, geo):
        assert geo.biot_savart.shielding.T0.shape == (geo.n_atoms, 1)

    def test_T1_shape(self, geo):
        assert geo.biot_savart.shielding.T1.shape == (geo.n_atoms, 3)

    def test_T2_shape(self, geo):
        assert geo.biot_savart.shielding.T2.shape == (geo.n_atoms, 5)

    def test_components_partition_data(self, geo):
        st = geo.biot_savart.shielding
        reconstructed = np.concatenate([st.T0, st.T1, st.T2], axis=-1)
        np.testing.assert_array_equal(reconstructed, st.data)

    def test_isotropic_is_T0_squeezed(self, geo):
        st = geo.biot_savart.shielding
        np.testing.assert_array_equal(st.isotropic, st.T0[:, 0])

    def test_T2_magnitude(self, geo):
        st = geo.biot_savart.shielding
        expected = np.linalg.norm(st.T2, axis=-1)
        np.testing.assert_allclose(st.T2_magnitude, expected)


# ── Per-ring-type T2 ────────────────────────────────────────────────


class TestPerRingTypeT2:

    def test_for_type_shape(self, geo):
        t2 = geo.biot_savart.per_type_T2
        for rt in RingType:
            assert t2.for_type(rt).shape == (geo.n_atoms, 5)

    def test_as_block_shape(self, geo):
        t2 = geo.biot_savart.per_type_T2
        assert t2.as_block().shape == (geo.n_atoms, 8, 5)

    def test_total_is_sum_of_types(self, geo):
        t2 = geo.biot_savart.per_type_T2
        manual_sum = sum(t2.for_type(rt) for rt in RingType)
        np.testing.assert_allclose(t2.total, manual_sum, atol=1e-14)

    def test_block_sum_matches_total(self, geo):
        t2 = geo.biot_savart.per_type_T2
        np.testing.assert_allclose(t2.as_block().sum(axis=1), t2.total, atol=1e-14)


# ── Per-ring-type T1 ────────────────────────────────────────────────


class TestPerRingTypeT1:

    @pytest.mark.parametrize("group_name", ["biot_savart", "haigh_mallion"])
    def test_wrapper_and_shapes(self, geo, group_name):
        t1 = getattr(geo, group_name).per_type_T1
        assert isinstance(t1, PerRingTypeT1)
        assert t1.data.shape == (geo.n_atoms, 24)
        for rt in RingType:
            assert t1.for_type(rt).shape == (geo.n_atoms, 3)
        assert t1.as_block().shape == (geo.n_atoms, 8, 3)

    @pytest.mark.parametrize("group_name", ["biot_savart", "haigh_mallion"])
    def test_components_reconstruct_shielding(self, geo, group_name):
        group = getattr(geo, group_name)
        t0 = group.per_type_T0
        t1 = group.per_type_T1
        t2 = group.per_type_T2

        assert isinstance(t0, PerRingTypeT0)
        assert isinstance(t1, PerRingTypeT1)
        assert isinstance(t2, PerRingTypeT2)
        assert t0.data.shape == (geo.n_atoms, 8)
        assert t1.data.shape == (geo.n_atoms, 24)
        assert t2.data.shape == (geo.n_atoms, 40)

        reconstructed = np.concatenate(
            [t0.total[:, None], t1.total, t2.total],
            axis=-1,
        )
        np.testing.assert_allclose(reconstructed, group.shielding.data, atol=1e-14)


# ── McConnell source channels ───────────────────────────────────────


class TestMcConnellChannels:

    def test_channel_shapes(self, geo):
        for name in MCCONNELL_CHANNELS:
            channel = getattr(geo.mcconnell, name)
            assert channel.data.shape == (geo.n_atoms, 9)
            assert channel.T0.shape == (geo.n_atoms, 1)
            assert channel.T1.shape == (geo.n_atoms, 3)
            assert channel.T2.shape == (geo.n_atoms, 5)

    def test_fixed_and_bo_views(self, geo):
        mc = geo.mcconnell
        assert set(mc.fixed) == set(MCCONNELL_CATEGORIES)
        assert set(mc.bo) == set(MCCONNELL_CATEGORIES)
        assert len(mc.fixed) == len(MCCONNELL_CATEGORIES) == 10
        assert len(mc.bo) == len(MCCONNELL_CATEGORIES) == 10
        for category in MCCONNELL_CATEGORIES:
            assert mc.fixed[category] is getattr(mc, f"{category}_fixed")
            assert mc.bo[category] is getattr(mc, f"{category}_bo")


# ── Ring contributions (sparse per-ring) ────────────────────────────


class TestRingContributions:

    def test_type(self, geo):
        assert isinstance(geo.ring_contributions, RingContributions)

    def test_atom_indices_valid(self, geo):
        rc = geo.ring_contributions
        assert rc.atom_index.min() >= 0
        assert rc.atom_index.max() < geo.n_atoms

    def test_ring_types_valid(self, geo):
        rc = geo.ring_contributions
        assert rc.ring_type.min() >= 0
        assert rc.ring_type.max() < 8

    def test_distances_positive(self, geo):
        rc = geo.ring_contributions
        assert (rc.distance > 0).all()

    def test_kernel_views_are_spherical_tensors(self, geo):
        rc = geo.ring_contributions
        assert rc.data.shape == (rc.n_pairs, 40)
        for kernel in [rc.bs, rc.hm_H, rc.hm]:
            assert hasattr(kernel, 'T2')
            assert kernel.data.shape == (rc.n_pairs, 9)
            assert kernel.T2.shape == (rc.n_pairs, 5)

    def test_for_atom_filters(self, geo):
        rc = geo.ring_contributions
        atom_0 = rc.for_atom(0)
        assert atom_0.n_pairs <= rc.n_pairs
        assert (atom_0.atom_index == 0).all()

    def test_for_ring_type_filters(self, geo):
        rc = geo.ring_contributions
        phe = rc.for_ring_type(RingType.PHE)
        assert (phe.ring_type == int(RingType.PHE)).all()

    def test_sparse_reconstructs_per_type_bs(self, geo):
        rc = geo.ring_contributions
        for rt in RingType:
            subset = rc.for_ring_type(rt)
            if subset.n_pairs == 0:
                continue
            sparse_sum = np.zeros((geo.n_atoms, 5))
            np.add.at(sparse_sum, subset.atom_index, subset.bs.T2)
            presummed = geo.biot_savart.per_type_T2.for_type(rt)
            np.testing.assert_allclose(sparse_sum, presummed, atol=1e-14)


# ── Ring geometry ────────────────────────────────────────────────────


class TestRingGeometry:

    def test_type(self, geo):
        assert isinstance(geo.ring_geometry, RingGeometry)

    def test_normals_unit_length(self, geo):
        norms = np.linalg.norm(geo.ring_geometry.normal, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_radius_positive(self, geo):
        assert (geo.ring_geometry.radius > 0).all()

    def test_ring_types_valid(self, geo):
        rt = geo.ring_geometry.ring_type
        assert rt.min() >= 0
        assert rt.max() < 8

    def test_ring_count_matches_contributions(self, geo):
        rc_rings = set(geo.ring_contributions.ring_index.tolist())
        rg_rings = set(geo.ring_geometry.ring_index.tolist())
        assert rc_rings.issubset(rg_rings)


# ── Calculator groups (geometry-only) ───────────────────────────────


class TestBiotSavart:

    def test_all_fields_present(self, geo):
        bs = geo.biot_savart
        assert bs.shielding is not None
        assert bs.per_type_T0 is not None
        assert bs.per_type_T1 is not None
        assert bs.per_type_T2 is not None
        assert bs.total_B is not None
        assert bs.ring_B_field is not None
        assert bs.ring_B_cylindrical is not None
        assert bs.ring_counts is not None

    def test_ring_counts_monotone(self, geo):
        rc = geo.biot_savart.ring_counts
        assert (rc.within_3A <= rc.within_5A + 1e-10).all()
        assert (rc.within_5A <= rc.within_8A + 1e-10).all()
        assert (rc.within_8A <= rc.within_12A + 1e-10).all()

    def test_total_B_is_magnetic_vector_field(self, geo):
        b = geo.biot_savart.total_B
        assert isinstance(b, MagneticVectorField)
        assert b.irreps == Irreps("1x1e")
        assert b.data.shape == (geo.n_atoms, 3)

    def test_ring_B_fields_are_sparse_magnetic_vectors(self, geo):
        n_pairs = len(geo.ring_contributions.atom_index)
        field = geo.biot_savart.ring_B_field
        assert isinstance(field, MagneticVectorField)
        assert field.irreps == Irreps("1x1e")
        assert field.data.shape == (n_pairs, 3)

        cylindrical = geo.biot_savart.ring_B_cylindrical
        assert isinstance(cylindrical, np.ndarray)
        assert cylindrical.shape == (n_pairs, 3)

    def test_ring_B_fields_row_align_with_contributions(self, geo):
        n_pairs = len(geo.ring_contributions.atom_index)
        assert len(geo.biot_savart.ring_B_field.data) == n_pairs
        assert len(geo.biot_savart.ring_B_cylindrical) == n_pairs


class TestHaighMallion:

    def test_all_fields_present(self, geo):
        hm = geo.haigh_mallion
        assert hm.shielding is not None
        assert hm.per_type_T0 is not None
        assert hm.per_type_T1 is not None
        assert hm.per_type_T2 is not None
        assert hm.ring_B_field is not None

    def test_ring_B_field_is_sparse_magnetic_vector(self, geo):
        n_pairs = len(geo.ring_contributions.atom_index)
        field = geo.haigh_mallion.ring_B_field
        assert isinstance(field, MagneticVectorField)
        assert field.irreps == Irreps("1x1e")
        assert field.data.shape == (n_pairs, 3)

    def test_ring_B_field_row_aligns_with_contributions(self, geo):
        assert len(geo.haigh_mallion.ring_B_field.data) == \
            len(geo.ring_contributions.atom_index)


class TestMcConnell:

    def test_all_fields_present(self, geo):
        mc = geo.mcconnell
        for name in MCCONNELL_CHANNELS:
            assert getattr(mc, name) is not None
        assert mc.nearfield_counts is not None

    def test_nearfield_counts_shape(self, geo):
        counts = geo.mcconnell.nearfield_counts
        assert isinstance(counts, McConnellNearFieldCounts)
        assert counts.data.shape == (geo.n_atoms, 2)


class TestCoulomb:

    def test_all_fields_present(self, geo):
        c = geo.coulomb
        assert c.efg is not None
        assert c.E is not None
        assert c.efg_backbone is not None
        assert c.efg_aromatic is not None
        assert c.scalars is not None

    def test_E_is_vectorfield(self, geo):
        assert isinstance(geo.coulomb.E, VectorField)

    def test_efg_is_tensor(self, geo):
        assert isinstance(geo.coulomb.efg_backbone, EFGTensor)
        assert isinstance(geo.coulomb.efg_aromatic, EFGTensor)

    def test_total_efg_t2_companion_when_fixture_reemitted(self, geo):
        assert isinstance(geo.coulomb.efg, ShieldingTensor)
        assert CATALOG["coulomb_efg_t2"].cols == 5
        assert CATALOG["mopac_coulomb_efg"].cols == 9
        if geo.coulomb.efg_t2 is None:
            pytest.skip("geo fixture predates coulomb_efg_t2.npy")
        assert isinstance(geo.coulomb.efg_t2, EFGTensor)
        np.testing.assert_allclose(geo.coulomb.efg_t2.data, geo.coulomb.efg.T2)


class TestHBond:

    def test_all_fields_present(self, geo):
        assert geo.hbond.scalars is not None
        assert geo.hbond.nearest_dir is not None
        assert geo.hbond.flags is not None

    def test_scalars_shape(self, geo):
        assert geo.hbond.scalars.data.shape == (geo.n_atoms, 4)


class TestDssp:

    def test_shape(self, geo):
        assert geo.dssp.data.shape == (geo.n_atoms, 5)


class TestRingKernelGroups:

    @pytest.mark.parametrize("group_name", [
        "haigh_mallion",
    ])
    def test_group_structure(self, geo, group_name):
        g = getattr(geo, group_name)
        N = geo.n_atoms
        assert g.shielding.data.shape == (N, 9)
        assert g.per_type_T0.data.shape == (N, 8)
        assert g.per_type_T2.data.shape == (N, 40)

    def test_scalar_ring_groups(self, geo):
        N = geo.n_atoms
        assert geo.pi_quadrupole.per_type_T0.data.shape == (N, 8)
        assert geo.dispersion.per_type_T0.data.shape == (N, 8)
        alias = geo.dispersion.aromatic_r6_proximity_per_type_T0
        if alias is None:
            pytest.skip(
                "geo fixture predates aromatic_r6_proximity_per_type_T0.npy"
            )
        assert alias.data.shape == (N, 8)
        assert geo.pi_quadrupole.quad_scalar.shape[0] == \
            geo.ring_contributions.n_pairs


# ── Optional groups: absent in geometry-only ────────────────────────


class TestOptionalAbsent:

    def test_mopac_absent(self, geo):
        assert geo.mopac is None

    def test_apbs_absent(self, geo):
        assert geo.apbs is None

    def test_orca_absent(self, geo):
        assert geo.orca is None

    def test_delta_absent(self, geo):
        assert geo.delta is None

    # test_aimnet2_absent retired 2026-05-04: AIMNet2 is no longer
    # optional per project_aimnet2_contract_20260426. AIMNet2 is part
    # of production output across all use cases on machines with CUDA
    # (every machine in scope). Positive assertion lives in
    # TestRequired.test_aimnet2_present below.


# ── Optional: MOPAC (only if baseline available) ───────────────────


class TestMopac:

    def test_present_in_baseline(self, baseline):
        assert baseline.mopac is not None

    def test_coulomb(self, baseline):
        mc = baseline.mopac.coulomb
        assert mc.efg.data.shape[-1] == 9
        assert mc.E.data.shape[-1] == 3


# ── Optional: Mutant delta ──────────────────────────────────────────


class TestDelta:

    def test_present_in_mutant(self, mutant):
        assert mutant.delta is not None

    def test_shielding(self, mutant):
        assert mutant.delta.shielding.data.shape == (mutant.n_atoms, 9)
        assert isinstance(mutant.delta.shielding, ShieldingTensor)


# ── Required: AIMNet2 (project_aimnet2_contract_20260426) ───────────


class TestRequired:

    def test_aimnet2_present(self, geo):
        """AIMNet2 must be in production output. The geo-only fixture is
        regenerated with AIMNet2 included; absence indicates either a
        stale fixture or a regression."""
        assert geo.aimnet2 is not None

    def test_aimnet2_charges_shape(self, geo):
        assert geo.aimnet2.charges.data.shape == (geo.n_atoms,)

    def test_aimnet2_aim_shape(self, geo):
        assert geo.aimnet2.aim.data.shape == (geo.n_atoms, 256)


# ── Catalog ──────────────────────────────────────────────────────────


class TestCatalog:

    def test_all_geo_only_files_registered(self):
        on_disk = {p.stem for p in Path(GEO_ONLY).glob("*.npy")}
        unregistered = on_disk - set(CATALOG.keys()) - RETIRED_NPY_STEMS
        assert not unregistered, f"Unregistered: {unregistered}"

    def test_required_count(self):
        required = [s for s in CATALOG.values() if s.required]
        assert len(required) >= 27  # Static catalog entries remain absence-compatible.

    def test_every_spec_has_description(self):
        for stem, spec in CATALOG.items():
            assert spec.description, f"{stem} has no description"

    def test_aimnet2_specs_registered(self):
        """The original AIMNet2 NPY files remain in the expanded catalog."""
        expected = {
            "aimnet2_charges", "aimnet2_aim", "aimnet2_efg",
            "aimnet2_efg_aromatic", "aimnet2_efg_backbone",
        }
        assert expected.issubset(set(CATALOG.keys())), \
            f"Missing: {expected - set(CATALOG.keys())}"
