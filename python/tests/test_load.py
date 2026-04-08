"""Tests for nmr_extract against real extraction data.

Two extraction directories:
  BASELINE — single-protein extraction (P84477), no mutation delta
  MUTANT   — mutant comparison (A0A7C5FAR6), includes delta arrays
"""

import numpy as np
import pytest
import torch
from e3nn.o3 import Irreps

from nmr_extract import (
    load,
    Protein,
    RingType,
    BondCategory,
    ShieldingTensor,
    EFGTensor,
    VectorField,
    PerRingTypeT0,
    PerRingTypeT2,
    PerBondCategoryT2,
    RingContributions,
    RingGeometry,
    CATALOG,
)


BASELINE = "/shared/2026Thesis/nmr-shielding/baseline_features/P84477"
MUTANT = "/tmp/sdk_mutant_test"


# ── Fixtures ────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def baseline():
    return load(BASELINE)


@pytest.fixture(scope="module")
def mutant():
    return load(MUTANT)


# ── Identity ────────────────────────────────────────────────────────


class TestIdentity:

    def test_protein_id(self, baseline):
        assert baseline.protein_id == "P84477"

    def test_n_atoms_consistent(self, baseline):
        N = baseline.n_atoms
        assert baseline.pos.data.shape == (N, 3)
        assert baseline.element.shape == (N,)
        assert baseline.residue_type.shape == (N,)
        assert baseline.residue_index.shape == (N,)

    def test_element_values(self, baseline):
        # Atomic numbers: H=1, C=6, N=7, O=8, S=16
        valid = {1, 6, 7, 8, 16}
        assert set(baseline.element.tolist()).issubset(valid)

    def test_pos_is_vectorfield(self, baseline):
        assert isinstance(baseline.pos, VectorField)
        assert baseline.pos.irreps == Irreps("1x1o")


# ── Irreps on every tensor type ─────────────────────────────────────


class TestIrreps:

    def test_spherical_tensor(self, baseline):
        st = baseline.biot_savart.shielding
        assert isinstance(st, ShieldingTensor)
        assert st.irreps == Irreps("1x0e + 1x1o + 1x2e")
        assert st.irreps.dim == 9

    def test_efg_tensor(self, baseline):
        efg = baseline.coulomb.efg_backbone
        assert isinstance(efg, EFGTensor)
        assert efg.irreps == Irreps("1x0e + 1x1o + 1x2e")

    def test_vector_field(self, baseline):
        v = baseline.coulomb.E
        assert isinstance(v, VectorField)
        assert v.irreps == Irreps("1x1o")

    def test_per_ring_type_T0(self, baseline):
        t0 = baseline.biot_savart.per_type_T0
        assert isinstance(t0, PerRingTypeT0)
        assert t0.irreps == Irreps("8x0e")

    def test_per_ring_type_T2(self, baseline):
        t2 = baseline.biot_savart.per_type_T2
        assert isinstance(t2, PerRingTypeT2)
        assert t2.irreps == Irreps("8x2e")
        assert t2.irreps.dim == 40

    def test_per_bond_category_T2(self, baseline):
        mc = baseline.mcconnell.category_T2
        assert isinstance(mc, PerBondCategoryT2)
        assert mc.irreps == Irreps("5x2e")
        assert mc.irreps.dim == 25


# ── Torch conversion ────────────────────────────────────────────────


class TestTorch:

    def test_shielding_to_torch(self, baseline):
        t = baseline.biot_savart.shielding.torch()
        assert isinstance(t, torch.Tensor)
        assert t.dtype == torch.float64
        assert t.shape == (baseline.n_atoms, 9)

    def test_vectorfield_to_torch(self, baseline):
        t = baseline.biot_savart.total_B.torch()
        assert isinstance(t, torch.Tensor)
        assert t.shape == (baseline.n_atoms, 3)

    def test_torch_shares_memory(self, baseline):
        t = baseline.biot_savart.shielding.torch()
        d = baseline.biot_savart.shielding.data
        # torch.from_numpy shares memory with the numpy array
        assert t.numpy().base is d or np.shares_memory(t.numpy(), d)


# ── SphericalTensor components ──────────────────────────────────────


class TestSphericalComponents:

    def test_T0_shape(self, baseline):
        assert baseline.biot_savart.shielding.T0.shape == (baseline.n_atoms, 1)

    def test_T1_shape(self, baseline):
        assert baseline.biot_savart.shielding.T1.shape == (baseline.n_atoms, 3)

    def test_T2_shape(self, baseline):
        assert baseline.biot_savart.shielding.T2.shape == (baseline.n_atoms, 5)

    def test_components_partition_data(self, baseline):
        """T0 + T1 + T2 = the full 9-component array."""
        st = baseline.biot_savart.shielding
        reconstructed = np.concatenate([st.T0, st.T1, st.T2], axis=-1)
        np.testing.assert_array_equal(reconstructed, st.data)

    def test_isotropic_is_T0_squeezed(self, baseline):
        st = baseline.biot_savart.shielding
        np.testing.assert_array_equal(st.isotropic, st.T0[:, 0])

    def test_T2_magnitude(self, baseline):
        st = baseline.biot_savart.shielding
        expected = np.linalg.norm(st.T2, axis=-1)
        np.testing.assert_allclose(st.T2_magnitude, expected)


# ── Per-ring-type T2 ────────────────────────────────────────────────


class TestPerRingTypeT2:

    def test_for_type_shape(self, baseline):
        t2 = baseline.biot_savart.per_type_T2
        for rt in RingType:
            assert t2.for_type(rt).shape == (baseline.n_atoms, 5)

    def test_as_block_shape(self, baseline):
        t2 = baseline.biot_savart.per_type_T2
        assert t2.as_block().shape == (baseline.n_atoms, 8, 5)

    def test_total_is_sum_of_types(self, baseline):
        t2 = baseline.biot_savart.per_type_T2
        manual_sum = sum(t2.for_type(rt) for rt in RingType)
        np.testing.assert_allclose(t2.total, manual_sum, atol=1e-14)

    def test_block_sum_matches_total(self, baseline):
        t2 = baseline.biot_savart.per_type_T2
        np.testing.assert_allclose(t2.as_block().sum(axis=1), t2.total, atol=1e-14)


# ── Per-bond-category T2 ────────────────────────────────────────────


class TestPerBondCategoryT2:

    def test_for_category_shape(self, baseline):
        mc = baseline.mcconnell.category_T2
        for cat in BondCategory:
            assert mc.for_category(cat).shape == (baseline.n_atoms, 5)

    def test_as_block_shape(self, baseline):
        mc = baseline.mcconnell.category_T2
        assert mc.as_block().shape == (baseline.n_atoms, 5, 5)


# ── Ring contributions (sparse per-ring) ────────────────────────────


class TestRingContributions:

    def test_type(self, baseline):
        assert isinstance(baseline.ring_contributions, RingContributions)

    def test_shape(self, baseline):
        rc = baseline.ring_contributions
        assert rc.data.shape == (rc.n_pairs, 57)

    def test_atom_indices_valid(self, baseline):
        rc = baseline.ring_contributions
        assert rc.atom_index.min() >= 0
        assert rc.atom_index.max() < baseline.n_atoms

    def test_ring_types_valid(self, baseline):
        rc = baseline.ring_contributions
        assert rc.ring_type.min() >= 0
        assert rc.ring_type.max() < 8

    def test_distances_positive(self, baseline):
        rc = baseline.ring_contributions
        assert (rc.distance > 0).all()

    def test_exp_decay_range(self, baseline):
        rc = baseline.ring_contributions
        assert (rc.exp_decay > 0).all()
        assert (rc.exp_decay <= 1).all()

    def test_kernel_views_are_spherical_tensors(self, baseline):
        rc = baseline.ring_contributions
        for kernel in [rc.bs, rc.hm_H, rc.hm, rc.pq, rc.chi]:
            assert hasattr(kernel, 'T2')
            assert kernel.data.shape == (rc.n_pairs, 9)
            assert kernel.T2.shape == (rc.n_pairs, 5)

    def test_for_atom_filters(self, baseline):
        rc = baseline.ring_contributions
        atom_0 = rc.for_atom(0)
        assert atom_0.n_pairs <= rc.n_pairs
        assert (atom_0.atom_index == 0).all()

    def test_for_ring_type_filters(self, baseline):
        rc = baseline.ring_contributions
        phe = rc.for_ring_type(RingType.PHE)
        assert (phe.ring_type == int(RingType.PHE)).all()

    def test_sparse_reconstructs_per_type_bs(self, baseline):
        """Summing sparse BS T2 by ring type must match per_type_T2."""
        rc = baseline.ring_contributions
        for rt in RingType:
            subset = rc.for_ring_type(rt)
            if subset.n_pairs == 0:
                continue
            sparse_sum = np.zeros((baseline.n_atoms, 5))
            np.add.at(sparse_sum, subset.atom_index, subset.bs.T2)
            presummed = baseline.biot_savart.per_type_T2.for_type(rt)
            np.testing.assert_allclose(sparse_sum, presummed, atol=1e-14)

    def test_sparse_reconstructs_per_type_hm(self, baseline):
        """hm (G, shielding kernel) sums must match per_type_T2."""
        rc = baseline.ring_contributions
        for rt in RingType:
            subset = rc.for_ring_type(rt)
            if subset.n_pairs == 0:
                continue
            sparse_sum = np.zeros((baseline.n_atoms, 5))
            np.add.at(sparse_sum, subset.atom_index, subset.hm.T2)
            presummed = baseline.haigh_mallion.per_type_T2.for_type(rt)
            np.testing.assert_allclose(sparse_sum, presummed, atol=1e-14)

    def test_hm_H_differs_from_hm_G(self, baseline):
        """hm_H (raw integral) and hm (shielding kernel) differ by intensity."""
        rc = baseline.ring_contributions
        if rc.n_pairs == 0:
            pytest.skip("No ring pairs")
        assert not np.allclose(rc.hm_H.T2, rc.hm.T2, atol=1e-6)


# ── Ring geometry ────────────────────────────────────────────────────


class TestRingGeometry:

    def test_type(self, baseline):
        assert isinstance(baseline.ring_geometry, RingGeometry)

    def test_shape(self, baseline):
        rg = baseline.ring_geometry
        assert rg.data.shape == (rg.n_rings, 10)

    def test_normals_unit_length(self, baseline):
        norms = np.linalg.norm(baseline.ring_geometry.normal, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_radius_positive(self, baseline):
        assert (baseline.ring_geometry.radius > 0).all()

    def test_ring_types_valid(self, baseline):
        rt = baseline.ring_geometry.ring_type
        assert rt.min() >= 0
        assert rt.max() < 8

    def test_ring_index_sequential(self, baseline):
        rg = baseline.ring_geometry
        np.testing.assert_array_equal(rg.ring_index, np.arange(rg.n_rings))

    def test_ring_count_matches_contributions(self, baseline):
        """Ring geometry table covers all rings referenced in contributions."""
        rc_rings = set(baseline.ring_contributions.ring_index.tolist())
        rg_rings = set(baseline.ring_geometry.ring_index.tolist())
        assert rc_rings.issubset(rg_rings)


# ── Calculator groups ────────────────────────────────────────────────


class TestBiotSavart:

    def test_all_fields_present(self, baseline):
        bs = baseline.biot_savart
        assert bs.shielding is not None
        assert bs.per_type_T0 is not None
        assert bs.per_type_T2 is not None
        assert bs.total_B is not None
        assert bs.ring_counts is not None

    def test_ring_counts_nonnegative(self, baseline):
        rc = baseline.biot_savart.ring_counts
        assert (rc.data >= 0).all()

    def test_ring_counts_monotone(self, baseline):
        """within_3A <= within_5A <= within_8A <= within_12A."""
        rc = baseline.biot_savart.ring_counts
        assert (rc.within_3A <= rc.within_5A + 1e-10).all()
        assert (rc.within_5A <= rc.within_8A + 1e-10).all()
        assert (rc.within_8A <= rc.within_12A + 1e-10).all()


class TestMcConnell:

    def test_all_fields_present(self, baseline):
        mc = baseline.mcconnell
        assert mc.shielding is not None
        assert mc.category_T2 is not None
        assert mc.scalars is not None

    def test_scalars_shape(self, baseline):
        assert baseline.mcconnell.scalars.data.shape == (baseline.n_atoms, 6)

    def test_nearest_distances_positive(self, baseline):
        s = baseline.mcconnell.scalars
        # Sentinel values may be large, but should not be negative
        assert (s.nearest_CO_dist >= 0).all()
        assert (s.nearest_CN_dist >= 0).all()


class TestCoulomb:

    def test_all_fields_present(self, baseline):
        c = baseline.coulomb
        assert c.shielding is not None
        assert c.E is not None
        assert c.efg_backbone is not None
        assert c.efg_aromatic is not None
        assert c.scalars is not None

    def test_E_is_vectorfield(self, baseline):
        assert isinstance(baseline.coulomb.E, VectorField)

    def test_efg_is_tensor(self, baseline):
        assert isinstance(baseline.coulomb.efg_backbone, EFGTensor)
        assert isinstance(baseline.coulomb.efg_aromatic, EFGTensor)


class TestHBond:

    def test_all_fields_present(self, baseline):
        assert baseline.hbond.shielding is not None
        assert baseline.hbond.scalars is not None

    def test_scalars_shape(self, baseline):
        assert baseline.hbond.scalars.data.shape == (baseline.n_atoms, 3)


class TestRingSusceptibility:

    def test_is_shielding_tensor(self, baseline):
        assert isinstance(baseline.ring_susceptibility, ShieldingTensor)
        assert baseline.ring_susceptibility.data.shape == (baseline.n_atoms, 9)


class TestDssp:

    def test_shape(self, baseline):
        assert baseline.dssp.data.shape == (baseline.n_atoms, 5)


# ── All four ring kernel groups ──────────────────────────────────────


class TestRingKernelGroups:
    """Haigh-Mallion, pi-quadrupole, dispersion share the same structure."""

    @pytest.mark.parametrize("group_name", [
        "haigh_mallion", "pi_quadrupole", "dispersion",
    ])
    def test_group_structure(self, baseline, group_name):
        g = getattr(baseline, group_name)
        N = baseline.n_atoms
        assert g.shielding.data.shape == (N, 9)
        assert g.per_type_T0.data.shape == (N, 8)
        assert g.per_type_T2.data.shape == (N, 40)


# ── Optional: MOPAC ──────────────────────────────────────────────────


class TestMopac:

    def test_present_in_baseline(self, baseline):
        assert baseline.mopac is not None

    def test_core(self, baseline):
        core = baseline.mopac.core
        assert core.charges.shape == (baseline.n_atoms,)
        assert core.scalars.data.shape == (baseline.n_atoms, 3)
        assert core.bond_orders.n_bonds > 0
        assert core.global_.data.shape == (4,)

    def test_coulomb(self, baseline):
        mc = baseline.mopac.coulomb
        assert mc.shielding.data.shape[-1] == 9
        assert mc.E.data.shape[-1] == 3

    def test_mcconnell(self, baseline):
        mc = baseline.mopac.mcconnell
        assert mc.shielding.data.shape[-1] == 9
        assert mc.category_T2.data.shape[-1] == 25


# ── Optional: Orca DFT ──────────────────────────────────────────────


class TestOrca:

    def test_present_in_baseline(self, baseline):
        assert baseline.orca is not None

    def test_all_three(self, baseline):
        o = baseline.orca
        N = baseline.n_atoms
        assert o.total.data.shape == (N, 9)
        assert o.diamagnetic.data.shape == (N, 9)
        assert o.paramagnetic.data.shape == (N, 9)


# ── Optional: APBS ──────────────────────────────────────────────────


class TestAPBS:

    def test_present_in_baseline(self, baseline):
        assert baseline.apbs is not None

    def test_fields(self, baseline):
        a = baseline.apbs
        assert a.E.data.shape[-1] == 3
        assert a.efg.data.shape[-1] == 9


# ── Mutant delta ─────────────────────────────────────────────────────


class TestDelta:

    def test_absent_in_baseline(self, baseline):
        assert baseline.delta is None

    def test_present_in_mutant(self, mutant):
        assert mutant.delta is not None

    def test_shielding(self, mutant):
        assert mutant.delta.shielding.data.shape == (mutant.n_atoms, 9)
        assert isinstance(mutant.delta.shielding, ShieldingTensor)

    def test_scalars(self, mutant):
        ds = mutant.delta.scalars
        assert ds.data.shape == (mutant.n_atoms, 6)
        n_matched = ds.matched_mask.sum()
        assert 0 < n_matched < mutant.n_atoms

    def test_ring_proximity(self, mutant):
        drp = mutant.delta.ring_proximity
        assert drp.n_removed_rings > 0
        d0 = drp.distance(0)
        assert d0.shape == (mutant.n_atoms,)
        assert (d0 >= 0).all()  # atoms on the mutated residue have distance 0

    def test_delta_apbs(self, mutant):
        if mutant.delta.apbs is None:
            pytest.skip("No APBS delta in this extraction")
        a = mutant.delta.apbs
        assert a.delta_E.data.shape[-1] == 3
        assert a.delta_efg.data.shape[-1] == 9

    def test_mutant_also_has_ring_contributions(self, mutant):
        assert mutant.ring_contributions is not None
        assert mutant.ring_contributions.n_pairs > 0


# ── Catalog ──────────────────────────────────────────────────────────


class TestCatalog:

    def test_all_baseline_files_registered(self):
        from pathlib import Path
        on_disk = {p.stem for p in Path(BASELINE).glob("*.npy")}
        unregistered = on_disk - set(CATALOG.keys())
        assert not unregistered, f"Unregistered: {unregistered}"

    def test_all_mutant_files_registered(self):
        from pathlib import Path
        on_disk = {p.stem for p in Path(MUTANT).glob("*.npy")}
        unregistered = on_disk - set(CATALOG.keys())
        assert not unregistered, f"Unregistered: {unregistered}"

    def test_required_count(self):
        required = [s for s in CATALOG.values() if s.required]
        assert len(required) >= 30

    def test_every_spec_has_description(self):
        for stem, spec in CATALOG.items():
            assert spec.description, f"{stem} has no description"
