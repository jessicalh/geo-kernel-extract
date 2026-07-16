"""SDK round-trip tests for LarsenHBondGroup.

Verifies that:
  - the catalog registers the complete Larsen per-atom and pair contract,
  - all are optional (host-dependent on LarsenHBondGrid availability),
  - `nmr_extract.load()` attaches a LarsenHBondGroup when the NPYs are
    present, and
  - each field surfaces with the expected dtype / shape.

Uses a synthetic optional-calculator round trip.
"""

import numpy as np
import pytest

from nmr_extract import (
    CATALOG,
    LarsenHBondGroup,
    ShieldingTensor,
    load,
)

from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


N_ATOMS = 24


def _required_identity_npys(out_dir, n_atoms):
    write_required_sdk_npys(out_dir, n_atoms)


def _required_calculator_npys(out_dir, n_atoms):
    write_required_sdk_npys(out_dir, n_atoms)


def _larsen_hbond_npys(out_dir, n_atoms):
    """Complete Larsen NPY contract with synthetic values:
       - shielding: total, atoms 0..7 = 1.0 (received contributions),
         8..n NaN.
       - 1pHB / 2pHB / 1pHaB / 2pHaB: partial per-class breakdown.
       - diagnostic_CB: should be near-zero (asserted at C++ smoke).
       - water_term: atoms 16..19 = 2.07 (solvent-exposed amide Hs).
       - count: int32 number of contributing pairs per atom.
    """
    sh = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    sh[:8] = 1.0
    np.save(out_dir / "larsen_hbond_shielding.npy", sh)

    one_pHB = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    one_pHB[:8] = 0.5
    np.save(out_dir / "larsen_hbond_1pHB_shielding.npy", one_pHB)

    two_pHB = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    two_pHB[:8] = 0.5
    np.save(out_dir / "larsen_hbond_2pHB_shielding.npy", two_pHB)

    # H-alpha term arrays are present independently of their values.
    one_pHaB = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    np.save(out_dir / "larsen_hbond_1pHaB_shielding.npy", one_pHaB)
    two_pHaB = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    np.save(out_dir / "larsen_hbond_2pHaB_shielding.npy", two_pHaB)

    diag_cb = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    np.save(out_dir / "larsen_hbond_diagnostic_CB_shielding.npy", diag_cb)

    water = np.full(n_atoms, np.nan, dtype=np.float64)
    water[16:20] = 2.07
    np.save(out_dir / "larsen_hbond_water_term.npy", water)

    count = np.zeros(n_atoms, dtype=np.int32)
    count[:8] = 2
    np.save(out_dir / "larsen_hbond_count.npy", count)

    np.save(out_dir / "larsen_corner_imputed.npy",
            np.zeros(n_atoms, dtype=np.int8))
    imputed_count = np.zeros(n_atoms, dtype=np.int32)
    imputed_count[0] = 1
    np.save(out_dir / "larsen_imputed_pair_count.npy", imputed_count)
    sidechain_count = np.zeros(n_atoms, dtype=np.int32)
    sidechain_count[1] = 1
    np.save(out_dir / "larsen_sidechain_carbonyl_pair_count.npy",
            sidechain_count)

    pair_index = np.zeros((2, 16), dtype=np.int32)
    pair_index[:, 5] = [0, 1]
    pair_index[:, 6] = [3, 5]
    pair_geometry = np.zeros((2, 6), dtype=np.float64)
    pair_geometry[:, 3:5] = [[1, 2], [0, 0]]
    pair_geometry[:, 5] = 1
    pair_isotropic = np.zeros((2, 6), dtype=np.float64)
    pair_isotropic[0, :4] = [1, 2, 3, 4]
    pair_isotropic[0, 5] = 10
    pair_compat = np.concatenate(
        [pair_index.astype(float), pair_geometry, pair_isotropic], axis=1)
    np.save(out_dir / "larsen_hbond_pairs_index.npy", pair_index)
    np.save(out_dir / "larsen_hbond_pairs_geometry.npy", pair_geometry)
    np.save(out_dir / "larsen_hbond_pairs_isotropic.npy", pair_isotropic)
    np.save(out_dir / "larsen_hbond_pairs.npy", pair_compat)

    donor_atoms = np.zeros((n_atoms, 6), dtype=np.int32)
    donor_atoms[3] = [1, 8, 2, 0, 1, 1]
    donor_candidates = np.zeros((1, 13), dtype=np.float64)
    np.save(out_dir / "larsen_sidechain_donor_atoms.npy", donor_atoms)
    np.save(out_dir / "larsen_sidechain_donor_candidates.npy",
            donor_candidates)


# ── Tests ────────────────────────────────────────────────────────────


class TestLarsenHBondCatalog:

    def test_all_specs_registered_with_pair_axis(self):
        expected = {
            "larsen_hbond_shielding",
            "larsen_hbond_1pHB_shielding",
            "larsen_hbond_2pHB_shielding",
            "larsen_hbond_1pHaB_shielding",
            "larsen_hbond_2pHaB_shielding",
            "larsen_hbond_diagnostic_CB_shielding",
            "larsen_hbond_water_term",
            "larsen_hbond_count",
            "larsen_corner_imputed",
            "larsen_imputed_pair_count",
            "larsen_sidechain_carbonyl_pair_count",
            "larsen_hbond_pairs_index",
            "larsen_hbond_pairs_geometry",
            "larsen_hbond_pairs_isotropic",
            "larsen_hbond_pairs",
        }
        missing = expected - set(CATALOG.keys())
        assert not missing, f"Missing larsen_hbond specs: {missing}"
        for stem in {
            "larsen_hbond_pairs_index",
            "larsen_hbond_pairs_geometry",
            "larsen_hbond_pairs_isotropic",
            "larsen_hbond_pairs",
        }:
            assert CATALOG[stem].native_axis == "larsen_hbond_pair"

    def test_all_specs_are_optional(self):
        for stem, spec in CATALOG.items():
            if stem.startswith("larsen_hbond_"):
                assert not spec.required, (
                    f"{stem} marked required; LarsenHBondGrid is "
                    "host-dependent (needs the dense.h5 grids) so "
                    "must be optional"
                )

    def test_shielding_outputs_keep_their_tensor_contract(self):
        for stem in (
            "larsen_hbond_shielding",
            "larsen_hbond_1pHB_shielding",
            "larsen_hbond_2pHB_shielding",
            "larsen_hbond_1pHaB_shielding",
            "larsen_hbond_2pHaB_shielding",
            "larsen_hbond_diagnostic_CB_shielding",
        ):
            spec = CATALOG[stem]
            assert spec.wrapper is ShieldingTensor
            assert spec.parity == "even"
            assert spec.irreps == CATALOG["bs_shielding"].irreps
            assert spec.e3nn_export == CATALOG["bs_shielding"].e3nn_export
            assert "under proper rotations" in spec.transformation
            assert "no improper-transform contract" in spec.transformation


class TestLarsenHBondLoad:

    @pytest.fixture
    def fake_extraction(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        _larsen_hbond_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        return tmp_path

    def test_larsen_hbond_group_attached(self, fake_extraction):
        p = load(fake_extraction)
        assert p.larsen_hbond is not None
        assert isinstance(p.larsen_hbond, LarsenHBondGroup)

    def test_shielding_shape_and_values(self, fake_extraction):
        p = load(fake_extraction)
        lh = p.larsen_hbond
        assert lh.shielding is not None
        data = lh.shielding
        assert isinstance(data, ShieldingTensor)
        assert data.data.shape == (N_ATOMS, 9)
        # First 8 atoms = 1.0 (received contributions), rest NaN.
        assert np.all(data.data[:8] == 1.0)
        assert np.all(np.isnan(data.data[8:]))

    def test_irrep_filtered_consumer_receives_all_six_tensors(
            self, fake_extraction):
        expected_fields = {
            "larsen_hbond_shielding": "shielding",
            "larsen_hbond_1pHB_shielding": "pHB_1",
            "larsen_hbond_2pHB_shielding": "pHB_2",
            "larsen_hbond_1pHaB_shielding": "pHaB_1",
            "larsen_hbond_2pHaB_shielding": "pHaB_2",
            "larsen_hbond_diagnostic_CB_shielding": "diagnostic_CB",
        }
        shielding_irreps = CATALOG["bs_shielding"].irreps
        selected = {
            stem for stem, spec in CATALOG.items()
            if spec.group == "larsen_hbond"
            and spec.irreps == shielding_irreps
        }
        assert selected == set(expected_fields)

        group = load(fake_extraction).larsen_hbond
        received = {
            stem: getattr(group, field_name)
            for stem, field_name in expected_fields.items()
            if stem in selected
        }
        assert set(received) == selected
        for datum in received.values():
            assert isinstance(datum, ShieldingTensor)
            assert datum.data.shape == (N_ATOMS, 9)

    def test_per_class_breakdown(self, fake_extraction):
        p = load(fake_extraction)
        lh = p.larsen_hbond
        assert lh.pHB_1 is not None
        assert lh.pHB_2 is not None
        # 1pHB + 2pHB should sum to total in the contributing region.
        a = lh.pHB_1.data[:8]
        b = lh.pHB_2.data[:8]
        total = lh.shielding.data[:8]
        np.testing.assert_allclose(a + b, total, rtol=1e-12)

    def test_water_term_isolated_to_solvent_exposed(self, fake_extraction):
        p = load(fake_extraction)
        water = p.larsen_hbond.water_term
        assert water is not None
        # Synthetic fixture: atoms 16..19 are unbound, rest NaN.
        assert np.all(water[16:20] == 2.07)
        assert np.all(np.isnan(water[:16]))
        assert np.all(np.isnan(water[20:]))

    def test_count_is_int32(self, fake_extraction):
        p = load(fake_extraction)
        count = p.larsen_hbond.count
        assert count.dtype == np.int32
        # First 8 atoms = 2 pairs each, rest 0.
        assert np.all(count[:8] == 2)
        assert np.all(count[8:] == 0)

    def test_pair_wrapper_and_audit_counts(self, fake_extraction):
        p = load(fake_extraction)
        pairs = p.larsen_hbond.pairs
        assert pairs.index.shape == (2, 16)
        assert pairs.geometry.shape == (2, 6)
        assert pairs.isotropic.shape == (2, 6)
        assert pairs.compatibility.shape == (2, 28)
        assert pairs.geometry[0, 4] == 2
        assert pairs.isotropic[0, 5] == pairs.isotropic[0, :4].sum()
        assert pairs.index[0, 6] == 3  # success
        assert pairs.index[1, 6] == 5  # carboxylate symmetry-filtered
        assert p.larsen_hbond.imputed_pair_count[0] == 1
        assert p.larsen_hbond.sidechain_carbonyl_pair_count[1] == 1
        assert p.larsen_sidechain_donor_audit.atoms.shape == (N_ATOMS, 6)
        assert p.larsen_sidechain_donor_audit.candidates.shape == (1, 13)

    def test_absent_npys_means_no_group(self, tmp_path):
        """If NO larsen_hbond_* NPY is present, the group is None."""
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # No larsen_hbond_* npys written.
        p = load(tmp_path)
        assert p.larsen_hbond is None
