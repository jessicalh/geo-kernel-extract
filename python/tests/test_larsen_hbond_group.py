"""SDK round-trip tests for LarsenHBondGroup.

Verifies that:
  - the catalog registers all 8 larsen_hbond_* NPYs,
  - all are optional (host-dependent on LarsenHBondGrid availability),
  - `nmr_extract.load()` attaches a LarsenHBondGroup when the NPYs are
    present, and
  - each field surfaces with the expected dtype / shape.

Mirrors test_tripeptide_group.py.
"""

import numpy as np
import pytest

from nmr_extract import (
    CATALOG,
    LarsenHBondGroup,
    load,
)

from _topology_fixture import write_minimal_topology_sidecar


N_ATOMS = 24


def _required_identity_npys(out_dir, n_atoms):
    pos = np.zeros((n_atoms, 3), dtype=np.float64)
    np.save(out_dir / "pos.npy", pos)
    np.save(out_dir / "element.npy", np.full(n_atoms, 6, dtype=np.int32))
    np.save(out_dir / "residue_index.npy",
            np.arange(n_atoms, dtype=np.int32) // 4)
    np.save(out_dir / "residue_type.npy",
            np.zeros(n_atoms, dtype=np.int32))
    np.save(out_dir / "ring_contributions.npy",
            np.zeros((0, 59), dtype=np.float64))
    np.save(out_dir / "ring_geometry.npy",
            np.zeros((0, 10), dtype=np.float64))


def _required_calculator_npys(out_dir, n_atoms):
    stems_9 = [
        "bs_shielding", "hm_shielding", "pq_shielding",
        "disp_shielding", "ringchi_shielding",
        "mc_shielding", "hbond_shielding",
    ]
    for s in stems_9:
        np.save(out_dir / f"{s}.npy",
                np.zeros((n_atoms, 9), dtype=np.float64))

    for stem in ("bs", "hm", "pq", "disp"):
        np.save(out_dir / f"{stem}_per_type_T0.npy",
                np.zeros((n_atoms, 8), dtype=np.float64))
        np.save(out_dir / f"{stem}_per_type_T2.npy",
                np.zeros((n_atoms, 40), dtype=np.float64))

    np.save(out_dir / "bs_total_B.npy",
            np.zeros((n_atoms, 3), dtype=np.float64))
    np.save(out_dir / "bs_ring_counts.npy",
            np.zeros((n_atoms, 4), dtype=np.float64))
    np.save(out_dir / "mc_category_T2.npy",
            np.zeros((n_atoms, 25), dtype=np.float64))
    np.save(out_dir / "mc_scalars.npy",
            np.zeros((n_atoms, 6), dtype=np.float64))
    np.save(out_dir / "hbond_scalars.npy",
            np.zeros((n_atoms, 3), dtype=np.float64))
    np.save(out_dir / "dssp_backbone.npy",
            np.zeros((n_atoms, 5), dtype=np.float64))
    np.save(out_dir / "aimnet2_charges.npy",
            np.zeros(n_atoms, dtype=np.float64))
    np.save(out_dir / "aimnet2_aim.npy",
            np.zeros((n_atoms, 256), dtype=np.float32))
    np.save(out_dir / "aimnet2_efg.npy",
            np.zeros((n_atoms, 9), dtype=np.float64))
    np.save(out_dir / "aimnet2_efg_aromatic.npy",
            np.zeros((n_atoms, 9), dtype=np.float64))
    np.save(out_dir / "aimnet2_efg_backbone.npy",
            np.zeros((n_atoms, 9), dtype=np.float64))


def _larsen_hbond_npys(out_dir, n_atoms):
    """8 NPYs with synthetic values:
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

    # Phase 2 stubs (zero — Phase 1 doesn't populate these).
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


# ── Tests ────────────────────────────────────────────────────────────


class TestLarsenHBondCatalog:

    def test_all_eight_specs_registered(self):
        expected = {
            "larsen_hbond_shielding",
            "larsen_hbond_1pHB_shielding",
            "larsen_hbond_2pHB_shielding",
            "larsen_hbond_1pHaB_shielding",
            "larsen_hbond_2pHaB_shielding",
            "larsen_hbond_diagnostic_CB_shielding",
            "larsen_hbond_water_term",
            "larsen_hbond_count",
        }
        missing = expected - set(CATALOG.keys())
        assert not missing, f"Missing larsen_hbond specs: {missing}"

    def test_all_specs_are_optional(self):
        for stem, spec in CATALOG.items():
            if stem.startswith("larsen_hbond_"):
                assert not spec.required, (
                    f"{stem} marked required; LarsenHBondGrid is "
                    "host-dependent (needs the dense.h5 grids) so "
                    "must be optional"
                )


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
        data = lh.shielding.data
        assert data.shape == (N_ATOMS, 9)
        # First 8 atoms = 1.0 (received contributions), rest NaN.
        assert np.all(data[:8] == 1.0)
        assert np.all(np.isnan(data[8:]))

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

    def test_absent_npys_means_no_group(self, tmp_path):
        """If NO larsen_hbond_* NPY is present, the group is None."""
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # No larsen_hbond_* npys written.
        p = load(tmp_path)
        assert p.larsen_hbond is None
