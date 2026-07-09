"""SDK round-trip tests for TripeptideGroup.

Verifies that:
  - the catalog registers the tripeptide_* NPYs owned by this producer,
  - `nmr_extract.load()` attaches a TripeptideGroup when the NPYs are
    present, and
  - each field surfaces with the expected dtype / shape, including the
    NaN-fill convention on per-direction neighbour residuals.

Uses a synthetic extraction directory (tmp_path) built from the catalog's
required-shape contract — avoids depending on a real extraction host
with the tensorcs15 Postgres replica configured.
"""

import numpy as np
import pytest

from nmr_extract import (
    CATALOG,
    TripeptideGroup,
    load,
)

from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


# Synthetic atom count — small enough to inspect by hand if a test fails,
# large enough that the (N, 9) and (N, 3) shapes are unambiguous.
N_ATOMS = 32


def _required_identity_npys(out_dir, n_atoms):
    write_required_sdk_npys(out_dir, n_atoms)


def _required_calculator_npys(out_dir, n_atoms):
    write_required_sdk_npys(out_dir, n_atoms)


def _tripeptide_npys(out_dir, n_atoms):
    """Write the original seven tripeptide_* NPYs with deterministic synthetic
    values. Half the atoms are tagged as having a BB match; one
    direction (prev) is sparse so the NaN-fill convention is exercised."""
    # bb_shielding — first 16 atoms have a match (row of ones); rest NaN.
    bb = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    bb[:16] = 1.0
    np.save(out_dir / "tripeptide_bb_shielding.npy", bb)

    bb_res = np.full((n_atoms, 3), np.nan, dtype=np.float64)
    bb_res[:16] = [0.01, 0.02, 0.03]
    np.save(out_dir / "tripeptide_bb_residual_vec.npy", bb_res)

    bb_dist = np.full(n_atoms, np.nan, dtype=np.float64)
    bb_dist[:16] = 0.0374
    np.save(out_dir / "tripeptide_bb_match_distance.npy", bb_dist)

    bb_tag = np.zeros(n_atoms, dtype=np.int8)
    bb_tag[:16] = 1
    np.save(out_dir / "tripeptide_bb_method_tag.npy", bb_tag)

    # neighbor_shielding — same 16 atoms get a Δσ sum.
    n_sh = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    n_sh[:16] = 0.5
    np.save(out_dir / "tripeptide_neighbor_shielding.npy", n_sh)

    # prev residual: only atoms 0..7 carry contribution; 8..15 stay NaN.
    # next residual: atoms 0..15 all contribute.
    prev = np.full((n_atoms, 3), np.nan, dtype=np.float64)
    prev[:8] = [-0.04, -0.05, -0.06]
    np.save(out_dir / "tripeptide_neighbor_residual_vec_prev.npy", prev)

    nxt = np.full((n_atoms, 3), np.nan, dtype=np.float64)
    nxt[:16] = [0.07, 0.08, 0.09]
    np.save(out_dir / "tripeptide_neighbor_residual_vec_next.npy", nxt)


def _tripeptide_new_npys(out_dir, n_atoms):
    """Write the four pure-emit arrays added after the original group."""
    match_atoms = np.zeros((n_atoms, 5), dtype=np.float64)
    match_atoms[:, 0] = np.arange(n_atoms) // 4
    match_atoms[:16, 1] = 1
    match_atoms[:16, 2] = np.arange(1, 17)
    match_atoms[:16, 3] = 0.0374
    match_atoms[:16, 4] = 1
    np.save(out_dir / "tripeptide_bb_match_atoms.npy", match_atoms)

    prev = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    prev[:8] = 0.25
    np.save(out_dir / "tripeptide_neighbor_shielding_prev.npy", prev)

    nxt = np.full((n_atoms, 9), np.nan, dtype=np.float64)
    nxt[:16] = 0.75
    np.save(out_dir / "tripeptide_neighbor_shielding_next.npy", nxt)

    ref = np.array([[1234, 1, -120, 140, 1]], dtype=np.float64)
    np.save(out_dir / "tripeptide_neighbor_reference.npy", ref)


# ── Tests ────────────────────────────────────────────────────────────


class TestTripeptideCatalog:

    def test_specs_registered(self):
        expected = {
            "tripeptide_bb_shielding",
            "tripeptide_bb_residual_vec",
            "tripeptide_bb_match_distance",
            "tripeptide_bb_method_tag",
            "tripeptide_bb_match_atoms",
            "tripeptide_neighbor_shielding",
            "tripeptide_neighbor_shielding_prev",
            "tripeptide_neighbor_shielding_next",
            "tripeptide_neighbor_residual_vec_prev",
            "tripeptide_neighbor_residual_vec_next",
            "tripeptide_neighbor_reference",
        }
        missing = expected - set(CATALOG.keys())
        assert not missing, f"Missing tripeptide specs: {missing}"

    def test_all_specs_are_optional(self):
        for stem, spec in CATALOG.items():
            if stem.startswith("tripeptide_"):
                assert not spec.required, (
                    f"{stem} marked required; tripeptide is host-dependent "
                    "(needs tensorcs15 DSN) so must be optional"
                )


class TestTripeptideLoad:

    @pytest.fixture
    def fake_extraction(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        _tripeptide_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        return tmp_path

    @pytest.fixture
    def fake_extraction_new(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        _tripeptide_npys(tmp_path, N_ATOMS)
        _tripeptide_new_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        return tmp_path

    def test_tripeptide_group_attached(self, fake_extraction):
        p = load(fake_extraction)
        assert p.tripeptide is not None
        assert isinstance(p.tripeptide, TripeptideGroup)

    def test_bb_shielding_shape_and_values(self, fake_extraction):
        p = load(fake_extraction)
        tp = p.tripeptide
        assert tp.bb_shielding is not None
        # Wrapped as ShieldingTensor; underlying numpy is .data.
        data = tp.bb_shielding.data
        assert data.shape == (N_ATOMS, 9)
        # First 16 atoms = 1.0, remainder = NaN (matched / unmatched).
        assert np.all(data[:16] == 1.0)
        assert np.all(np.isnan(data[16:]))

    def test_bb_method_tag_dtype(self, fake_extraction):
        p = load(fake_extraction)
        tag = p.tripeptide.bb_method_tag
        # int8 per the C++ writer.
        assert tag.dtype == np.int8
        assert tag[0] == 1
        assert tag[-1] == 0

    def test_neighbor_residuals_use_nan_for_absent_direction(
            self, fake_extraction):
        """M4 contract: NaN at per-atom positions where a direction's
        contribution was absent. The synthetic fixture sets atoms 8..15
        with no i-1 contribution but a present i+1 contribution. The
        NPY load must surface NaN on prev for those atoms while keeping
        next finite."""
        p = load(fake_extraction)
        tp = p.tripeptide
        prev = tp.neighbor_residual_vec_prev.data
        nxt = tp.neighbor_residual_vec_next.data
        assert prev.shape == (N_ATOMS, 3)
        assert nxt.shape == (N_ATOMS, 3)

        # 0..7: both directions present.
        assert np.all(np.isfinite(prev[:8]))
        assert np.all(np.isfinite(nxt[:8]))

        # 8..15: only i+1 direction contributed.
        assert np.all(np.isnan(prev[8:16]))
        assert np.all(np.isfinite(nxt[8:16]))

        # 16..31: no match in either direction.
        assert np.all(np.isnan(prev[16:]))
        assert np.all(np.isnan(nxt[16:]))

    def test_original_seven_leave_new_fields_none(self, fake_extraction):
        p = load(fake_extraction)
        tp = p.tripeptide
        assert tp.bb_match_atoms is None
        assert tp.neighbor_shielding_prev is None
        assert tp.neighbor_shielding_next is None
        assert tp.neighbor_reference is None

    def test_new_pure_emit_fields_load(self, fake_extraction_new):
        p = load(fake_extraction_new)
        tp = p.tripeptide
        assert tp.bb_match_atoms.shape == (N_ATOMS, 5)
        assert tp.neighbor_shielding_prev.data.shape == (N_ATOMS, 9)
        assert tp.neighbor_shielding_next.data.shape == (N_ATOMS, 9)
        assert tp.neighbor_reference.shape == (1, 5)
        assert tp.bb_match_atoms[0, 2] == 1
        assert np.all(np.isfinite(tp.neighbor_shielding_prev.data[:8]))
        assert np.all(np.isnan(tp.neighbor_shielding_prev.data[8:]))
        np.testing.assert_array_equal(
            tp.neighbor_reference[0],
            np.array([1234, 1, -120, 140, 1], dtype=np.float64),
        )


class TestTripeptideAbsent:
    """When no tripeptide NPYs are on disk, the group stays None."""

    @pytest.fixture
    def fake_extraction_no_tripeptide(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Intentionally do NOT write any tripeptide_* NPYs.
        return tmp_path

    def test_tripeptide_is_none(self, fake_extraction_no_tripeptide):
        p = load(fake_extraction_no_tripeptide)
        assert p.tripeptide is None
