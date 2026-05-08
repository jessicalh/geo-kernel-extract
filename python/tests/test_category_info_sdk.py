"""SDK round-trip on atoms_category_info.npy.

Loads a real mutant-pair extraction directory, walks the CategoryInfo
wrapper, and asserts the convenience accessors agree with the legacy
flat NPYs (element, residue_type, residue_index). Also exercises the
DeltaGroup dia/para component fields land cleanly via the loader.

Set ``NMR_MUTANT_SMOKE_DIR`` to point at an existing mutant extraction
to run; otherwise the tests skip. The default path is the one our
session uses for /tmp/mutant_smoke.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from nmr_extract import load, CategoryInfo


SMOKE_DIR = Path(os.environ.get("NMR_MUTANT_SMOKE_DIR", "/tmp/mutant_smoke"))


@pytest.fixture(scope="module")
def protein():
    if not SMOKE_DIR.exists():
        pytest.skip(
            f"Mutant smoke extraction not found at {SMOKE_DIR}. Run "
            "build/nmr_extract --mutant --wt tests/data/orca/A0A7C5FAR6_WT "
            "--ala tests/data/orca/A0A7C5FAR6_ALA --output /tmp/mutant_smoke "
            "or set NMR_MUTANT_SMOKE_DIR.")
    return load(SMOKE_DIR)


# ── CategoryInfo wrapper ────────────────────────────────────────────


def test_category_info_loads(protein):
    """CategoryInfo wrapper is non-None and reports n_atoms."""
    assert protein.category_info is not None
    assert isinstance(protein.category_info, CategoryInfo)
    assert protein.category_info.n_atoms == protein.n_atoms


def test_element_matches_legacy_npy(protein):
    """The structured-NPY element column equals element.npy atomic numbers."""
    ci = protein.category_info
    assert (ci.element.astype(np.int32) == protein.element).all(), (
        "CategoryInfo.element should be atomic number, matching the "
        "legacy element.npy convention. If this fails, the C++ writer "
        "is emitting Element-enum index instead of AtomicNumberForElement."
    )


def test_residue_index_matches_legacy_npy(protein):
    """The structured-NPY residue_index column matches residue_index.npy."""
    ci = protein.category_info
    assert (ci.residue_index.astype(np.int32) == protein.residue_index).all()


def test_atom_index_is_contiguous(protein):
    """atom_index column is 0..N-1."""
    ci = protein.category_info
    assert (ci.atom_index == np.arange(ci.n_atoms, dtype=ci.atom_index.dtype)).all()


def test_decoded_strings_nonempty(protein):
    """All atom-name and residue-3letter columns decode to non-empty bytes."""
    ci = protein.category_info
    assert (ci.amber_atom_name != b"").all()
    assert (ci.iupac_residue_3letter != b"").all()


def test_stratification_masks(protein):
    """Convenience boolean masks add up consistently with raw element column."""
    ci = protein.category_info
    # Atom counts by element should match atomic-number unique count.
    elements, counts = np.unique(ci.element, return_counts=True)
    # H + C + N + O + S = total
    assert sum(counts) == ci.n_atoms

    # Backbone count is non-zero and at most ~7 atoms × n_residues.
    bb = int(ci.is_backbone.sum())
    assert bb > 0
    n_res = len(np.unique(ci.residue_index))
    assert bb <= 7 * n_res, (
        f"backbone count {bb} > 7 × {n_res} residues — implausible.")

    # Aromatic count is non-zero (this fixture has PHE/TYR/HIS/TRP).
    assert int(ci.is_aromatic.sum()) > 0


def test_phe_aromatic_carbon_stratification(protein):
    """PHE has 6 aromatic carbons (Cγ, Cδ1, Cδ2, Cε1, Cε2, Cζ) per residue."""
    ci = protein.category_info
    n_phe_residues = len(
        np.unique(ci.residue_index[ci.iupac_residue_3letter == b"PHE"]))
    if n_phe_residues == 0:
        pytest.skip("Fixture has no PHE residues")
    mask = (
        ci.is_aromatic
        & (ci.element == 6)
        & (ci.iupac_residue_3letter == b"PHE")
    )
    assert int(mask.sum()) == 6 * n_phe_residues


# ── DeltaGroup dia/para fields ─────────────────────────────────────


def test_delta_dia_para_present(protein):
    """All six dia/para shielding NPYs are wrapped on DeltaGroup."""
    if protein.delta is None:
        pytest.skip("Fixture is not a mutation pair")
    d = protein.delta
    for f in (
        "wt_shielding_diamagnetic",
        "wt_shielding_paramagnetic",
        "mut_shielding_diamagnetic",
        "mut_shielding_paramagnetic",
        "delta_shielding_diamagnetic",
        "delta_shielding_paramagnetic",
    ):
        assert getattr(d, f) is not None, f"DeltaGroup.{f} missing"
        assert getattr(d, f).data.shape == (protein.n_atoms, 9)


def test_delta_components_sum_to_total_within_orca_precision(protein):
    """delta_dia + delta_para should equal delta_total at ORCA's output precision.

    ORCA renders total / dia / para as separate matrices with limited
    decimal places, so the sum identity holds at ~1e-3 ppm rather than
    machine epsilon. This test pins that expected discrepancy.
    """
    if protein.delta is None or protein.delta.delta_shielding_diamagnetic is None:
        pytest.skip("No mutant delta dia/para data")
    d = protein.delta
    total = d.shielding.data
    dia = d.delta_shielding_diamagnetic.data
    para = d.delta_shielding_paramagnetic.data
    diff = np.abs(total - (dia + para))
    assert diff.max() < 1e-2, (
        f"delta_total - (delta_dia + delta_para) max = {diff.max():.3e}; "
        f"expected ~ORCA output precision (1e-3 ppm). Larger discrepancy "
        f"indicates the components were not sourced from the same ORCA run "
        f"or the sign convention diverged.")


def test_sides_minus_delta_consistent_on_matched(protein):
    """On bound atoms, (wt - mut) should equal the delta tensor.

    Unmatched atoms have all-zero deltas and zero on the unfilled side,
    so the identity only holds where the atom was actually matched.
    """
    if protein.delta is None or protein.delta.delta_shielding_diamagnetic is None:
        pytest.skip("No mutant delta dia/para data")
    d = protein.delta
    # Matched atoms have nonzero values in delta_shielding (modulo the
    # zero-tensor case where chemistry happens to be identical).
    delta_dia = d.delta_shielding_diamagnetic.data
    matched = np.abs(delta_dia).sum(axis=1) > 1e-12
    if matched.sum() == 0:
        pytest.skip("No bound atoms in this fixture")
    wt = d.wt_shielding_diamagnetic.data[matched]
    mut = d.mut_shielding_diamagnetic.data[matched]
    diff = wt - mut
    assert np.allclose(diff, delta_dia[matched]), (
        "On matched atoms, wt_dia - mut_dia must equal delta_dia.")
