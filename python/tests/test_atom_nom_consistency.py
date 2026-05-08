"""Deep-consistency tests for atom_nom.tbl and the structured NPY it feeds.

Two layers:

1. atom_nom.tbl deep parse — reimplements parsing in Python (independent
   of the C++ parser) and validates expected coverage, atom presence,
   stereo labels. Catches file truncation, format regressions, missing
   residues, missing prochiral stereo entries before any fleet run sees
   them.

2. atoms_category_info.npy round-trip — if an extraction is available,
   loads the structured NPY with numpy and asserts every field is
   readable, dtypes match the C++ writer, atom_index is contiguous,
   and string columns are decodable.

Failure here means the projection pipeline is broken in a way that would
otherwise surface mid-fleet-run, much later, much more expensively.
"""

from pathlib import Path
import os
import re

import numpy as np
import pytest


# ── Inputs ────────────────────────────────────────────────────────────

ATOM_NOM_TBL = Path(
    "/shared/2026Thesis/nmr-shielding/references/bmrb_data/atom_nom.tbl"
)

# Optional: if a recent extraction has atoms_category_info.npy, point at
# the directory here and the round-trip tests will run against it. The
# environment override lets CI swap in any extraction without code edits.
CATEGORY_INFO_NPY = os.environ.get(
    "NMR_CATEGORY_INFO_NPY",
    "/shared/2026Thesis/nmr-shielding/tests/data/sdk_geo_only/"
    "1Q8K/1Q8K_10023_01_boltzmann_minimum_8.86e-02_frame001/"
    "atoms_category_info.npy",
)


# ── atom_nom.tbl pure-Python parser ───────────────────────────────────

def _parse_atom_nom(path: Path) -> tuple[dict, dict, list[str]]:
    """Parse atom_nom.tbl into (per_residue, cap_atoms, raw_rows).

    per_residue: {(one_letter, atom_name): {'bmrb', 'sc'}}
    cap_atoms:   {atom_name: {'bmrb', 'sc'}}
    raw_rows:    list of stripped non-comment lines (for diagnostic use)
    """
    per_residue: dict[tuple[str, str], dict] = {}
    cap_atoms: dict[str, dict] = {}
    raw: list[str] = []

    with path.open() as f:
        for line in f:
            stripped = line.split("#", 1)[0].rstrip("\n")
            if not stripped.strip():
                continue
            raw.append(stripped)

            # Tab-split exactly mirrors the C++ parser's column logic.
            # The file uses tabs as the canonical separator; cells are
            # nominally indexed 0..10:
            #   [0] one-letter code (A..V or X)
            #   [1] empty (legacy padding column)
            #   [2] BMRB / IUPAC atom name
            #   [3] SC stereo (pro-R / pro-S / Z / E / "")
            #   [4..10] UCSF / XPLOR / MSI / PDB / SYBYL / MIDAS
            cols = stripped.split("\t")
            # Whitespace-split fallback if no tabs (shouldn't happen but
            # be defensive).
            if len(cols) < 3:
                cols = stripped.split()
                if len(cols) < 2:
                    continue
                # Whitespace-only format: collapse to (res, bmrb, sc?, ...)
                res = cols[0]
                bmrb = cols[1] if len(cols) > 1 else ""
                sc = cols[2] if len(cols) > 2 else ""
            else:
                res = cols[0].strip()
                bmrb = cols[2].strip() if len(cols) >= 3 else ""
                sc = cols[3].strip() if len(cols) >= 4 else ""

            if not bmrb:
                continue
            if len(res) != 1:
                continue

            row = {"bmrb": bmrb, "sc": sc}
            if res == "X":
                cap_atoms[bmrb] = row
            elif res.isupper():
                per_residue[(res, bmrb)] = row

    return per_residue, cap_atoms, raw


# ── 20 single-letter residue codes (canonical IUPAC) ──────────────────

STANDARD_20 = list("ARNDCQEGHILKMFPSTWYV")


# ── Layer 1: atom_nom.tbl deep consistency ────────────────────────────

@pytest.fixture(scope="module")
def atom_nom():
    if not ATOM_NOM_TBL.exists():
        pytest.skip(f"atom_nom.tbl not found at {ATOM_NOM_TBL}")
    return _parse_atom_nom(ATOM_NOM_TBL)


def test_file_exists():
    assert ATOM_NOM_TBL.exists(), (
        f"atom_nom.tbl missing at {ATOM_NOM_TBL}. The CategoryInfoProjection "
        "Configure step will fail loud at startup if this is unset; "
        "make sure the file is committed and ~/.nmr_tools.toml has "
        "bmrb_atom_nom pointing at it."
    )


def test_minimum_row_counts(atom_nom):
    """Floor on row counts. Catches truncation or parser regression."""
    per_residue, cap_atoms, raw = atom_nom
    assert len(per_residue) >= 280, (
        f"per_residue rows = {len(per_residue)}; expected >= 280. "
        "Fewer than this indicates a truncated table or a parser bug."
    )
    assert len(cap_atoms) >= 4, (
        f"cap_atoms rows = {len(cap_atoms)}; expected >= 4 "
        "(at minimum H1, H2, H3, OXT)."
    )


def test_all_standard_20_present(atom_nom):
    """Every standard 1-letter residue must have entries."""
    per_residue, _, _ = atom_nom
    seen_codes = {res for (res, _) in per_residue}
    missing = set(STANDARD_20) - seen_codes
    assert not missing, (
        f"atom_nom.tbl is missing rows for residue codes: {sorted(missing)}. "
        "Every standard-20 residue must have at least one row."
    )


def test_backbone_atoms_present_for_all_residues(atom_nom):
    """Every standard residue must have backbone N, CA, C, O.

    Backbone amide H is NOT required (Pro has none) — checked
    separately.
    """
    per_residue, _, _ = atom_nom
    required = ["N", "CA", "C", "O"]
    missing: list[str] = []
    for code in STANDARD_20:
        for atom in required:
            if (code, atom) not in per_residue:
                missing.append(f"{code}:{atom}")
    assert not missing, (
        f"Missing backbone atoms in atom_nom.tbl: {missing[:10]}"
        + (" ..." if len(missing) > 10 else "")
    )


def test_backbone_amide_h_present_except_pro(atom_nom):
    """Every non-Pro residue must have a backbone amide H entry."""
    per_residue, _, _ = atom_nom
    for code in STANDARD_20:
        if code == "P":
            continue
        assert (code, "H") in per_residue, (
            f"Residue {code}: backbone amide H missing from atom_nom.tbl"
        )


def test_cap_atoms_present_under_bmrb_names(atom_nom):
    """Cap atoms appear in atom_nom.tbl X-rows under their BMRB names.

    Note: atom_nom.tbl uses BMRB conventions for cap atoms — `O''` for the
    C-terminal extra oxygen (AMBER calls it OXT), `H''` for the C-terminal
    proton (AMBER HXT), and H1/H2/H3 for the N-terminal protons (AMBER
    matches BMRB here). The C++ projection layer must therefore carry an
    AMBER→BMRB cap-name map for OXT/HXT — atom_nom.tbl alone does NOT
    contain those AMBER strings as keys.
    """
    _, cap_atoms, _ = atom_nom
    for atom in ("H1", "H2", "H3", "O''", "H''"):
        assert atom in cap_atoms, (
            f"BMRB cap atom {atom!r} missing from atom_nom.tbl X-rows"
        )


def test_amber_cap_names_NOT_in_table(atom_nom):
    """AMBER OXT/HXT are intentionally absent from atom_nom.tbl as keys.

    This is the structural fact that motivates the AMBER→BMRB cap-name
    map in CategoryInfoProjection. If a future atom_nom.tbl revision adds
    OXT/HXT directly, the projection should prefer the table entry; this
    test will fail and remind us to revisit the C++ map.
    """
    _, cap_atoms, _ = atom_nom
    for amber_only in ("OXT", "HXT"):
        assert amber_only not in cap_atoms, (
            f"atom_nom.tbl now contains {amber_only!r} as a BMRB name. "
            "Revisit the AMBER→BMRB cap-name map in "
            "CategoryInfoProjection.cpp — the projection may now be "
            "able to use the table directly."
        )


def test_methylene_stereo_pair_consistency(atom_nom):
    """For every prochiral HB2/HB3 pair, one is pro-R and the other pro-S.

    The mapping is NOT uniform across residues — CIP priority depends on
    side-chain γ-atom chemistry, so HB2 is sometimes pro-R, sometimes
    pro-S. What IS invariant is that each pair carries opposite labels.

    Excludes residues whose Cβ carries three equivalent methyl Hs (ALA's
    HB1/HB2/HB3) — those are NOT prochiral methylenes; SC is correctly
    empty there. The discriminator is the presence of HB1 alongside HB2/
    HB3.

    Catches a malformed row that has only one of the prochiral pair
    labeled, or both blank, or both the same.
    """
    per_residue, _, _ = atom_nom
    bad: list[str] = []
    for res in STANDARD_20:
        if (res, "HB2") not in per_residue or (res, "HB3") not in per_residue:
            continue
        # Skip methyl groups (three Hs on Cβ): ALA-style HB1+HB2+HB3.
        if (res, "HB1") in per_residue:
            continue
        sc2 = per_residue[(res, "HB2")]["sc"]
        sc3 = per_residue[(res, "HB3")]["sc"]
        valid = {("pro-R", "pro-S"), ("pro-S", "pro-R")}
        if (sc2, sc3) not in valid:
            bad.append(f"{res}: HB2={sc2!r} HB3={sc3!r}")
    assert not bad, (
        "Prochiral HB2/HB3 pairs must be pro-R/pro-S (in some order). "
        "Bad pairs: " + "; ".join(bad)
    )


def test_ala_methyl_has_empty_stereo(atom_nom):
    """ALA HB1/HB2/HB3 are equivalent methyl Hs, not prochiral.

    SC column should be empty for all three. Catches a regression where
    someone "fixed" the empty SC assuming it was a missing label.
    """
    per_residue, _, _ = atom_nom
    if (("A", "HB1") not in per_residue
            or ("A", "HB2") not in per_residue
            or ("A", "HB3") not in per_residue):
        pytest.skip("ALA HB1/HB2/HB3 not all present in atom_nom.tbl")
    for atom in ("HB1", "HB2", "HB3"):
        sc = per_residue[("A", atom)]["sc"]
        assert sc == "", (
            f"ALA {atom} has SC={sc!r}; expected empty (methyl Hs are "
            "not prochiral)."
        )


# Per-residue CIP-priority-driven stereo assignments (verified against
# atom_nom.tbl as it sits today). The list is the chemistry, not a
# convention; if it changes, atom_nom.tbl has changed under us.
EXPECTED_HB_STEREO = {
    # HB2 = pro-R
    "R": ("pro-R", "pro-S"),
    "E": ("pro-R", "pro-S"),
    "Q": ("pro-R", "pro-S"),
    "L": ("pro-R", "pro-S"),
    "K": ("pro-R", "pro-S"),
    "F": ("pro-R", "pro-S"),
    "P": ("pro-R", "pro-S"),
    "W": ("pro-R", "pro-S"),
    "Y": ("pro-R", "pro-S"),
    # HB2 = pro-S
    "D": ("pro-S", "pro-R"),
    "N": ("pro-S", "pro-R"),
    "C": ("pro-S", "pro-R"),
    "H": ("pro-S", "pro-R"),
    "M": ("pro-S", "pro-R"),
    "S": ("pro-S", "pro-R"),
}


def test_methylene_stereo_pinned_per_residue(atom_nom):
    """Pin the per-residue HB2/HB3 stereo against atom_nom.tbl as it stands.

    Drift here means atom_nom.tbl was edited or replaced — likely
    intentional, but we want the test to fail loudly so we can audit.
    """
    per_residue, _, _ = atom_nom
    bad: list[str] = []
    for res, (sc2_exp, sc3_exp) in EXPECTED_HB_STEREO.items():
        if (res, "HB2") not in per_residue or (res, "HB3") not in per_residue:
            continue
        sc2 = per_residue[(res, "HB2")]["sc"]
        sc3 = per_residue[(res, "HB3")]["sc"]
        if (sc2, sc3) != (sc2_exp, sc3_exp):
            bad.append(f"{res}: HB2={sc2} HB3={sc3} (expected {sc2_exp}/{sc3_exp})")
    assert not bad, "Per-residue HB stereo drift: " + "; ".join(bad)


def test_glycine_ha_pair_stereo(atom_nom):
    """Gly HA2/HA3 follows the same suffix-2 = pro-R convention."""
    per_residue, _, _ = atom_nom
    assert (("G", "HA2") in per_residue), "Gly HA2 missing from atom_nom.tbl"
    assert (("G", "HA3") in per_residue), "Gly HA3 missing from atom_nom.tbl"
    assert per_residue[("G", "HA2")]["sc"] == "pro-R", \
        "Gly HA2 should be pro-R per atom_nom.tbl"
    assert per_residue[("G", "HA3")]["sc"] == "pro-S", \
        "Gly HA3 should be pro-S per atom_nom.tbl"


def test_arg_guanidinium_z_e_labels(atom_nom):
    """Arg guanidinium 1H1/1H2 and 2H1/2H2 carry Z / E labels."""
    per_residue, _, _ = atom_nom
    expected = [
        ("R", "HH11", "Z"),
        ("R", "HH12", "E"),
        ("R", "HH21", "Z"),
        ("R", "HH22", "E"),
    ]
    for res, atom, sc_expected in expected:
        if (res, atom) not in per_residue:
            continue  # tolerate naming variation; just don't check absent
        sc = per_residue[(res, atom)]["sc"]
        assert sc == sc_expected, (
            f"{res}:{atom} stereo = '{sc}'; expected '{sc_expected}'"
        )


def test_no_stereo_label_is_garbage(atom_nom):
    """SC column values are confined to a small known set."""
    per_residue, cap_atoms, _ = atom_nom
    valid = {"", "pro-R", "pro-S", "Z", "E"}
    bad: list[str] = []
    for (res, atom), row in per_residue.items():
        if row["sc"] not in valid:
            bad.append(f"{res}:{atom}={row['sc']!r}")
    for atom, row in cap_atoms.items():
        if row["sc"] not in valid:
            bad.append(f"X:{atom}={row['sc']!r}")
    assert not bad, (
        f"Unexpected SC labels in atom_nom.tbl: {bad[:20]}"
        + (" ..." if len(bad) > 20 else "")
    )


def test_atom_names_fit_s8(atom_nom):
    """No atom name exceeds 7 chars (S8 column width less null margin)."""
    per_residue, cap_atoms, _ = atom_nom
    too_long: list[str] = []
    for (res, atom) in per_residue:
        if len(atom) > 7:
            too_long.append(f"{res}:{atom}")
    for atom in cap_atoms:
        if len(atom) > 7:
            too_long.append(f"X:{atom}")
    assert not too_long, (
        f"Atom names exceed S8 width 7: {too_long}. "
        "The C++ writer would abort at runtime; widen the column or "
        "rename the offending entries."
    )


# ── Layer 2: atoms_category_info.npy round-trip ───────────────────────

@pytest.fixture(scope="module")
def npy_path():
    p = Path(CATEGORY_INFO_NPY)
    if not p.exists():
        pytest.skip(
            f"atoms_category_info.npy not found at {p}. "
            "Run nmr_extract --pdb on a fixture to generate one, or set "
            "NMR_CATEGORY_INFO_NPY to an existing extraction."
        )
    return p


@pytest.fixture(scope="module")
def npy_array(npy_path):
    return np.load(npy_path)


def test_structured_dtype(npy_array):
    """numpy can load the file with a structured dtype, no errors."""
    assert npy_array.dtype.fields is not None, (
        "Loaded array is not a structured array. The NPY header's "
        "'descr' field probably did not produce a list-of-tuples dtype."
    )


EXPECTED_FIELDS = {
    "atom_index": "<i4",
    "residue_index": "<i4",
    "element": "i1",
    "amber_atom_name": "|S8",
    "iupac_atom_name": "|S8",
    "bmrb_atom_name": "|S8",
    "amber_residue_3letter": "|S4",
    "iupac_residue_3letter": "|S4",
    "bmrb_residue_3letter": "|S4",
    "residue_1letter": "|S1",
    "residue_type": "i1",
    "residue_variant_index": "i1",
    "terminal_state": "i1",
    "locant": "i1",
    "branch_outer": "i1",
    "branch_inner": "i1",
    "di_index": "i1",
    "backbone_role": "i1",
    "prochiral": "i1",
    "planar_group": "i1",
    "planar_stereo": "i1",
    "polar_h_kind": "i1",
    "ring_position_primary": "i1",
    "ring_position_secondary": "i1",
    "pseudoatom_kind": "i1",
    "in_super_group": "i1",
    "aromatic": "i1",
    "formal_charge": "i1",
    "is_exchangeable": "i1",
    "iupac_naming_provenance": "i1",
    "bmrb_naming_provenance": "i1",
}


def test_all_fields_present(npy_array):
    """Every expected field name appears in the structured dtype."""
    fields = set(npy_array.dtype.names or ())
    missing = set(EXPECTED_FIELDS) - fields
    extra = fields - set(EXPECTED_FIELDS)
    assert not missing, f"Missing fields in NPY: {sorted(missing)}"
    # `extra` is informational; we don't fail on additions but log them.
    if extra:
        print(f"NPY has extra fields not in expected schema: {sorted(extra)}")


def test_field_dtypes_match(npy_array):
    """Each field's stored dtype matches the C++ writer's declaration.

    numpy normalises some byte-order codes; we tolerate `i1` vs `|i1`
    and equivalent S<N> spellings.
    """
    bad: list[str] = []
    for name, expected in EXPECTED_FIELDS.items():
        if name not in npy_array.dtype.fields:
            continue  # caught by test_all_fields_present
        actual = npy_array.dtype.fields[name][0].str
        # Normalise: 'i1' and '|i1' are equivalent; '<i4' and '=i4' on
        # little-endian hosts are equivalent.
        norm_a = actual.lstrip("|").lstrip("=")
        norm_e = expected.lstrip("|").lstrip("=")
        if norm_a != norm_e and not (norm_a == "i4" and norm_e == "i4"):
            bad.append(f"{name}: actual={actual!r} expected={expected!r}")
    assert not bad, "Field dtype mismatches: " + "; ".join(bad)


def test_atom_index_is_contiguous(npy_array):
    """atom_index is 0..N-1, in order, no gaps."""
    n = len(npy_array)
    assert (npy_array["atom_index"] == np.arange(n, dtype=np.int32)).all(), (
        "atom_index column is not 0..N-1 in row order. "
        "Either records were written out of order, or the column was "
        "not populated."
    )


def test_residue_index_in_range(npy_array):
    """residue_index values are in [0, max+1) and non-decreasing."""
    ri = npy_array["residue_index"]
    assert ri.min() == 0, f"residue_index min = {ri.min()} (expected 0)"
    # Non-decreasing because atoms within a residue are contiguous.
    diff = np.diff(ri.astype(np.int64))
    assert (diff >= 0).all(), (
        "residue_index decreases between adjacent atoms — atoms within "
        "the same residue are not contiguous."
    )


def test_amber_atom_name_nonempty(npy_array):
    """Every atom has a non-empty AMBER atom name."""
    names = npy_array["amber_atom_name"]
    empty = (names == b"").sum()
    assert empty == 0, f"{empty} atoms have empty amber_atom_name"


def test_residue_three_letters_nonempty(npy_array):
    """Every atom has all three residue 3-letter codes populated."""
    for col in ("amber_residue_3letter",
                "iupac_residue_3letter",
                "bmrb_residue_3letter"):
        empty = (npy_array[col] == b"").sum()
        assert empty == 0, f"{empty} atoms have empty {col}"


def test_provenance_is_zero_or_one(npy_array):
    """Provenance enum values are confined to {0=Match, 1=MissLogged}."""
    for col in ("iupac_naming_provenance", "bmrb_naming_provenance"):
        vals = set(np.unique(npy_array[col]).tolist())
        bad = vals - {0, 1}
        assert not bad, f"{col} has unexpected values: {sorted(bad)}"


def test_residue_one_letter_is_a_letter(npy_array):
    """residue_1letter is an uppercase ASCII letter."""
    letters = np.unique(npy_array["residue_1letter"])
    valid = re.compile(rb"^[A-Z]$")
    bad = [b for b in letters.tolist() if not valid.match(b)]
    assert not bad, f"residue_1letter has non-letter values: {bad}"


def test_strings_are_decodable(npy_array):
    """Every S-typed column decodes cleanly to str (no embedded garbage)."""
    string_cols = [
        "amber_atom_name", "iupac_atom_name", "bmrb_atom_name",
        "amber_residue_3letter", "iupac_residue_3letter",
        "bmrb_residue_3letter", "residue_1letter",
    ]
    bad: list[str] = []
    for col in string_cols:
        for i, b in enumerate(npy_array[col]):
            try:
                s = b.decode("ascii")
            except UnicodeDecodeError:
                bad.append(f"{col}[{i}] non-ASCII: {b!r}")
                continue
            # Must be printable.
            if any(ord(c) < 0x20 or ord(c) > 0x7e for c in s):
                bad.append(f"{col}[{i}] non-printable: {b!r}")
    assert not bad, "String columns have garbage bytes: " + "; ".join(bad[:20])
