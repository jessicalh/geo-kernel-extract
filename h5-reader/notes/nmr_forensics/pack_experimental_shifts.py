#!/usr/bin/env python3
"""
pack_experimental_shifts.py — bind BMRB/RefDB shifts to H5 atom_index
and produce per-protein experimental_shifts/ bundles.

Scope (2026-04-20): 10 calibration proteins. Fleet deferred
(memory:feedback_nmr_scope_10_first).

Design: no string-translation logic. Four typed data constants cover
every known divergence between BMRB/RefDB (IUPAC) and the H5 `/atoms`
naming (IUPAC for most, CHARMM-in-H5 at six gap positions, one
library bug on ALA β-methyl).

  - H5_NAME_FROM_IUPAC          — (residue, IUPAC atom) → H5 atom name.
                                  Covers six registry-gap positions where
                                  H5 carries CHARMM spellings.
  - BMRB_METHYL_VIA_PARENT      — ALA β-methyl: bind 3 BMRB rows to the
                                  3-atom set via /topology/parent_atom_index.
                                  Handles the library wildcard-β bug.
  - REFDB_PSEUDO_PARENTS        — (residue, pseudo name) → (parent, count).
                                  RefDB 2.1 methyl pseudos bound via parent.
  - H5_RESIDUE_TO_CANONICAL     — CHARMM/AMBER residue variants → canonical.

All heavy lifting is Python dict indexing of
`(residue_number, atom_name) → atom_index` plus
`/topology/parent_atom_index` for methyl set expansion.

Outputs per protein (all under <protein_dir>/experimental_shifts/):
  - experimental_shifts_atom_view.csv    one row per CHARMM atom
  - experimental_shifts_bmrb_view.csv    one row per BMRB shift
  - experimental_shifts_refdb_view.csv   one row per RefDB shift
  - provenance.toml                       hashes, precheck status, counts
  - <pdb>_<bmrb>.bmrb.toml                per-protein translation (copy)
  - <pdb>_<bmrb>.refdb.toml               per-protein translation (copy)
  - <pdb>_<bmrb>.audit.json               audit output (copy)
  - <pdb>_<bmrb>.audit.md                 audit output (copy)
  - raw/bmr<id>_3.str                     BMRB NMR-STAR (copy)
  - raw/bmr<id>.str.corr                  RefDB .str.corr (copy)
  - raw/bmrb_<id>.csv                     parsed BMRB CSV (copy)
  - raw/refdb_<id>.csv                    parsed RefDB CSV (copy)
  - raw/refdb_<id>_header.toml            RefDB header metadata (copy)
  - raw/raw_provenance.toml               source hashes (copy)

Usage:
  ./.venv/bin/python pack_experimental_shifts.py --tree working
  ./.venv/bin/python pack_experimental_shifts.py --one 1DV0_4757 --tree working
"""

import argparse
import csv
import hashlib
import json
import shutil
import sys
import tomllib
from datetime import datetime, timezone
from pathlib import Path

import h5py

HERE = Path(__file__).parent
SOURCE_DATA = HERE / "source_data"

PROTEINS = [
    ("1DV0", 4757), ("1HS5", 4934), ("1HD6", 4820), ("1B1V", 4292),
    ("1HA9", 5028), ("1G26", 4656), ("1CBH", 192),   ("1I2V", 4976),
    ("1I8X", 4351), ("1BWX", 3449),
]

CALIBRATION_TREES = {
    "working": Path("/shared/2026Thesis/fleet_calibration-working"),
    "stats":   Path("/shared/2026Thesis/fleet_calibration-stats"),
    "backup":  Path("/shared/2026Thesis/fleet_calibration-backup"),
}

ELEMENT_LETTER = {1: "H", 6: "C", 7: "N", 8: "O", 15: "P", 16: "S"}

# ============================================================================
# Typed data constants. No string parsing in code — only dict lookups.
# ============================================================================

# H5 /atoms/atom_name carries CHARMM-style spellings at several positions
# because the library's NamingRegistry is missing rules for these
# CHARMM↔IUPAC divergences. See spec/ChangesRequiredBeforeProductionH5Run.md
# (2026-04-20 section) for the full story and the commented-out fix in
# src/NamingRegistry.cpp. Becomes unnecessary when the library fix
# activates AND the 10 proteins are re-extracted (neither currently planned).
# Keyed on (residue_type, IUPAC atom name). Value is the H5 atom name.
H5_NAME_FROM_IUPAC = {
    # ILE γ-carbon + δ-methyl — no library rules at all for ILE
    ("ILE", "CD1"):  "CD",
    ("ILE", "HD11"): "HD1",
    ("ILE", "HD12"): "HD2",
    ("ILE", "HD13"): "HD3",
    # ILE γ1-methylene
    ("ILE", "HG12"): "HG11",
    ("ILE", "HG13"): "HG12",
    # δ-methylene (ARG, LYS, PRO) — no library rule for δ
    ("ARG", "HD2"): "HD1",
    ("ARG", "HD3"): "HD2",
    ("LYS", "HD2"): "HD1",
    ("LYS", "HD3"): "HD2",
    ("PRO", "HD2"): "HD1",
    ("PRO", "HD3"): "HD2",
    # ε-methylene (LYS) — no library rule for ε
    ("LYS", "HE2"): "HE1",
    ("LYS", "HE3"): "HE2",
    # α-methylene (GLY) — no library rule for α
    ("GLY", "HA2"): "HA1",
    ("GLY", "HA3"): "HA2",
}

# BMRB rows that must bind via /topology/parent_atom_index because the H5
# lacks unique atom_names. Currently just ALA β-methyl: the library's
# wildcard β-methylene rule over-fires on ALA's 3-H methyl, producing
# an H5 with two atoms named HB3 and one named HB2 at every ALA residue.
# Parent-based binding preserves all three atom_indices (all H children
# of CB) and flags the resolution AMBIGUOUS_METHYL_ORDER — methyl Hs
# are chemically indistinguishable, so the individual HB1/HB2/HB3 → atom_index
# correspondence was never physically meaningful.
# Keyed on (residue_type, BMRB atom name); value is the parent heavy-atom
# IUPAC name, resolved to H5 via H5_NAME_FROM_IUPAC if applicable.
BMRB_METHYL_VIA_PARENT = {
    ("ALA", "HB1"): "CB",
    ("ALA", "HB2"): "CB",
    ("ALA", "HB3"): "CB",
}

# RefDB 2.1 methyl pseudo-atoms. IUPAC-1998 conventions where three Hs
# on a methyl are collapsed to a single pseudo name.
# Keyed by (canonical_residue_name, pseudo_atom_name). Value is
# (parent heavy-atom name in IUPAC, expected child H count).
# At resolution time, the parent atom is looked up in the H5 via
# (residue_number, parent_name-with-ILE-override-if-applicable), and
# all H atoms whose /topology/parent_atom_index equals that parent's
# index are returned as the expansion set.
REFDB_PSEUDO_PARENTS = {
    ("ALA", "HB"):  ("CB",  3),
    ("VAL", "HG1"): ("CG1", 3),
    ("VAL", "HG2"): ("CG2", 3),
    ("THR", "HG2"): ("CG2", 3),
    ("LEU", "HD1"): ("CD1", 3),
    ("LEU", "HD2"): ("CD2", 3),
    ("ILE", "HG2"): ("CG2", 3),
    ("ILE", "HD1"): ("CD1", 3),  # parent CD1 resolved to H5 "CD" via ILE override
    ("MET", "HE"):  ("CE",  3),
    ("LYS", "HZ"):  ("NZ",  3),
}

# Residue canonicalisation at the H5 ↔ BMRB/RefDB boundary. Belt-and-suspenders.
# If the H5 /residues/residue_name carries a CHARMM/AMBER variant, canonicalise
# before using in the constants above (BMRB/RefDB always use canonical HIS, CYS,
# etc.). Mirrors nmr::NamingRegistry::to_canonical_ (src/NamingRegistry.cpp).
H5_RESIDUE_TO_CANONICAL = {
    "HSP": "HIS", "HSD": "HIS", "HSE": "HIS",
    "HIE": "HIS", "HID": "HIS", "HIP": "HIS",
    "CYX": "CYS", "CYM": "CYS", "CYS2": "CYS",
    "ASH": "ASP", "ASPP": "ASP",
    "GLH": "GLU", "GLUP": "GLU",
    "LYN": "LYS", "ARN": "ARG", "TYM": "TYR", "MSE": "MET",
}


def canonical_residue(name):
    return H5_RESIDUE_TO_CANONICAL.get(name, name)


# ============================================================================
# Exceptions and small helpers
# ============================================================================

class PrecheckFailure(RuntimeError):
    pass


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_toml(path):
    with path.open("rb") as f:
        return tomllib.load(f)


def load_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


# ============================================================================
# H5 topology load
# ============================================================================

def load_h5_topology(h5_path):
    """
    Returns a dict of typed arrays + derived lookup dicts:
      atom_count, residue_count                (int scalars)
      atom_names, atom_elements,
      atom_residue_indices, parent_atom_index  (arrays, length atom_count)
      residue_numbers, residue_names           (arrays, length residue_count)
      name_lookup                              {(residue_number, atom_name): atom_index}
      children_by_parent                       {parent_atom_index: [child_atom_indices]}
      meta_attrs                               {...}
    """
    with h5py.File(h5_path, "r") as f:
        atom_names = [x.decode() if isinstance(x, bytes) else x
                      for x in f["/atoms/atom_name"][:]]
        atom_residue_indices = f["/atoms/residue_index"][:].tolist()
        atom_elements = f["/atoms/element"][:].tolist()
        parent_atom_index = f["/topology/parent_atom_index"][:].tolist()
        residue_numbers = f["/residues/residue_number"][:].tolist()
        residue_names = [x.decode() if isinstance(x, bytes) else x
                         for x in f["/residues/residue_name"][:]]
        meta_attrs = dict(f["/meta"].attrs)

    name_lookup = {}
    for atom_idx, aname in enumerate(atom_names):
        res_idx = atom_residue_indices[atom_idx]
        res_num = int(residue_numbers[res_idx])
        name_lookup[(res_num, aname)] = atom_idx

    children_by_parent = {}
    for atom_idx, p in enumerate(parent_atom_index):
        if p is None or p < 0:
            continue
        children_by_parent.setdefault(int(p), []).append(atom_idx)

    return {
        "atom_count": len(atom_names),
        "residue_count": len(residue_numbers),
        "atom_names": atom_names,
        "atom_elements": atom_elements,
        "atom_residue_indices": atom_residue_indices,
        "parent_atom_index": parent_atom_index,
        "residue_numbers": residue_numbers,
        "residue_names": residue_names,
        "name_lookup": name_lookup,
        "children_by_parent": children_by_parent,
        "meta_attrs": meta_attrs,
    }


# ============================================================================
# Resolution (H5 atom_index from BMRB/RefDB row)
# ============================================================================

def resolve_h5_atom_name(residue_name, atom_name_iupac):
    """
    Apply H5_NAME_FROM_IUPAC to convert an IUPAC atom name to the H5's
    atom_name. Identity for anything not in the table.
    """
    key = (residue_name, atom_name_iupac)
    return H5_NAME_FROM_IUPAC.get(key, atom_name_iupac)


def _expand_via_parent(h5, residue_number, parent_h5_name, expected_count,
                        label, diagnostic_ref):
    """
    Shared helper: given a parent heavy atom name, find the H children
    via /topology/parent_atom_index. Returns (child_indices, status, detail).
    """
    parent_key = (residue_number, parent_h5_name)
    if parent_key not in h5["name_lookup"]:
        return [], f"UNMATCHED_{label}_PARENT", (
            f"no parent atom {parent_h5_name} at residue {residue_number} "
            f"for {diagnostic_ref}"
        )
    parent_idx = h5["name_lookup"][parent_key]
    children = h5["children_by_parent"].get(parent_idx, [])
    h_children = [c for c in children if h5["atom_elements"][c] == 1]
    if expected_count is not None and len(h_children) != expected_count:
        return h_children, f"{label}_COUNT_MISMATCH", (
            f"parent {parent_h5_name} has {len(h_children)} H children, "
            f"expected {expected_count}"
        )
    return h_children, None, parent_h5_name


def resolve_bmrb_row(h5, residue_number, residue_name_h5, atom_name_bmrb):
    """
    Return (atom_indices_list, status, detail) for a single BMRB row.

    Status values:
      MATCHED                      — single-atom bind via name lookup
      AMBIGUOUS_METHYL_ORDER       — bound via parent (ALA β-methyl case)
      UNMATCHED_RESIDUE            — residue number outside topology
      UNMATCHED_NAME               — atom name not in H5 after translation
      UNMATCHED_METHYL_PARENT      — ALA CB parent missing (unlikely)
    """
    canonical_h5 = canonical_residue(residue_name_h5)

    # 1. ALA β-methyl library-bug handling: bind via parent.
    if (canonical_h5, atom_name_bmrb) in BMRB_METHYL_VIA_PARENT:
        parent_name_iupac = BMRB_METHYL_VIA_PARENT[(canonical_h5, atom_name_bmrb)]
        parent_h5_name = resolve_h5_atom_name(canonical_h5, parent_name_iupac)
        children, err_status, err_detail = _expand_via_parent(
            h5, residue_number, parent_h5_name, None,
            "METHYL", f"BMRB {atom_name_bmrb}"
        )
        if err_status:
            return children, err_status, err_detail
        return children, "AMBIGUOUS_METHYL_ORDER", err_detail

    # 2. Single-atom lookup with H5_NAME_FROM_IUPAC translation.
    h5_name = resolve_h5_atom_name(canonical_h5, atom_name_bmrb)
    if residue_number not in h5["residue_numbers"]:
        return [], "UNMATCHED_RESIDUE", f"residue {residue_number} not in topology"
    key = (residue_number, h5_name)
    if key not in h5["name_lookup"]:
        return [], "UNMATCHED_NAME", f"no atom named {h5_name} at residue {residue_number}"
    return [h5["name_lookup"][key]], "MATCHED", h5_name


def resolve_refdb_row(h5, residue_number, residue_name_h5, atom_name_refdb):
    """
    Return (atom_indices_list, status, detail) for a single RefDB row.

    Status values:
      MATCHED, PSEUDO_EXPANDED, UNMATCHED_RESIDUE, UNMATCHED_NAME,
      UNMATCHED_PSEUDO_PARENT, PSEUDO_COUNT_MISMATCH
    """
    canonical_h5 = canonical_residue(residue_name_h5)
    pseudo_key = (canonical_h5, atom_name_refdb)

    if pseudo_key in REFDB_PSEUDO_PARENTS:
        parent_name_iupac, expected_count = REFDB_PSEUDO_PARENTS[pseudo_key]
        parent_h5_name = resolve_h5_atom_name(canonical_h5, parent_name_iupac)
        children, err_status, err_detail = _expand_via_parent(
            h5, residue_number, parent_h5_name, expected_count,
            "PSEUDO", f"RefDB pseudo {atom_name_refdb}"
        )
        if err_status:
            return children, err_status, err_detail
        return children, "PSEUDO_EXPANDED", err_detail

    # Not a pseudo — single-atom lookup with H5_NAME_FROM_IUPAC translation.
    h5_name = resolve_h5_atom_name(canonical_h5, atom_name_refdb)
    if residue_number not in h5["residue_numbers"]:
        return [], "UNMATCHED_RESIDUE", f"residue {residue_number} not in topology"
    key = (residue_number, h5_name)
    if key not in h5["name_lookup"]:
        return [], "UNMATCHED_NAME", f"no atom named {h5_name} at residue {residue_number}"
    return [h5["name_lookup"][key]], "MATCHED", h5_name


# ============================================================================
# Prechecks (11 invariants, refuse-to-write on any failure)
# ============================================================================

def run_prechecks(ctx):
    """
    Returns list of (check_id, name, status, detail) tuples.
    Raises PrecheckFailure on any failure.
    """
    results = []
    tag = ctx["tag"]
    h5 = ctx["h5_topology"]
    raw_prov = ctx["raw_provenance"]
    audit = ctx["audit_json"]
    bmrb_rows = ctx["bmrb_rows"]
    refdb_rows = ctx["refdb_rows"]
    offset = ctx["bmrb_offset"]

    # --- Check 1: source-file hashes match raw_provenance
    for src in raw_prov["sources"]:
        expected = src["sha256"]
        actual_path = SOURCE_DATA / tag / "raw" / src["filename"]
        if not actual_path.is_file():
            raise PrecheckFailure(f"[{tag}] check1: missing source {src['filename']}")
        actual = sha256_file(actual_path)
        if actual != expected:
            raise PrecheckFailure(
                f"[{tag}] check1: {src['filename']} hash drift "
                f"(expected {expected[:12]}..., got {actual[:12]}...)"
            )
    results.append((1, "raw source hashes",
                    "PASS", f"{len(raw_prov['sources'])} sources"))

    # --- Check 2: reference.pdb exists and is hashable
    ref_pdb = ctx["protein_dir"] / "reference.pdb"
    if not ref_pdb.is_file():
        raise PrecheckFailure(f"[{tag}] check2: reference.pdb missing at {ref_pdb}")
    ctx["reference_pdb_sha256"] = sha256_file(ref_pdb)
    results.append((2, "reference.pdb present",
                    "PASS", f"sha256 {ctx['reference_pdb_sha256'][:12]}..."))

    # --- Check 3: H5 /meta attrs consistent with /atoms and /residues,
    # and H5 actually belongs to the protein we're packaging (guards
    # against accidentally-renamed or symlinked directories).
    m = h5["meta_attrs"]
    if int(m["n_atoms"]) != h5["atom_count"]:
        raise PrecheckFailure(
            f"[{tag}] check3: /meta.n_atoms={m['n_atoms']} != "
            f"len(/atoms)={h5['atom_count']}")
    if int(m["n_residues"]) != h5["residue_count"]:
        raise PrecheckFailure(
            f"[{tag}] check3: /meta.n_residues={m['n_residues']} != "
            f"len(/residues)={h5['residue_count']}")
    h5_protein_id = m.get("protein_id", "")
    if h5_protein_id != tag:
        raise PrecheckFailure(
            f"[{tag}] check3: /meta.protein_id={h5_protein_id!r} does not "
            f"match protein_dir tag {tag!r} — wrong H5 at this location?")
    results.append((3, "H5 /meta consistent and protein_id matches tag",
                    "PASS", f"{h5['atom_count']} atoms, {h5['residue_count']} residues, protein_id={h5_protein_id!r}"))

    # --- Check 4: reference.pdb residue sequence matches H5 /residues
    pdb_residues = {}
    for line in ref_pdb.read_text().split("\n"):
        if not line.startswith("ATOM"):
            continue
        res_name = line[17:20].strip()
        res_seq = int(line[22:26])
        pdb_residues[res_seq] = res_name
    if sorted(pdb_residues.keys()) != sorted(int(n) for n in h5["residue_numbers"]):
        raise PrecheckFailure(
            f"[{tag}] check4: reference.pdb residue_seq set differs from H5 /residues")
    for res_num, res_name in zip(h5["residue_numbers"], h5["residue_names"]):
        pdb_name = pdb_residues.get(int(res_num))
        # Compare after canonicalisation — H5 may carry HSP while reference.pdb has HSP too,
        # but we want to catch unexpected divergence.
        if canonical_residue(pdb_name) != canonical_residue(res_name):
            raise PrecheckFailure(
                f"[{tag}] check4: residue_name mismatch at seq={res_num}: "
                f"pdb={pdb_name} h5={res_name}")
    results.append((4, "reference.pdb matches H5 /residues",
                    "PASS", f"{h5['residue_count']} residues matched"))

    # --- Check 5: every BMRB row resolves (or is documented UNMATCHED in audit)
    resolved_statuses = ("MATCHED", "AMBIGUOUS_METHYL_ORDER", "PSEUDO_EXPANDED")
    unresolved_bmrb = []
    bmrb_resolved = []
    for row in bmrb_rows:
        try:
            residue_number = int(row["residue_seq"]) + offset
            idxs, status, detail = resolve_bmrb_row(
                h5, residue_number, row["residue_name"], row["atom_name"]
            )
            bmrb_resolved.append((row, residue_number, idxs, status, detail))
            if status not in resolved_statuses:
                unresolved_bmrb.append((row["residue_seq"], row["residue_name"],
                                          row["atom_name"], status))
        except (KeyError, ValueError) as exc:
            unresolved_bmrb.append((row.get("residue_seq"), row.get("residue_name"),
                                    row.get("atom_name"), f"parse_fail:{exc}"))
    # UNMATCHED rows are allowed only if audit classified same coord UNMATCHED_ROW
    # (document-and-carry). Use int keys on both sides — the audit stores
    # residue_seq as int; the BMRB CSV stores it as str.
    audit_unmatched_by_coord = {
        (int(e["row"]["residue_seq"]), e["row"]["atom_name"]): e
        for e in audit["bmrb_audit"] if e["status"] == "UNMATCHED_ROW"
    }
    unexpected = []
    for seq, resname, atom, status in unresolved_bmrb:
        if (int(seq), atom) not in audit_unmatched_by_coord:
            unexpected.append((seq, resname, atom, status))
    if unexpected:
        raise PrecheckFailure(
            f"[{tag}] check5: {len(unexpected)} BMRB rows failed to resolve and "
            f"are NOT in the audit's UNMATCHED list; samples: {unexpected[:5]}")
    n_matched = sum(1 for r in bmrb_resolved if r[3] == "MATCHED")
    n_ambig = sum(1 for r in bmrb_resolved if r[3] == "AMBIGUOUS_METHYL_ORDER")
    results.append((5, "every BMRB row resolves or is documented UNMATCHED",
                    "PASS", f"{n_matched} matched + {n_ambig} methyl-ambiguous "
                    f"+ {len(unresolved_bmrb)} documented UNMATCHED"))
    ctx["bmrb_resolved"] = bmrb_resolved

    # --- Check 6: element consistency for every resolved BMRB row
    elem_mismatches = []
    for row, residue_number, idxs, status, _ in bmrb_resolved:
        if status not in ("MATCHED", "AMBIGUOUS_METHYL_ORDER"):
            continue
        row_type = row.get("atom_type", "")
        for idx in idxs:
            h5_elem = h5["atom_elements"][idx]
            h5_letter = ELEMENT_LETTER.get(int(h5_elem), str(h5_elem))
            if row_type and row_type != h5_letter:
                elem_mismatches.append((row["residue_seq"], row["atom_name"],
                                        row_type, h5_letter))
    if elem_mismatches:
        raise PrecheckFailure(
            f"[{tag}] check6: {len(elem_mismatches)} BMRB-H5 element mismatches; "
            f"samples: {elem_mismatches[:5]}")
    n_elem_checked = sum(len(r[2]) for r in bmrb_resolved
                         if r[3] in ("MATCHED", "AMBIGUOUS_METHYL_ORDER"))
    results.append((6, "BMRB atom_type matches H5 element",
                    "PASS", f"{n_elem_checked} atom bindings type-checked"))

    # --- Check 7: HSP/HSD/HSE residue types agree between audit and H5
    # The library canonicalises /residues/residue_name — all HIS variants
    # are recorded as "HIS". Protonation is encoded by atom inventory:
    #   HSP: HD1 present AND HE2 present (doubly protonated)
    #   HSD: HD1 present, HE2 absent     (Nδ1 protonated only)
    #   HSE: HD1 absent,  HE2 present    (Nε2 protonated only)
    # We check:
    #   (a) H5 residue_name canonicalises to "HIS"
    #   (b) inferred protonation (from HD1/HE2 presence) matches the
    #       CHARMM variant the audit asserts
    h5_res_by_num = dict(zip([int(n) for n in h5["residue_numbers"]], h5["residue_names"]))
    has_hd1 = {}
    has_he2 = {}
    for atom_idx, aname in enumerate(h5["atom_names"]):
        rn = int(h5["residue_numbers"][h5["atom_residue_indices"][atom_idx]])
        if aname == "HD1":
            has_hd1[rn] = True
        if aname == "HE2":
            has_he2[rn] = True

    hsp_checks = []
    for entry in audit["bmrb_audit"]:
        charmm_variant = None
        for a in entry["anomalies"]:
            if "RESIDUE_ALIAS_APPLIED: CHARMM 'HSP'" in a:
                charmm_variant = "HSP"; break
            if "RESIDUE_ALIAS_APPLIED: CHARMM 'HSD'" in a:
                charmm_variant = "HSD"; break
            if "RESIDUE_ALIAS_APPLIED: CHARMM 'HSE'" in a:
                charmm_variant = "HSE"; break
        if not charmm_variant:
            continue
        seq = int(entry["row"].get("charmm_seq", entry["row"]["residue_seq"]))
        h5_name = h5_res_by_num.get(seq, "?")
        if canonical_residue(h5_name) != "HIS":
            raise PrecheckFailure(
                f"[{tag}] check7: audit asserts HIS at seq={seq}, "
                f"H5 records {h5_name!r}")
        hd1 = has_hd1.get(seq, False)
        he2 = has_he2.get(seq, False)
        if hd1 and he2:
            inferred = "HSP"
        elif hd1 and not he2:
            inferred = "HSD"
        elif he2 and not hd1:
            inferred = "HSE"
        else:
            inferred = "NO_RING_NH"
        if inferred != charmm_variant:
            raise PrecheckFailure(
                f"[{tag}] check7: audit claims {charmm_variant} at seq={seq}, "
                f"H5 atom inventory (HD1={hd1}, HE2={he2}) implies {inferred}")
        hsp_checks.append((seq, charmm_variant))
    results.append((7, "audit HSP/HSD/HSE variant matches H5 atom inventory",
                    "PASS", f"{len(hsp_checks)} HIS residues checked "
                    f"(variants: {sorted(set(v for _, v in hsp_checks))})"))

    # --- Check 8: RefDB pseudo expansion completeness
    pseudo_mismatches = []
    refdb_resolved = []
    for row in refdb_rows:
        try:
            residue_number = int(row["residue_seq"]) + ctx["refdb_offset"]
            idxs, status, detail = resolve_refdb_row(
                h5, residue_number, row["residue_name"], row["atom_name"]
            )
            refdb_resolved.append((row, residue_number, idxs, status, detail))
            if status == "PSEUDO_COUNT_MISMATCH":
                pseudo_mismatches.append((row["residue_seq"], row["residue_name"],
                                          row["atom_name"], detail))
        except (KeyError, ValueError):
            pass
    if pseudo_mismatches:
        raise PrecheckFailure(
            f"[{tag}] check8: {len(pseudo_mismatches)} RefDB pseudo count mismatches; "
            f"samples: {pseudo_mismatches[:5]}")
    results.append((8, "RefDB pseudo expansion matches H5 child count",
                    "PASS", f"{sum(1 for r in refdb_resolved if r[3]=='PSEUDO_EXPANDED')} pseudos OK"))
    ctx["refdb_resolved"] = refdb_resolved

    # --- Check 9: no duplicate BMRB row bindings to the same atom_index
    # AMBIGUOUS_METHYL_ORDER is expected to have N rows binding to N atoms
    # in the same set (ALA β-methyl: 3 rows × 3 atoms). That's not a
    # collision — it's the design. Collision detection is only for MATCHED
    # single-atom bindings.
    bmrb_seen = {}
    collisions = []
    for row, _, idxs, status, _ in bmrb_resolved:
        if status != "MATCHED":
            continue
        idx = idxs[0]
        coord = (row["residue_seq"], row["atom_name"])
        if idx in bmrb_seen and bmrb_seen[idx] != coord:
            collisions.append((idx, bmrb_seen[idx], coord))
        bmrb_seen[idx] = coord
    if collisions:
        raise PrecheckFailure(
            f"[{tag}] check9: BMRB binding collisions: {collisions[:5]}")
    results.append((9, "no BMRB atom_index collisions (single-atom bindings)",
                    "PASS", f"{len(bmrb_seen)} unique single-atom bindings"))

    # --- Check 10: cross-db agreement
    # For single-atom bindings (MATCHED on both sides): require same atom_index.
    # For AMBIGUOUS_METHYL_ORDER vs PSEUDO_EXPANDED on the same methyl group:
    # both should resolve to the same set of atoms. Compare as sets.
    bmrb_by_triple = {}        # (rn, canonical, atom) -> atom_index (MATCHED only)
    bmrb_methyl_sets = {}      # (rn, canonical, parent_CB_name) -> frozenset(indices)
    for row, rn, idxs, status, detail in bmrb_resolved:
        canonical = canonical_residue(row["residue_name"])
        if status == "MATCHED":
            bmrb_by_triple[(rn, canonical, row["atom_name"])] = idxs[0]
        elif status == "AMBIGUOUS_METHYL_ORDER":
            bmrb_methyl_sets[(rn, canonical)] = frozenset(idxs)
    xdb_mismatches = []
    for row, rn, idxs, status, detail in refdb_resolved:
        canonical = canonical_residue(row["residue_name"])
        if status == "MATCHED":
            key = (rn, canonical, row["atom_name"])
            if key in bmrb_by_triple and idxs[0] != bmrb_by_triple[key]:
                xdb_mismatches.append((key, bmrb_by_triple[key], idxs[0]))
        elif status == "PSEUDO_EXPANDED":
            # Compare set against BMRB's methyl-set for the same residue.
            key = (rn, canonical)
            if key in bmrb_methyl_sets and frozenset(idxs) != bmrb_methyl_sets[key]:
                xdb_mismatches.append((key, bmrb_methyl_sets[key], frozenset(idxs)))
    if xdb_mismatches:
        raise PrecheckFailure(
            f"[{tag}] check10: cross-db atom_index mismatches: {xdb_mismatches[:5]}")
    results.append((10, "BMRB and RefDB agree on atom binding at shared keys",
                    "PASS", f"{len(bmrb_by_triple)} single + {len(bmrb_methyl_sets)} methyl-set bindings"))

    # --- Check 11: for every residue with a CHARMM-vs-IUPAC gap in the
    # library, the H5 must carry a uniform style — never a mix.
    # Detects partial library activation that would silently corrupt
    # downstream binding at a specific residue.
    #
    # CHARMM and IUPAC name sets often overlap on gap-residues (e.g.,
    # for ARG δ-methylene, CHARMM uses {HD1, HD2} and IUPAC uses
    # {HD2, HD3} — HD2 is in both). We can only decide the H5's style
    # from the DISAMBIGUATING names (charmm-only and iupac-only).
    gap_residues = sorted({k[0] for k in H5_NAME_FROM_IUPAC})
    per_residue_style = {}
    mixed = []
    for residue in gap_residues:
        iupac_names  = {k[1] for k in H5_NAME_FROM_IUPAC if k[0] == residue}
        charmm_names = {v for (r, _), v in H5_NAME_FROM_IUPAC.items() if r == residue}
        charmm_only = charmm_names - iupac_names
        iupac_only  = iupac_names  - charmm_names
        cc = 0   # atoms with a charmm-only marker name
        ic = 0   # atoms with an iupac-only marker name
        for atom_idx, aname in enumerate(h5["atom_names"]):
            res_idx = h5["atom_residue_indices"][atom_idx]
            if canonical_residue(h5["residue_names"][res_idx]) != residue:
                continue
            if aname in charmm_only:
                cc += 1
            elif aname in iupac_only:
                ic += 1
        if cc and ic:
            mixed.append((residue, cc, ic))
        per_residue_style[residue] = (
            "charmm" if cc and not ic else
            "iupac"  if ic and not cc else
            "absent" if not cc and not ic else
            "MIXED"
        )
    if mixed:
        raise PrecheckFailure(
            f"[{tag}] check11: H5 naming is mixed CHARMM/IUPAC at gap "
            f"residues {mixed} — partial library activation suspected")
    summary = ", ".join(f"{r}={per_residue_style[r]}"
                         for r in gap_residues
                         if per_residue_style[r] != "absent")
    results.append((11, "H5 naming style uniform at every gap-residue",
                    "PASS", summary or "no gap residues present in this protein"))

    return results


# ============================================================================
# View builders
# ============================================================================

def build_views(ctx):
    h5 = ctx["h5_topology"]
    bmrb_resolved = ctx["bmrb_resolved"]
    refdb_resolved = ctx["refdb_resolved"]

    # Index per atom_index: list of (row, status, detail) for each db.
    # Both BMRB and RefDB may bind a row to multiple atom_indices:
    # BMRB via AMBIGUOUS_METHYL_ORDER, RefDB via PSEUDO_EXPANDED.
    bmrb_by_atom = {}
    for row, rn, idxs, status, detail in bmrb_resolved:
        if not idxs:
            continue
        for idx in idxs:
            bmrb_by_atom.setdefault(idx, []).append((row, status, detail))

    refdb_by_atom = {}
    for row, rn, idxs, status, detail in refdb_resolved:
        for idx in idxs:
            refdb_by_atom.setdefault(idx, []).append((row, status, detail))

    # ---- atom-view ----
    atom_view = []
    for ai in range(h5["atom_count"]):
        res_idx = h5["atom_residue_indices"][ai]
        res_num = int(h5["residue_numbers"][res_idx])
        res_name_h5 = h5["residue_names"][res_idx]
        atom_name_h5 = h5["atom_names"][ai]
        element_num = int(h5["atom_elements"][ai])
        element = ELEMENT_LETTER.get(element_num, str(element_num))

        bmrb_hits = bmrb_by_atom.get(ai, [])
        refdb_hits = refdb_by_atom.get(ai, [])

        bmrb_row = bmrb_hits[0][0] if bmrb_hits else None
        refdb_row = refdb_hits[0][0] if refdb_hits else None
        refdb_status = refdb_hits[0][1] if refdb_hits else ""

        if not bmrb_hits and not refdb_hits:
            binding = "UNMATCHED_TOPOLOGY"  # atom has no experimental shift
        elif bmrb_hits and refdb_hits:
            binding = f"BMRB+REFDB{'_PSEUDO' if refdb_status == 'PSEUDO_EXPANDED' else ''}"
        elif bmrb_hits:
            binding = "BMRB_ONLY"
        else:
            binding = f"REFDB_ONLY{'_PSEUDO' if refdb_status == 'PSEUDO_EXPANDED' else ''}"

        atom_view.append({
            "atom_index": ai,
            "charmm_residue_seq": res_num,
            "charmm_residue_name": res_name_h5,
            "charmm_residue_canonical": canonical_residue(res_name_h5),
            "h5_atom_name": atom_name_h5,
            "element": element,
            "bmrb_residue_seq": bmrb_row["residue_seq"] if bmrb_row else "",
            "bmrb_residue_name": bmrb_row["residue_name"] if bmrb_row else "",
            "bmrb_atom_name": bmrb_row["atom_name"] if bmrb_row else "",
            "bmrb_shift": bmrb_row.get("shift_value", "") if bmrb_row else "",
            "bmrb_ambig_code": bmrb_row.get("ambiguity_code", "") if bmrb_row else "",
            "refdb_residue_seq": refdb_row["residue_seq"] if refdb_row else "",
            "refdb_atom_name": refdb_row["atom_name"] if refdb_row else "",
            "refdb_shift": refdb_row.get("shift_value", "") if refdb_row else "",
            "refdb_ambig_code": refdb_row.get("ambiguity_code", "") if refdb_row else "",
            "binding": binding,
            "bmrb_row_count": len(bmrb_hits),
            "refdb_row_count": len(refdb_hits),
        })

    # ---- bmrb-view ----
    bmrb_view = []
    for row, rn, idxs, status, detail in bmrb_resolved:
        bound = ";".join(str(i) for i in idxs) if idxs else ""
        h5_names = ";".join(h5["atom_names"][i] for i in idxs) if idxs else ""
        bmrb_view.append({
            "bmrb_source_row_id": row.get("source_row_id", ""),
            "bmrb_residue_seq": row["residue_seq"],
            "bmrb_residue_name": row["residue_name"],
            "bmrb_atom_name": row["atom_name"],
            "bmrb_atom_type": row.get("atom_type", ""),
            "bmrb_shift": row.get("shift_value", ""),
            "bmrb_ambig_code": row.get("ambiguity_code", ""),
            "charmm_residue_seq": rn,
            "bound_atom_indices": bound,
            "h5_atom_names": h5_names,
            "status": status,
            "status_detail": "" if status in ("MATCHED", "AMBIGUOUS_METHYL_ORDER") else detail,
        })

    # ---- refdb-view ----
    refdb_view = []
    for row, rn, idxs, status, detail in refdb_resolved:
        bound = ";".join(str(i) for i in idxs) if idxs else ""
        h5_names = ";".join(h5["atom_names"][i] for i in idxs) if idxs else ""
        refdb_view.append({
            "refdb_source_row_id": row.get("source_row_id", ""),
            "refdb_residue_seq": row["residue_seq"],
            "refdb_residue_name": row["residue_name"],
            "refdb_atom_name": row["atom_name"],
            "refdb_atom_type": row.get("atom_type", ""),
            "refdb_shift": row.get("shift_value", ""),
            "refdb_ambig_code": row.get("ambiguity_code", ""),
            "charmm_residue_seq": rn,
            "bound_atom_indices": bound,
            "h5_atom_names": h5_names,
            "status": status,
            "status_detail": "" if status in ("MATCHED", "PSEUDO_EXPANDED") else detail,
        })

    return atom_view, bmrb_view, refdb_view


# ============================================================================
# Writers
# ============================================================================

def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


ATOM_VIEW_FIELDS = [
    "atom_index",
    "charmm_residue_seq", "charmm_residue_name", "charmm_residue_canonical",
    "h5_atom_name", "element",
    "bmrb_residue_seq", "bmrb_residue_name", "bmrb_atom_name",
    "bmrb_shift", "bmrb_ambig_code",
    "refdb_residue_seq", "refdb_atom_name",
    "refdb_shift", "refdb_ambig_code",
    "binding", "bmrb_row_count", "refdb_row_count",
]
BMRB_VIEW_FIELDS = [
    "bmrb_source_row_id",
    "bmrb_residue_seq", "bmrb_residue_name", "bmrb_atom_name",
    "bmrb_atom_type", "bmrb_shift", "bmrb_ambig_code",
    "charmm_residue_seq", "bound_atom_indices", "h5_atom_names",
    "status", "status_detail",
]
REFDB_VIEW_FIELDS = [
    "refdb_source_row_id",
    "refdb_residue_seq", "refdb_residue_name", "refdb_atom_name",
    "refdb_atom_type", "refdb_shift", "refdb_ambig_code",
    "charmm_residue_seq", "bound_atom_indices", "h5_atom_names",
    "status", "status_detail",
]


def write_provenance(out_dir, ctx, precheck_results):
    path = out_dir / "provenance.toml"
    tag = ctx["tag"]
    h5 = ctx["h5_topology"]
    now_iso = datetime.now(timezone.utc).isoformat()
    packager_path = Path(__file__)
    packager_sha = sha256_file(packager_path) if packager_path.is_file() else "unknown"

    lines = [
        "# experimental_shifts bundle — provenance",
        "# Generated by pack_experimental_shifts.py (2026-04-20 pipeline).",
        "# Scope: 10-protein calibration; fleet replication deferred.",
        "",
        f'protein_tag = "{tag}"',
        f'bmrb_id = {ctx["bmrb_id"]}',
        f'generated_at = "{now_iso}"',
        f'packager_script = "{packager_path.name}"',
        f'packager_script_sha256 = "{packager_sha}"',
        f'reference_pdb_sha256 = "{ctx["reference_pdb_sha256"]}"',
        f'h5_file = "{ctx["h5_path"].relative_to(ctx["protein_dir"])}"',
        f'h5_n_atoms = {h5["atom_count"]}',
        f'h5_n_residues = {h5["residue_count"]}',
        f'bmrb_row_count = {len(ctx["bmrb_rows"])}',
        f'refdb_row_count = {len(ctx["refdb_rows"])}',
        f'bmrb_offset_applied = {ctx["bmrb_offset"]}',
        f'refdb_offset_applied = {ctx["refdb_offset"]}',
        "",
        "[precheck_results]",
    ]
    for cid, name, status, detail in precheck_results:
        lines.append(f'# check {cid:2d} {status}: {name} — {detail}')
    lines.append("")
    lines.append("[raw_sources]")
    for src in ctx["raw_provenance"]["sources"]:
        lines.append("[[raw_sources.entries]]")
        for k, v in src.items():
            if isinstance(v, int):
                lines.append(f'{k} = {v}')
            else:
                lines.append(f'{k} = "{v}"')
        lines.append("")

    # The three data constants (recorded verbatim for auditability — a
    # future reader can verify the code matched this provenance).
    lines.append("[constants]")
    lines.append("# H5_NAME_FROM_IUPAC: H5 carries CHARMM names at six registry-gap")
    lines.append("# positions (ILE quirks, δ/ε/α-methylenes). Temporary workaround;")
    lines.append("# see spec/ChangesRequiredBeforeProductionH5Run.md (2026-04-20).")
    for k, v in sorted(H5_NAME_FROM_IUPAC.items()):
        lines.append(f'# H5_NAME_FROM_IUPAC[{k}] = "{v}"')
    lines.append("# BMRB_METHYL_VIA_PARENT: ALA β-methyl bound via parent_atom_index")
    lines.append("# because library wildcard-β bug made atom_names non-unique in H5.")
    for k, v in sorted(BMRB_METHYL_VIA_PARENT.items()):
        lines.append(f'# BMRB_METHYL_VIA_PARENT[{k}] = "{v}"')
    lines.append("# REFDB_PSEUDO_PARENTS: RefDB 2.1 methyl pseudos → parent heavy atom.")
    for k, v in sorted(REFDB_PSEUDO_PARENTS.items()):
        lines.append(f'# REFDB_PSEUDO_PARENTS[{k}] = {v}')
    lines.append("")

    path.write_text("\n".join(lines))


# ============================================================================
# Per-protein pipeline
# ============================================================================

def pack_one(protein_dir, pdb_id, bmrb_id):
    tag = f"{pdb_id}_{bmrb_id}"
    out_dir = protein_dir / "experimental_shifts"
    out_dir.mkdir(exist_ok=True)
    (out_dir / "raw").mkdir(exist_ok=True)

    src_dir = SOURCE_DATA / tag
    raw_prov_path = src_dir / "raw_provenance.toml"
    bmrb_toml_path = HERE / f"{tag}.bmrb.toml"
    refdb_toml_path = HERE / f"{tag}.refdb.toml"
    audit_json_path = HERE / f"{tag}.audit.json"
    audit_md_path = HERE / f"{tag}.audit.md"
    bmrb_csv_path = src_dir / f"bmrb_{bmrb_id}.csv"
    refdb_csv_path = src_dir / f"refdb_{bmrb_id}.csv"
    refdb_header_toml = src_dir / f"refdb_{bmrb_id}_header.toml"
    h5_path = protein_dir / "analysis_output" / f"{tag}_analysis.h5"

    for p in (raw_prov_path, bmrb_toml_path, refdb_toml_path,
              audit_json_path, audit_md_path,
              bmrb_csv_path, refdb_csv_path, h5_path):
        if not p.is_file():
            raise PrecheckFailure(f"[{tag}] missing input: {p}")

    bmrb_toml = load_toml(bmrb_toml_path)
    refdb_toml = load_toml(refdb_toml_path)
    bmrb_offset = bmrb_toml.get("residue_seq_offset", {}).get("offset", 0)
    refdb_offset = refdb_toml.get("residue_seq_offset", {}).get("offset", 0)

    # Guard: BMRB and RefDB translation tables must declare the same
    # residue_seq_offset for a single protein — both apply to the same
    # raw residue numbering. Disagreement is a reviewer error that would
    # silently misalign cross-db binding.
    if bmrb_offset != refdb_offset:
        raise PrecheckFailure(
            f"[{tag}] BMRB TOML residue_seq_offset={bmrb_offset} but "
            f"RefDB TOML residue_seq_offset={refdb_offset} — must agree "
            f"for a single protein (both apply to raw depositor residue_seq)."
        )

    ctx = {
        "tag": tag,
        "pdb_id": pdb_id,
        "bmrb_id": bmrb_id,
        "protein_dir": protein_dir,
        "h5_path": h5_path,
        "raw_provenance": load_toml(raw_prov_path),
        "audit_json": json.loads(audit_json_path.read_text()),
        "bmrb_rows": load_csv(bmrb_csv_path),
        "refdb_rows": load_csv(refdb_csv_path),
        "bmrb_offset": bmrb_offset,
        "refdb_offset": refdb_offset,
        "h5_topology": load_h5_topology(h5_path),
    }

    precheck_results = run_prechecks(ctx)
    atom_view, bmrb_view, refdb_view = build_views(ctx)

    write_csv(out_dir / "experimental_shifts_atom_view.csv",
              ATOM_VIEW_FIELDS, atom_view)
    write_csv(out_dir / "experimental_shifts_bmrb_view.csv",
              BMRB_VIEW_FIELDS, bmrb_view)
    write_csv(out_dir / "experimental_shifts_refdb_view.csv",
              REFDB_VIEW_FIELDS, refdb_view)
    write_provenance(out_dir, ctx, precheck_results)

    # Copies
    for src_path in (bmrb_toml_path, refdb_toml_path,
                     audit_json_path, audit_md_path):
        shutil.copy2(src_path, out_dir / src_path.name)
    for src_path in (src_dir / "raw" / f"bmr{bmrb_id}_3.str",
                     src_dir / "raw" / f"bmr{bmrb_id}.str.corr",
                     raw_prov_path, bmrb_csv_path, refdb_csv_path):
        shutil.copy2(src_path, out_dir / "raw" / src_path.name)
    if refdb_header_toml.is_file():
        shutil.copy2(refdb_header_toml, out_dir / "raw" / refdb_header_toml.name)

    return {
        "tag": tag,
        "prechecks": precheck_results,
        "atom_view_rows": len(atom_view),
        "bmrb_view_rows": len(bmrb_view),
        "refdb_view_rows": len(refdb_view),
        "out_dir": out_dir,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", choices=list(CALIBRATION_TREES), default="working")
    ap.add_argument("--one", metavar="TAG", help="pack only the named protein")
    args = ap.parse_args()

    base = CALIBRATION_TREES[args.tree]
    protein_list = PROTEINS
    if args.one:
        protein_list = [(p, b) for (p, b) in PROTEINS if f"{p}_{b}" == args.one]
        if not protein_list:
            print(f"unknown --one target: {args.one}", file=sys.stderr)
            sys.exit(2)

    successes, failures = [], []
    for pdb_id, bmrb_id in protein_list:
        protein_dir = base / f"{pdb_id}_{bmrb_id}"
        if not protein_dir.is_dir():
            failures.append((f"{pdb_id}_{bmrb_id}", f"protein dir missing: {protein_dir}"))
            continue
        try:
            r = pack_one(protein_dir, pdb_id, bmrb_id)
            successes.append(r)
            print(f"[{r['tag']}] PACKED → {r['out_dir']}  "
                  f"(atoms={r['atom_view_rows']}, bmrb={r['bmrb_view_rows']}, "
                  f"refdb={r['refdb_view_rows']})")
        except PrecheckFailure as e:
            failures.append((f"{pdb_id}_{bmrb_id}", str(e)))
            print(f"PRECHECK FAILED: {e}", file=sys.stderr)

    print(f"\n{len(successes)} packed, {len(failures)} failed")
    if failures:
        for tag, msg in failures:
            print(f"  FAIL {tag}: {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
