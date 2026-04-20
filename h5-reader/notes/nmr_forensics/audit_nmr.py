#!/usr/bin/env python3
"""
audit_nmr.py — per-protein, per-database forensic audit.

Methodology (Huxley principle, 2026-04-19):
  - No filtering of BMRB/RefDB input based on traditional knowledge.
  - Every row (BMRB or RefDB) is classified; every anomaly complains loudly.
  - CHARMM topology is trusted (one rigorous pipeline); databases are not.
  - Pseudo-atom expansion is NOT performed at audit time; ONE_TO_MANY flagged.
  - Ambiguity codes ≥ 2 surface as STEREO_AMBIGUOUS.
  - Residue-name aliases (HSP/HSD/HSE → HIS) fire RESIDUE_ALIAS_APPLIED
    on every affected row — loud, not silent.

Per-protein-per-database translation tables (2026-04-19 refactor):
  - BMRB and RefDB entries are not assumed to share naming conventions.
  - Each protein gets two tables: <pdb>_<bmrb>.bmrb.toml and
    <pdb>_<bmrb>.refdb.toml, living alongside this script.
  - A starter at _starter_translation.toml records the default modern-IUPAC
    conventions; per-protein tables start as copies and diverge per-entry.
  - No silent fallback: if a per-protein table is missing, audit refuses.

Outputs:
  <pdb>_<bmrb>.audit.json — complete machine-readable record
  <pdb>_<bmrb>.audit.md   — human-reviewable narrative + per-row tables

Usage:
  python3 audit_nmr.py <protein_dir> <bmrb_id>
"""

import csv
import json
import sys
import tomllib
from collections import defaultdict
from pathlib import Path

# Per-protein input CSVs live under `source_data/<pdb>_<bmrb>/` and are
# produced by `fetch_bmrb.py` + `extract_refdb.py` + `parse_shifts.py`.
# The old consolidated `chemical_shifts_685.csv` / `refdb_corrected_shifts_685.csv`
# were RefDB-derived (not BMRB raw) and dropped ~10-16% of the actual
# BMRB deposition rows. Replaced 2026-04-20. See SUMMARY.md.
HERE = Path(__file__).parent
SOURCE_DATA = HERE / "source_data"


# ----------------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------------

def parse_pdb(pdb_path):
    residues = {}
    with pdb_path.open() as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_seq = int(line[22:26])
            element_field = line[76:78].strip() if len(line) >= 78 else ""
            element = element_field if element_field else atom_name[0]
            residues.setdefault(res_seq, {"name": res_name, "atoms": {}})
            residues[res_seq]["atoms"][atom_name] = element
    return residues


def load_csv_rows(csv_path):
    """Read all rows from a per-protein CSV. No filtering — the file is already
    per-protein by construction (see parse_shifts.py)."""
    rows = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def build_value_lookup(rows):
    """Build {(residue_seq_int, atom_name): shift_value_float} for cross-db
    comparison. Skip rows with missing/unparseable values."""
    lookup = {}
    for r in rows:
        try:
            seq = int(r["residue_seq"])
            atom = r["atom_name"]
            val_s = (r.get("shift_value") or "").strip()
            if val_s and val_s != ".":
                lookup[(seq, atom)] = float(val_s)
        except (KeyError, ValueError):
            continue
    return lookup


def load_translation(path):
    with path.open("rb") as f:
        return tomllib.load(f)


REQUIRED_PROVENANCE_KEYS = (
    "protein", "database", "reviewed_date", "reviewer",
    "metadata_sources_consulted", "findings_summary",
)


def validate_provenance(translation, path):
    """
    Refuse to run unless the per-protein TOML carries a [provenance]
    section with all required keys filled in. This is the structural
    check that catches a table still in starter-derived state — no
    silent acceptance of unreviewed tables.

    Raises ValueError with a specific, actionable message on failure.
    """
    prov = translation.get("provenance")
    if prov is None:
        raise ValueError(
            f"{path.name}: missing [provenance] section. Per the 2026-04-19"
            f" discipline (PRINCIPLES.md §7), every per-protein translation"
            f" table must carry a [provenance] section with keys: "
            + ", ".join(REQUIRED_PROVENANCE_KEYS)
            + f". See 1DV0_4757.bmrb.toml for the worked example."
        )
    missing = [k for k in REQUIRED_PROVENANCE_KEYS if not prov.get(k)]
    if missing:
        raise ValueError(
            f"{path.name}: [provenance] section missing or empty for keys:"
            f" {missing}. All six required keys must be populated before"
            f" audit can run. Consult BMRB deposition metadata at"
            f" bmrb.io/data_library/summary/?bmrbId=<id> for provenance"
            f" inputs. See 1DV0_4757.bmrb.toml for the worked example."
        )


# ----------------------------------------------------------------------------
# Translation / reverse map
# ----------------------------------------------------------------------------

def apply_translation(translation, res_name, charmm_atom):
    """CHARMM atom in a residue → expected database atom. Identity fallback."""
    if charmm_atom in translation.get("backbone", {}):
        return translation["backbone"][charmm_atom]
    res_map = translation.get("residue", {}).get(res_name, {})
    if charmm_atom in res_map:
        return res_map[charmm_atom]
    return charmm_atom


def build_reverse_map(translation, topology):
    """Per residue_seq: {db_atom_name: [charmm_atom_names]}."""
    reverse = {}
    for seq, info in topology.items():
        res_name = info["name"]
        rev = defaultdict(list)
        for charmm_atom in info["atoms"]:
            db_atom = apply_translation(translation, res_name, charmm_atom)
            if db_atom:
                rev[db_atom].append(charmm_atom)
        reverse[seq] = dict(rev)
    return reverse


# ----------------------------------------------------------------------------
# Single-database audit pass
# ----------------------------------------------------------------------------

def audit_rows(db_rows, topology, translation, residue_aliases, db_label,
               value_col, include_refdb_delta=False, bmrb_lookup=None):
    """
    Run one audit pass against db_rows using translation table.

    Returns (row_audit_list, matched_charmm_set, offset_applied).
    row_audit_list entries: {row, status, matched_charmm, anomalies}.

    If the translation table declares `[residue_seq_offset] offset = N`,
    every database residue_seq is shifted by N before topology lookup.
    Each affected row's `row.charmm_seq` records the post-offset value;
    `row.residue_seq` preserves the raw database number.
    """
    offset = translation.get("residue_seq_offset", {}).get("offset", 0)
    reverse = build_reverse_map(translation, topology)
    row_audit = []
    matched_charmm = set()

    for r in db_rows:
        db_seq_raw = int(r["residue_seq"])
        bmrb_seq = db_seq_raw + offset
        db_resname = r["residue_name"]
        db_atomname = r["atom_name"]
        db_nucleus = r["atom_type"]
        value_str = (r.get(value_col) or "").strip()
        try:
            value = float(value_str) if value_str else None
        except ValueError:
            value = None

        ambig = None
        if "ambiguity_code" in r:
            s = (r["ambiguity_code"] or "").strip()
            try:
                ambig = int(s) if s else None
            except ValueError:
                ambig = None

        row_dict = {
            "residue_seq": db_seq_raw,
            "residue_name": db_resname,
            "atom_name": db_atomname,
            "atom_type": db_nucleus,
            "value": value,
            "value_col": value_col,
            "ambiguity_code": ambig,
        }
        if offset != 0:
            row_dict["charmm_seq"] = bmrb_seq

        entry = {
            "row": row_dict,
            "matched_charmm": None,
            "status": None,
            "anomalies": [],
        }

        if include_refdb_delta:
            # RefDB .str.corr files store only the corrected value. Original
            # (BMRB raw) comes from a cross-database lookup keyed on
            # (residue_seq, atom_name). Pre- vs post-offset: RefDB and BMRB
            # both use the raw deposited residue_seq, so the key matches
            # directly without offset.
            key = (db_seq_raw, db_atomname)
            if bmrb_lookup is not None and key in bmrb_lookup and value is not None:
                bmrb_value = bmrb_lookup[key]
                entry["row"]["bmrb_value"] = bmrb_value
                entry["row"]["refdb_value"] = value
                entry["row"]["delta"] = round(value - bmrb_value, 6)
                if abs(value - bmrb_value) > 1e-4:
                    entry["anomalies"].append(
                        f"REFDB_CORRECTED_VALUE_DIFFERS: bmrb_raw={bmrb_value}"
                        f" refdb_corrected={value} Δ={value-bmrb_value:+.4f}"
                    )
            elif bmrb_lookup is not None and value is not None:
                # RefDB row exists but no matching BMRB row — unexpected
                # since RefDB is derived from BMRB. Flag loudly.
                entry["anomalies"].append(
                    f"REFDB_NO_MATCHING_BMRB_ROW: no BMRB row at"
                    f" (residue_seq={db_seq_raw}, atom={db_atomname})"
                )

        if bmrb_seq not in topology:
            entry["status"] = "UNMATCHED_ROW"
            rng = f"{min(topology)}..{max(topology)}"
            entry["anomalies"].append(
                f"RESIDUE_SEQ_OUT_OF_RANGE: residue_seq={bmrb_seq} not in"
                f" topology (range {rng})"
            )
            row_audit.append(entry)
            continue

        charmm_resname = topology[bmrb_seq]["name"]
        canonical = residue_aliases.get(charmm_resname, charmm_resname)

        if charmm_resname in residue_aliases:
            entry["anomalies"].append(
                f"RESIDUE_ALIAS_APPLIED: CHARMM '{charmm_resname}' → canonical"
                f" '{canonical}' (protonation-state alias; {db_label} does not"
                f" specify state)"
            )
        if canonical != db_resname:
            entry["anomalies"].append(
                f"RESIDUE_NAME_MISMATCH: {db_label}={db_resname}"
                f" CHARMM={charmm_resname} (canonical_alias={canonical})"
            )

        candidates = reverse[bmrb_seq].get(db_atomname, [])
        if not candidates:
            entry["status"] = "UNMATCHED_ROW"
            entry["anomalies"].append(
                f"UNMATCHED_NAME: {db_label} atom '{db_atomname}' has no CHARMM"
                f" counterpart in residue {bmrb_seq} ({charmm_resname}) under"
                f" the per-protein translation table"
            )
            row_audit.append(entry)
            continue

        elements = {topology[bmrb_seq]["atoms"][c] for c in candidates}
        if len(elements) > 1:
            entry["anomalies"].append(
                f"ELEMENT_HETEROGENEOUS_CANDIDATES: {sorted(elements)}"
            )
        charmm_element = next(iter(elements))
        if charmm_element != db_nucleus:
            entry["anomalies"].append(
                f"NUCLEUS_MISMATCH: {db_label} atom_type={db_nucleus}"
                f" CHARMM element={charmm_element}"
            )

        if len(candidates) > 1:
            entry["status"] = "ONE_TO_MANY"
            entry["matched_charmm"] = sorted(candidates)
            entry["anomalies"].append(
                f"ONE_TO_MANY: {db_label} '{db_atomname}' → {len(candidates)}"
                f" CHARMM atoms {sorted(candidates)} (pseudo-atom, expansion"
                f" policy deferred)"
            )
            for c in candidates:
                matched_charmm.add((bmrb_seq, c))
        else:
            charmm_atom = candidates[0]
            entry["matched_charmm"] = charmm_atom
            matched_charmm.add((bmrb_seq, charmm_atom))
            if charmm_atom != db_atomname:
                entry["anomalies"].append(
                    f"NAME_TRANSLATION_APPLIED: {db_label} '{db_atomname}' →"
                    f" CHARMM '{charmm_atom}'"
                )
            if ambig is not None and ambig >= 2:
                entry["status"] = "STEREO_AMBIGUOUS"
                entry["anomalies"].append(
                    f"STEREO_AMBIGUOUS: ambiguity_code={ambig}"
                )
            elif entry["anomalies"]:
                entry["status"] = "MATCHED_WITH_ANOMALIES"
            else:
                entry["status"] = "MATCHED"

        row_audit.append(entry)

    return row_audit, matched_charmm, offset


# ----------------------------------------------------------------------------
# Top-level audit
# ----------------------------------------------------------------------------

def audit(protein_dir, bmrb_id):
    pdb_path = Path(protein_dir) / "reference.pdb"
    if not pdb_path.is_file():
        raise FileNotFoundError(f"reference.pdb not found at {pdb_path}")

    protein_tag = f"{Path(protein_dir).name.split('_')[0]}_{bmrb_id}"
    bmrb_table_path = HERE / f"{protein_tag}.bmrb.toml"
    refdb_table_path = HERE / f"{protein_tag}.refdb.toml"

    missing = [p for p in (bmrb_table_path, refdb_table_path) if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            f"Per-protein translation tables missing: {[str(p) for p in missing]}."
            f"\nDerive each from _starter_translation.toml, add a provenance"
            f" block, and rerun. No silent fallback is permitted (2026-04-19)."
        )

    bmrb_translation = load_translation(bmrb_table_path)
    refdb_translation = load_translation(refdb_table_path)
    validate_provenance(bmrb_translation, bmrb_table_path)
    validate_provenance(refdb_translation, refdb_table_path)
    bmrb_aliases = bmrb_translation.get("residue_name_aliases", {})
    refdb_aliases = refdb_translation.get("residue_name_aliases", {})

    topology = parse_pdb(pdb_path)

    # Per-protein BMRB + RefDB CSVs produced upstream by fetch_bmrb.py,
    # extract_refdb.py, and parse_shifts.py.
    source_dir = SOURCE_DATA / protein_tag
    bmrb_csv = source_dir / f"bmrb_{bmrb_id}.csv"
    refdb_csv = source_dir / f"refdb_{bmrb_id}.csv"
    for p, stage in ((bmrb_csv, "fetch_bmrb.py + parse_shifts.py"),
                     (refdb_csv, "extract_refdb.py + parse_shifts.py")):
        if not p.is_file():
            raise FileNotFoundError(
                f"Per-protein CSV missing: {p}. Run {stage} first."
            )
    bmrb_rows = load_csv_rows(bmrb_csv)
    refdb_rows = load_csv_rows(refdb_csv)

    bmrb_audit, bmrb_matched, bmrb_offset = audit_rows(
        bmrb_rows, topology, bmrb_translation, bmrb_aliases,
        db_label="BMRB", value_col="shift_value"
    )
    # Cross-db compare: REFDB_CORRECTED_VALUE_DIFFERS now fires when the
    # RefDB shift_value disagrees with the BMRB shift_value at the same
    # (residue_seq, atom_name). This is the proper semantic (BMRB raw vs
    # RefDB re-referenced), not the old intra-row original-vs-corrected
    # comparison that the previous consolidated CSV required.
    bmrb_lookup = build_value_lookup(bmrb_rows)
    refdb_audit, refdb_matched, refdb_offset = audit_rows(
        refdb_rows, topology, refdb_translation, refdb_aliases,
        db_label="RefDB", value_col="shift_value",
        include_refdb_delta=True, bmrb_lookup=bmrb_lookup,
    )

    # CHARMM side — recorded once, annotated from both audits.
    charmm_audit = []
    for res_seq in sorted(topology):
        info = topology[res_seq]
        for atom, element in info["atoms"].items():
            charmm_audit.append({
                "residue_seq": res_seq,
                "residue_name": info["name"],
                "atom_name": atom,
                "element": element,
                "matched_by_bmrb": (res_seq, atom) in bmrb_matched,
                "matched_by_refdb": (res_seq, atom) in refdb_matched,
            })

    # Per-db summaries
    def summarize(name, row_audit, matched_set, offset):
        status_counts = defaultdict(int)
        anomaly_counts = defaultdict(int)
        for e in row_audit:
            status_counts[e["status"] or "NONE"] += 1
            for a in e["anomalies"]:
                anomaly_counts[a.split(":")[0]] += 1
        by_el = defaultdict(lambda: {"topology": 0, "matched": 0})
        for c in charmm_audit:
            by_el[c["element"]]["topology"] += 1
            if (c["residue_seq"], c["atom_name"]) in matched_set:
                by_el[c["element"]]["matched"] += 1
        return {
            "row_count": len(row_audit),
            "residue_seq_offset_applied": offset,
            "status_counts": dict(status_counts),
            "anomaly_counts": dict(anomaly_counts),
            "charmm_coverage_by_element": {el: dict(v) for el, v in by_el.items()},
        }

    summary = {
        "bmrb_id": bmrb_id,
        "protein_dir": str(protein_dir),
        "translation_tables": {
            "bmrb": str(bmrb_table_path.relative_to(HERE)),
            "refdb": str(refdb_table_path.relative_to(HERE)),
        },
        "topology": {
            "residue_count": len(topology),
            "charmm_atom_count": sum(len(r["atoms"]) for r in topology.values()),
        },
        "bmrb": summarize("bmrb", bmrb_audit, bmrb_matched, bmrb_offset),
        "refdb": summarize("refdb", refdb_audit, refdb_matched, refdb_offset),
    }

    return {
        "summary": summary,
        "bmrb_audit": bmrb_audit,
        "refdb_audit": refdb_audit,
        "charmm_audit": charmm_audit,
    }


# ----------------------------------------------------------------------------
# Markdown rendering
# ----------------------------------------------------------------------------

def _row_table(db_label, rows):
    lines = [
        f"## Every {db_label} row (one per line, in source order)",
        "",
        f"| Seq ({db_label}) | CHARMM seq | Res | Atom | Nuc | Value | Ambig | Status | CHARMM atom | Anomalies |",
        f"|---|---|---|---|---|---|---|---|---|---|",
    ]
    for e in rows:
        r = e["row"]
        m = e["matched_charmm"]
        m_cell = "{" + ",".join(m) + "}" if isinstance(m, list) else (m or "")
        anoms = "<br>".join(e["anomalies"])
        v = r.get("value", "")
        if "delta" in r:
            v = f"{r.get('original_value')}→{r.get('corrected_value')}" if r.get("delta") else f"{r.get('corrected_value')}"
        amb = r["ambiguity_code"] if r["ambiguity_code"] is not None else ""
        charmm_seq_cell = r.get("charmm_seq", "") if "charmm_seq" in r else r["residue_seq"]
        lines.append(
            f"| {r['residue_seq']} | {charmm_seq_cell} | {r['residue_name']} |"
            f" {r['atom_name']} | {r['atom_type']} | {v} | {amb} |"
            f" `{e['status']}` | {m_cell} | {anoms} |"
        )
    lines.append("")
    return lines


def write_markdown(result, md_path):
    s = result["summary"]
    lines = []
    lines.append(f"# NMR Audit — BMRB {s['bmrb_id']}")
    lines.append("")
    lines.append(f"- Protein dir: `{s['protein_dir']}`")
    lines.append(f"- CHARMM topology: **{s['topology']['residue_count']} residues**,"
                 f" **{s['topology']['charmm_atom_count']} atoms**")
    lines.append(f"- Translation tables:")
    lines.append(f"  - BMRB: `{s['translation_tables']['bmrb']}`")
    lines.append(f"  - RefDB: `{s['translation_tables']['refdb']}`")
    lines.append("")

    for db in ("bmrb", "refdb"):
        label = db.upper()
        ds = s[db]
        lines.append(f"## {label} summary")
        lines.append("")
        lines.append(f"- Row count: **{ds['row_count']}**")
        offset = ds.get("residue_seq_offset_applied", 0)
        if offset:
            lines.append(
                f"- **residue_seq offset applied: {offset:+d}** (every"
                f" {label} row's residue_seq was shifted by {offset:+d}"
                f" before topology lookup per the per-protein table's"
                f" `[residue_seq_offset]` declaration)"
            )
        lines.append(f"- Status:")
        for k, v in sorted(ds["status_counts"].items()):
            lines.append(f"  - `{k}` — {v}")
        lines.append(f"- Anomalies:")
        if not ds["anomaly_counts"]:
            lines.append("  - _(none)_")
        else:
            for k, v in sorted(ds["anomaly_counts"].items(), key=lambda kv: (-kv[1], kv[0])):
                lines.append(f"  - `{k}` — {v}")
        lines.append(f"- CHARMM coverage:")
        for el, cov in sorted(ds["charmm_coverage_by_element"].items()):
            lines.append(f"  - {el}: {cov['matched']}/{cov['topology']}")
        lines.append("")

    lines.extend(_row_table("BMRB", result["bmrb_audit"]))
    lines.extend(_row_table("RefDB", result["refdb_audit"]))

    lines.append("## CHARMM atoms — matched by BMRB and/or RefDB")
    lines.append("")
    lines.append("| Seq | Res | Atom | El | BMRB | RefDB |")
    lines.append("|---|---|---|---|---|---|")
    for c in result["charmm_audit"]:
        lines.append(
            f"| {c['residue_seq']} | {c['residue_name']} | {c['atom_name']} |"
            f" {c['element']} | {'✓' if c['matched_by_bmrb'] else '—'} |"
            f" {'✓' if c['matched_by_refdb'] else '—'} |"
        )
    lines.append("")

    md_path.write_text("\n".join(lines))


# ----------------------------------------------------------------------------
# Entry
# ----------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print("Usage: audit_nmr.py <protein_dir> <bmrb_id>", file=sys.stderr)
        sys.exit(1)
    protein_dir = Path(sys.argv[1])
    bmrb_id = int(sys.argv[2])

    result = audit(protein_dir, bmrb_id)
    pdb_id = protein_dir.name.split("_")[0]
    json_path = HERE / f"{pdb_id}_{bmrb_id}.audit.json"
    md_path = HERE / f"{pdb_id}_{bmrb_id}.audit.md"
    with json_path.open("w") as f:
        json.dump(result, f, indent=2)
    write_markdown(result, md_path)

    print(f"Wrote {json_path} and {md_path}")
    print()
    print("Summary:")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
