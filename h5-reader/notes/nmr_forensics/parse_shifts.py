#!/usr/bin/env python3
"""
parse_shifts.py — parse per-protein BMRB and RefDB NMR-STAR files to CSV.

For each of the 10 calibration proteins, reads:
  source_data/<pdb>_<bmrb>/raw/bmr<id>_3.str       (BMRB, NMR-STAR 3.x)
  source_data/<pdb>_<bmrb>/raw/bmr<id>.str.corr    (RefDB, NMR-STAR 2.1.1)

Writes:
  source_data/<pdb>_<bmrb>/bmrb_<id>.csv           — every BMRB _Atom_chem_shift row
  source_data/<pdb>_<bmrb>/refdb_<id>.csv          — every RefDB _Atom_shift_assign_ID row
  source_data/<pdb>_<bmrb>/refdb_<id>_header.toml  — RefDB offset corrections + outliers

No filtering. No value transformations. All fields preserved as-deposited.

Run from nmr_forensics/ in the venv:
    ./.venv/bin/python parse_shifts.py
"""

import csv
import re
import sys
from pathlib import Path

import pynmrstar

HERE = Path(__file__).parent
SOURCE_DATA = HERE / "source_data"

PROTEINS = [
    ("1DV0", 4757),
    ("1HS5", 4934),
    ("1HD6", 4820),
    ("1B1V", 4292),
    ("1HA9", 5028),
    ("1G26", 4656),
    ("1CBH", 192),
    ("1I2V", 4976),
    ("1I8X", 4351),
    ("1BWX", 3449),
]

# ----------------------------------------------------------------------------
# BMRB NMR-STAR 3.x parsing (via pynmrstar)
# ----------------------------------------------------------------------------

BMRB_CANONICAL_COLUMNS = [
    ("source_row_id", "ID"),
    ("residue_seq", "Comp_index_ID"),
    ("residue_name", "Comp_ID"),
    ("atom_name", "Atom_ID"),
    ("atom_type", "Atom_type"),
    ("isotope", "Atom_isotope_number"),
    ("shift_value", "Val"),
    ("shift_error", "Val_err"),
    ("ambiguity_code", "Ambiguity_code"),
    ("auth_residue_seq", "Auth_seq_ID"),
    ("auth_residue_name", "Auth_comp_ID"),
    ("auth_atom_name", "Auth_atom_ID"),
    ("details", "Details"),
]


def parse_bmrb(tag, bmrb_id):
    star_path = SOURCE_DATA / tag / "raw" / f"bmr{bmrb_id}_3.str"
    entry = pynmrstar.Entry.from_file(str(star_path))

    saveframes = entry.get_saveframes_by_category("assigned_chemical_shifts")
    if not saveframes:
        raise RuntimeError(f"[{tag}] no assigned_chemical_shifts saveframe")
    if len(saveframes) > 1:
        # Document and continue with all of them; should be 1 for our 10 but
        # don't silently drop if future data has multiple.
        print(
            f"[{tag}] WARNING: {len(saveframes)} assigned_chemical_shifts saveframes; "
            f"concatenating all in output CSV",
            file=sys.stderr,
        )

    all_rows = []
    for sf in saveframes:
        for loop in sf:
            if loop.category != "_Atom_chem_shift":
                continue
            tag_index = {name: loop.tags.index(name) for name in loop.tags}
            for data_row in loop.data:
                out_row = {}
                for out_name, src_name in BMRB_CANONICAL_COLUMNS:
                    if src_name in tag_index:
                        out_row[out_name] = data_row[tag_index[src_name]]
                    else:
                        out_row[out_name] = ""
                # Also tag which saveframe it came from for audit
                out_row["saveframe"] = sf.name
                all_rows.append(out_row)

    csv_path = SOURCE_DATA / tag / f"bmrb_{bmrb_id}.csv"
    with csv_path.open("w", newline="") as f:
        fieldnames = [n for n, _ in BMRB_CANONICAL_COLUMNS] + ["saveframe"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)

    return len(all_rows)


# ----------------------------------------------------------------------------
# RefDB NMR-STAR 2.1.1 parsing (hand-roll, since pynmrstar 3.x refuses 2.1)
# ----------------------------------------------------------------------------

REFDB_CANONICAL_COLUMNS = [
    ("source_row_id", "_Atom_shift_assign_ID"),
    ("residue_seq", "_Residue_seq_code"),
    ("residue_name", "_Residue_label"),
    ("atom_name", "_Atom_name"),
    ("atom_type", "_Atom_type"),
    ("shift_value", "_Chem_shift_value"),
    ("shift_error", "_Chem_shift_value_error"),
    ("ambiguity_code", "_Chem_shift_ambiguity_code"),
]


def parse_refdb_loop(text):
    """
    Find the loop whose first tag is _Atom_shift_assign_ID and return its
    rows as list of dicts keyed by out_name from REFDB_CANONICAL_COLUMNS.
    """
    # Find the loop_ ... stop_ block where the first tag is _Atom_shift_assign_ID.
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            tags = []
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if "_Atom_shift_assign_ID" in tags:
                # Read data rows until `stop_`
                data_rows = []
                k = j
                while k < len(lines) and lines[k].strip() != "stop_":
                    line = lines[k].strip()
                    if line and not line.startswith("#"):
                        # tokens are whitespace-separated; no quoted strings in
                        # the shift loop for RefDB 2.1.1
                        toks = line.split()
                        if len(toks) == len(tags):
                            data_rows.append(toks)
                    k += 1
                return tags, data_rows
            else:
                i = j
                continue
        i += 1
    raise RuntimeError("No _Atom_shift_assign_ID loop found")


def parse_refdb_header(text):
    """
    Extract the per-nucleus offset corrections, average differences, and
    persistent outliers from the RefDB file's `#` comment header.
    """
    # The structure is:
    #   #bmr<id>.str.corr chemical shifts have been re-referenced with the
    #   #following offsets (these values have been added to the original
    #   #bmr<id>.str file):
    #   #HA       CA       CB      CO       N      HN
    #   #N/A    -0.11   -0.11    N/A    -0.72     N/A
    #
    # Also captures the "average CS difference" table and the persistent-outlier
    # tables. Stored as structured TOML under sections.

    reference_pdb = None
    m = re.search(r"#Corrected using PDB structure:\s*(\S+)", text)
    if m:
        reference_pdb = m.group(1)

    offsets_applied = None
    m = re.search(
        r"#bmr\d+\.str\.corr chemical shifts have been re-referenced with the following\s*\n"
        r"#offsets .*?\n"
        r"#(HA\s+CA\s+CB\s+CO\s+N\s+HN)\s*\n"
        r"#(.*?)\s*\n",
        text,
        re.DOTALL,
    )
    if m:
        headers = m.group(1).split()
        values = m.group(2).split()
        offsets_applied = dict(zip(headers, values))

    avg_diff = None
    m = re.search(
        r"#The average CS difference between predicted and observed:\s*\n"
        r"#(HA\s+CA\s+CB\s+CO\s+N\s+HN)\s*\n"
        r"#(.*?)\s*\n",
        text,
        re.DOTALL,
    )
    if m:
        headers = m.group(1).split()
        values = m.group(2).split()
        avg_diff = dict(zip(headers, values))

    # Outlier tables (residues still showing large deltas after correction).
    # Each block starts with "#After reference correction, the following
    # residues still have a <NUC> chemical shift difference (obs*-pred) greater
    # than <THR>ppm:" followed by "#NUM AA CS Observed* Predicted" and rows.
    outlier_blocks = []
    pattern = re.compile(
        r"#After reference correction, the following residues still\s*\n"
        r"#have a (\S+) chemical shift difference \(obs\*-pred\) greater than\s*(\S+)ppm:\s*\n"
        r"#NUM\s+AA\s+CS\s+Observed\*\s+Predicted\s*\n"
        r"((?:#\s*\d+.*?\n)+)",
    )
    for nuc, thr, rows_text in pattern.findall(text):
        rows = []
        for row_line in rows_text.strip().split("\n"):
            toks = row_line.lstrip("#").strip().split()
            if len(toks) >= 5:
                rows.append(
                    {
                        "residue_seq": toks[0],
                        "residue_name_1letter": toks[1],
                        "atom_name": toks[2],
                        "observed_corrected": toks[3],
                        "predicted": toks[4],
                    }
                )
        outlier_blocks.append(
            {"nucleus": nuc, "threshold_ppm": thr, "rows": rows}
        )

    return {
        "reference_pdb": reference_pdb,
        "offsets_applied": offsets_applied,
        "avg_diff": avg_diff,
        "outlier_blocks": outlier_blocks,
    }


def write_refdb_header_toml(tag, bmrb_id, header):
    path = SOURCE_DATA / tag / f"refdb_{bmrb_id}_header.toml"
    lines = [
        "# RefDB .str.corr file header metadata (parsed from the # comment block).",
        "# These fields describe what SHIFTCOR/SHIFTX applied to the original BMRB",
        "# deposition to produce the re-referenced values in refdb_<id>.csv.",
        "",
        f'bmrb_id = {bmrb_id}',
    ]
    if header["reference_pdb"]:
        lines.append(f'reference_pdb = "{header["reference_pdb"]}"')
    lines.append("")

    if header["offsets_applied"]:
        lines.append("[offsets_applied]")
        lines.append("# Offsets that were ADDED to original BMRB shifts to produce RefDB values.")
        lines.append("# Values: number in ppm, or 'N/A' if no correction applied to that nucleus.")
        for k, v in header["offsets_applied"].items():
            lines.append(f'{k} = "{v}"')
        lines.append("")

    if header["avg_diff"]:
        lines.append("[avg_difference_obs_minus_pred]")
        lines.append("# Average chemical-shift difference between corrected observation and")
        lines.append("# SHIFTX prediction, per nucleus (ppm).")
        for k, v in header["avg_diff"].items():
            lines.append(f'{k} = "{v}"')
        lines.append("")

    for block in header["outlier_blocks"]:
        lines.append("[[persistent_outliers]]")
        lines.append(f'nucleus = "{block["nucleus"]}"')
        lines.append(f'threshold_ppm = "{block["threshold_ppm"]}"')
        for r in block["rows"]:
            lines.append("[[persistent_outliers.rows]]")
            for k, v in r.items():
                lines.append(f'{k} = "{v}"')
        lines.append("")

    path.write_text("\n".join(lines))


def parse_refdb(tag, bmrb_id):
    corr_path = SOURCE_DATA / tag / "raw" / f"bmr{bmrb_id}.str.corr"
    text = corr_path.read_text(errors="replace")

    tags, data_rows = parse_refdb_loop(text)
    tag_index = {t: tags.index(t) for t in tags}

    out_rows = []
    for row in data_rows:
        out_row = {}
        for out_name, src_name in REFDB_CANONICAL_COLUMNS:
            if src_name in tag_index:
                out_row[out_name] = row[tag_index[src_name]]
            else:
                out_row[out_name] = ""
        out_rows.append(out_row)

    csv_path = SOURCE_DATA / tag / f"refdb_{bmrb_id}.csv"
    with csv_path.open("w", newline="") as f:
        fieldnames = [n for n, _ in REFDB_CANONICAL_COLUMNS]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    header = parse_refdb_header(text)
    write_refdb_header_toml(tag, bmrb_id, header)

    return len(out_rows)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    rows = []
    for pdb_id, bmrb_id in PROTEINS:
        tag = f"{pdb_id}_{bmrb_id}"
        bmrb_n = parse_bmrb(tag, bmrb_id)
        refdb_n = parse_refdb(tag, bmrb_id)
        rows.append((tag, bmrb_n, refdb_n))
        print(f"[{tag}] BMRB {bmrb_n} rows, RefDB {refdb_n} rows")

    print("\nSummary:")
    print(f"  {'Protein':14s}  {'BMRB':>6s}  {'RefDB':>6s}")
    for tag, bmrb_n, refdb_n in rows:
        print(f"  {tag:14s}  {bmrb_n:6d}  {refdb_n:6d}")


if __name__ == "__main__":
    main()
