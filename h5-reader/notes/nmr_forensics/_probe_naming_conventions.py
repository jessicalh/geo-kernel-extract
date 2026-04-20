#!/usr/bin/env python3
"""
One-shot empirical probe: what naming convention does each BMRB / RefDB
deposition actually use? Answers the question 'is every entry IUPAC, or
did they each just try?'

Reads the already-parsed per-protein CSVs under source_data/<tag>/ and
tabulates, per deposition, which form each structural class uses:

  1. ALA β-methyl                  pseudo 'HB'        vs. expanded HB1/HB2/HB3
  2. THR γ-methyl                  pseudo 'HG2'       vs. expanded HG21/HG22/HG23
  3. LEU δ-methyls                 pseudo HD1/HD2     vs. expanded HD11-HD13 / HD21-HD23
  4. ILE δ-methyl                  pseudo 'HD1'       vs. CHARMM HD1/HD2/HD3       vs. IUPAC HD11/HD12/HD13
  5. ILE γ-carbon label            'CD'               vs. 'CD1'
  6. β-methylene (any residue)     HB1/HB2 (pre-1998) vs. HB2/HB3 (IUPAC-1998)
  7. LYS ζ-NH3+                    pseudo 'HZ'        vs. expanded HZ1/HZ2/HZ3

Output: one row per deposition, 7 columns. No source-of-truth claim — just
what the deposition deposits. Decisions about convention naming live with
the reviewer.
"""

import csv
from pathlib import Path

HERE = Path(__file__).parent
SOURCE_DATA = HERE / "source_data"

PROTEINS = [
    ("1DV0", 4757), ("1HS5", 4934), ("1HD6", 4820), ("1B1V", 4292),
    ("1HA9", 5028), ("1G26", 4656), ("1CBH", 192),   ("1I2V", 4976),
    ("1I8X", 4351), ("1BWX", 3449),
]


def load(csv_path):
    rows = []
    with csv_path.open() as f:
        r = csv.DictReader(f)
        rows.extend(r)
    return rows


def classify(rows):
    names_by_res = {}
    for r in rows:
        names_by_res.setdefault(r["residue_name"], set()).add(r["atom_name"])

    def has(res, atom):
        return atom in names_by_res.get(res, set())

    results = {}

    # 1. ALA β-methyl
    if "ALA" in names_by_res:
        pseudo = has("ALA", "HB")
        expanded = any(has("ALA", f"HB{i}") for i in (1, 2, 3))
        results["ALA_HB"] = "pseudo" if pseudo and not expanded else (
            "expanded" if expanded and not pseudo else (
                "both" if pseudo and expanded else "absent"
            )
        )
    else:
        results["ALA_HB"] = "no-ALA"

    # 2. THR γ-methyl
    if "THR" in names_by_res:
        pseudo = has("THR", "HG2")
        expanded = any(has("THR", f"HG2{i}") for i in (1, 2, 3))
        results["THR_HG"] = ("pseudo" if pseudo and not expanded else
                              ("expanded" if expanded and not pseudo else
                               ("both" if pseudo and expanded else "absent")))
    else:
        results["THR_HG"] = "no-THR"

    # 3. LEU δ-methyls
    if "LEU" in names_by_res:
        pseudo = has("LEU", "HD1") or has("LEU", "HD2")
        expanded_1 = any(has("LEU", f"HD1{i}") for i in (1, 2, 3))
        expanded_2 = any(has("LEU", f"HD2{i}") for i in (1, 2, 3))
        expanded = expanded_1 or expanded_2
        results["LEU_HD"] = ("pseudo" if pseudo and not expanded else
                              ("expanded" if expanded and not pseudo else
                               ("both" if pseudo and expanded else "absent")))
    else:
        results["LEU_HD"] = "no-LEU"

    # 4. ILE δ-methyl
    if "ILE" in names_by_res:
        pseudo = has("ILE", "HD1") and not has("ILE", "HD2") and not has("ILE", "HD3")
        charmm = has("ILE", "HD1") and has("ILE", "HD2") and has("ILE", "HD3")
        iupac = any(has("ILE", f"HD1{i}") for i in (1, 2, 3))
        kinds = []
        if pseudo and not iupac and not (has("ILE", "HD2") or has("ILE", "HD3")):
            kinds.append("pseudo-HD1")
        if charmm:
            kinds.append("charmm-HD1/HD2/HD3")
        if iupac:
            kinds.append("iupac-HD11/HD12/HD13")
        results["ILE_HDm"] = ",".join(kinds) if kinds else "absent"
    else:
        results["ILE_HDm"] = "no-ILE"

    # 5. ILE γ-carbon locant
    if "ILE" in names_by_res:
        has_cd = has("ILE", "CD")
        has_cd1 = has("ILE", "CD1")
        results["ILE_Cd"] = ("CD" if has_cd and not has_cd1 else
                              ("CD1" if has_cd1 and not has_cd else
                               ("both" if has_cd and has_cd1 else "absent-or-H-only")))
    else:
        results["ILE_Cd"] = "no-ILE"

    # 6. β-methylene form (examine residues that have a β-methylene: not GLY,
    #    not ALA, not VAL, not THR, not ILE — keep it simple: ASP is a good probe)
    candidate_residues = ["ASP", "ASN", "SER", "CYS", "GLU", "GLN", "PHE",
                          "TYR", "HIS", "LEU", "LYS", "ARG", "MET", "TRP", "PRO"]
    any_res_with_methylene = next((r for r in candidate_residues if r in names_by_res), None)
    if any_res_with_methylene:
        r = any_res_with_methylene
        pre1998 = has(r, "HB1") and has(r, "HB2")
        post1998 = has(r, "HB2") and has(r, "HB3")
        if pre1998 and post1998:
            results["bMethyl"] = f"both-{r}"
        elif pre1998:
            results["bMethyl"] = f"pre1998-{r}(HB1/HB2)"
        elif post1998:
            results["bMethyl"] = f"iupac-{r}(HB2/HB3)"
        else:
            results["bMethyl"] = f"{r}-neither"
    else:
        results["bMethyl"] = "no-probe-residue"

    # 7. LYS ζ-NH3+
    if "LYS" in names_by_res:
        pseudo = has("LYS", "HZ")
        expanded = any(has("LYS", f"HZ{i}") for i in (1, 2, 3))
        if pseudo and expanded:
            results["LYS_HZ"] = "both"
        elif pseudo:
            results["LYS_HZ"] = "pseudo"
        elif expanded:
            results["LYS_HZ"] = "expanded"
        else:
            results["LYS_HZ"] = "absent"
    else:
        results["LYS_HZ"] = "no-LYS"

    return results


def main():
    print(f"{'Deposition':14s}  {'ALA-HB':10s} {'THR-HG':10s} {'LEU-HD':10s} "
          f"{'ILE-HDm':30s} {'ILE-Cd':15s} {'bMet':35s} {'LYS-HZ':10s}")
    for pdb, bmrb in PROTEINS:
        tag = f"{pdb}_{bmrb}"
        for db in ("bmrb", "refdb"):
            p = SOURCE_DATA / tag / f"{db}_{bmrb}.csv"
            rows = load(p)
            r = classify(rows)
            label = f"{tag}.{db}"
            print(f"{label:20s}  {r['ALA_HB']:10s} {r['THR_HG']:10s} {r['LEU_HD']:10s} "
                  f"{r['ILE_HDm']:30s} {r['ILE_Cd']:15s} {r['bMethyl']:35s} {r['LYS_HZ']:10s}")


if __name__ == "__main__":
    main()
