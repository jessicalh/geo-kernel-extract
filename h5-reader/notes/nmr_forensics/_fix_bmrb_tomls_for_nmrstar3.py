#!/usr/bin/env python3
"""
One-shot fix (2026-04-20): BMRB-side translation tables were authored against
the consolidated CSV, which turned out to be RefDB-style (pseudo-atom names
for methyls: ALA HB, ILE HD1, THR HG2, etc.). BMRB 3.x NMR-STAR files
(what we now fetch directly) use **expanded** methyl names (HB1/HB2/HB3,
HD11/HD12/HD13, HG21/HG22/HG23, etc.). This script updates the 10 BMRB
per-protein TOMLs to match BMRB 3.x conventions.

Changes per TOML (BMRB side only — RefDB TOMLs unchanged):
  - Remove methyl-pseudo fanning for ALA, VAL, LEU (δ methyls), ILE
    (γ2+δ methyls), MET, THR, LYS (NH3+ ζ).
  - Keep methylene translations (HB1→HB2, HB2→HB3, etc.).
  - Keep ILE CD→CD1 (no-locant vs locant divergence in CARBON naming).
  - Keep ILE γ1-methylene HG11→HG12, HG12→HG13.
  - Update ILE δ-methyl: CHARMM HD1/HD2/HD3 → BMRB HD11/HD12/HD13
    (not the pseudo "HD1").
  - GLY HA1/HA2 → HA2/HA3 kept (BMRB uses 2/3 convention on methylenes).

After running, re-run audit_nmr.py to confirm UNMATCHED rates drop to
the expected N-terminal residual only.
"""

from pathlib import Path

HERE = Path(__file__).parent

# Map of (residue_section_name → new_content_lines). Each entry is a list
# of lines that should appear under `[residue.X]`. If the entry is an
# empty list, the section is removed entirely (since BMRB 3.x uses
# identity for everything in that residue).

BMRB_RESIDUE_BLOCKS = {
    "GLY": ['HA1 = "HA2"', 'HA2 = "HA3"'],
    "ALA": [],  # identity for HB1/HB2/HB3
    "VAL": [],  # identity for γ1/γ2 methyls
    "LEU": ['HB1 = "HB2"', 'HB2 = "HB3"'],  # methylene only
    "ILE": [
        "# BMRB 3.x uses expanded methyl names (HG21/HG22/HG23,",
        "# HD11/HD12/HD13); only methylene renumbering and carbon-locant",
        "# need mapping.",
        'CD = "CD1"',
        "# γ1-methylene: CHARMM HG11/HG12 → BMRB HG12/HG13",
        'HG11 = "HG12"',
        'HG12 = "HG13"',
        "# δ-methyl: CHARMM uses HD1/HD2/HD3 (1-based), BMRB uses HD11/HD12/HD13",
        'HD1 = "HD11"',
        'HD2 = "HD12"',
        'HD3 = "HD13"',
    ],
    "MET": ['HB1 = "HB2"', 'HB2 = "HB3"', 'HG1 = "HG2"', 'HG2 = "HG3"'],
    "PRO": ['HB1 = "HB2"', 'HB2 = "HB3"', 'HG1 = "HG2"', 'HG2 = "HG3"',
            'HD1 = "HD2"', 'HD2 = "HD3"'],
    "PHE": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "TRP": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "TYR": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "SER": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "THR": [],  # γ-methyl identity HG21/HG22/HG23
    "CYS": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "ASP": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "ASN": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "GLU": ['HB1 = "HB2"', 'HB2 = "HB3"', 'HG1 = "HG2"', 'HG2 = "HG3"'],
    "GLN": ['HB1 = "HB2"', 'HB2 = "HB3"', 'HG1 = "HG2"', 'HG2 = "HG3"'],
    "LYS": ['HB1 = "HB2"', 'HB2 = "HB3"', 'HG1 = "HG2"', 'HG2 = "HG3"',
            'HD1 = "HD2"', 'HD2 = "HD3"', 'HE1 = "HE2"', 'HE2 = "HE3"'],
    # LYS NH3+ ζ protons: not translated — BMRB typically omits (fast exchange).
    "ARG": ['HB1 = "HB2"', 'HB2 = "HB3"', 'HG1 = "HG2"', 'HG2 = "HG3"',
            'HD1 = "HD2"', 'HD2 = "HD3"'],
    "HSP": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "HSD": ['HB1 = "HB2"', 'HB2 = "HB3"'],
    "HSE": ['HB1 = "HB2"', 'HB2 = "HB3"'],
}

PROTEINS = [
    "1DV0_4757", "1HS5_4934", "1HD6_4820", "1B1V_4292", "1HA9_5028",
    "1G26_4656", "1CBH_192", "1I2V_4976", "1I8X_4351", "1BWX_3449",
]


def rewrite_one(path):
    text = path.read_text()

    # Split the file into:
    #   header = everything before the first `[residue.` section
    #   residue_blocks = the [residue.X] sections (to be rewritten)
    #   tail = the [residue_name_aliases] section and anything after

    lines = text.split("\n")

    # Find section boundaries
    first_residue_idx = None
    aliases_idx = None
    for i, line in enumerate(lines):
        if first_residue_idx is None and line.startswith("[residue."):
            first_residue_idx = i
        if line.strip() == "[residue_name_aliases]":
            aliases_idx = i
            break

    if first_residue_idx is None or aliases_idx is None:
        print(f"  SKIP {path.name}: could not locate section boundaries")
        return

    header = "\n".join(lines[:first_residue_idx]).rstrip()
    tail = "\n".join(lines[aliases_idx:])

    # Rebuild residue blocks per BMRB_RESIDUE_BLOCKS.
    # Special note for HIS blocks: only include if the protein's TOML had
    # them (i.e., only if the protein contains HIS). Check by looking at
    # whether the original tail contained HSP/HSD/HSE references. For
    # simplicity, always emit all blocks; the audit ignores residues not
    # in topology.

    residue_body = []
    for res, entries in BMRB_RESIDUE_BLOCKS.items():
        if not entries:
            continue  # no section needed; all identity
        residue_body.append(f"[residue.{res}]")
        for e in entries:
            residue_body.append(e)
        residue_body.append("")

    new_text = header + "\n\n" + "\n".join(residue_body) + tail + ("\n" if not text.endswith("\n") else "")
    path.write_text(new_text)
    print(f"  rewrote {path.name}")


def main():
    for tag in PROTEINS:
        path = HERE / f"{tag}.bmrb.toml"
        if not path.is_file():
            print(f"MISSING {path.name}")
            continue
        rewrite_one(path)


if __name__ == "__main__":
    main()
