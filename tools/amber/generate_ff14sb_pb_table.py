#!/usr/bin/env python3
"""Generate the nmr-shielding ff14SB charge/PB-radius table.

The table is derived from AmberTools `tleap` with `leaprc.protein.ff14SB`
and `set default PBRadii mbondi2`. It records the force-field partial
charge and Amber PB radius for each residue/atom template used by the
flat `ParamFileChargeSource` path.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import subprocess
import tempfile
from dataclasses import dataclass


AMBER_CHARGE_FACTOR = 18.2223

STANDARD_UNITS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]

# Variant names from AminoAcidType. The AmberTools library files decide which
# of these have ff14SB templates; this list is the project-side filter.
MODEL_VARIANT_UNITS = [
    "ARN", "ASH", "CYM", "CYX", "GLH", "HID", "HIE", "HIP", "LYN", "TYM",
]


@dataclass(frozen=True)
class TemplateSpec:
    terminal_state: str
    residue_name: str
    sequence_units: tuple[str, ...]
    extract_residue_index: int  # 1-based Amber residue index


@dataclass(frozen=True)
class TableRow:
    terminal_state: str
    residue_name: str
    atom_name: str
    partial_charge: float
    pb_radius: float


def infer_amberhome(tleap: pathlib.Path, explicit: pathlib.Path | None) -> pathlib.Path:
    if explicit is not None:
        return explicit
    env = os.environ.get("AMBERHOME")
    if env:
        return pathlib.Path(env)
    return tleap.resolve().parent.parent


def read_unit_names(lib_path: pathlib.Path) -> set[str]:
    text = lib_path.read_text()
    return set(re.findall(r"^!entry\.([A-Za-z0-9]+)\.unit\.atoms table", text, re.M))


def amber_unit_sets(amberhome: pathlib.Path) -> tuple[set[str], set[str], set[str]]:
    lib_dir = amberhome / "dat" / "leap" / "lib"
    amino = read_unit_names(lib_dir / "amino12.lib")
    aminont = read_unit_names(lib_dir / "aminont12.lib")
    aminoct = read_unit_names(lib_dir / "aminoct12.lib")
    return amino, aminont, aminoct


def strip_terminal_prefix(units: set[str], prefix: str) -> set[str]:
    return {unit[1:] for unit in units if unit.startswith(prefix) and len(unit) > 1}


def parse_prmtop_strings(lines: list[str], width: int) -> list[str]:
    values: list[str] = []
    for line in lines:
        for pos in range(0, len(line), width):
            item = line[pos:pos + width].strip()
            if item:
                values.append(item)
    return values


def parse_prmtop_numbers(lines: list[str], typ: type[float] | type[int]) -> list:
    values = []
    for line in lines:
        for token in line.split():
            values.append(typ(token))
    return values


def read_prmtop_flag(path: pathlib.Path, flag: str) -> list[str]:
    lines = path.read_text().splitlines()
    out: list[str] = []
    in_flag = False
    saw_format = False
    for line in lines:
        if line.startswith("%FLAG "):
            if in_flag:
                break
            in_flag = line.strip() == f"%FLAG {flag}"
            saw_format = False
            continue
        if not in_flag:
            continue
        if line.startswith("%FORMAT"):
            saw_format = True
            continue
        if saw_format:
            out.append(line)
    if not out:
        raise RuntimeError(f"{path}: missing %FLAG {flag}")
    return out


def parse_prmtop(path: pathlib.Path) -> tuple[list[str], list[str], list[int], list[float], list[float]]:
    atom_names = parse_prmtop_strings(read_prmtop_flag(path, "ATOM_NAME"), 4)
    residue_labels = parse_prmtop_strings(read_prmtop_flag(path, "RESIDUE_LABEL"), 4)
    residue_pointers = parse_prmtop_numbers(read_prmtop_flag(path, "RESIDUE_POINTER"), int)
    charges = [
        value / AMBER_CHARGE_FACTOR
        for value in parse_prmtop_numbers(read_prmtop_flag(path, "CHARGE"), float)
    ]
    radii = parse_prmtop_numbers(read_prmtop_flag(path, "RADII"), float)
    return atom_names, residue_labels, residue_pointers, charges, radii


def run_tleap(tleap: pathlib.Path, spec: TemplateSpec, workdir: pathlib.Path) -> pathlib.Path:
    unit_expr = " ".join(spec.sequence_units)
    stem = f"{spec.terminal_state}_{spec.residue_name}"
    leap_in = workdir / f"{stem}.leap.in"
    prmtop = workdir / f"{stem}.prmtop"
    inpcrd = workdir / f"{stem}.inpcrd"
    leap_in.write_text(
        "\n".join([
            "source leaprc.protein.ff14SB",
            "set default PBRadii mbondi2",
            f"mol = sequence {{ {unit_expr} }}",
            f"saveamberparm mol {prmtop} {inpcrd}",
            "quit",
            "",
        ])
    )
    proc = subprocess.run(
        [str(tleap), "-f", str(leap_in)],
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if proc.returncode != 0 or not prmtop.exists():
        raise RuntimeError(
            f"tleap failed for {spec.terminal_state} {spec.residue_name}\n"
            f"sequence: {unit_expr}\n{proc.stdout}"
        )
    return prmtop


def rows_for_template(tleap: pathlib.Path, spec: TemplateSpec, workdir: pathlib.Path) -> list[TableRow]:
    prmtop = run_tleap(tleap, spec, workdir)
    atom_names, residue_labels, residue_pointers, charges, radii = parse_prmtop(prmtop)

    if len(charges) != len(atom_names) or len(radii) != len(atom_names):
        raise RuntimeError(f"{prmtop}: inconsistent atom/charge/radius counts")
    if spec.extract_residue_index < 1 or spec.extract_residue_index > len(residue_labels):
        raise RuntimeError(f"{prmtop}: residue index {spec.extract_residue_index} unavailable")

    start = residue_pointers[spec.extract_residue_index - 1] - 1
    if spec.extract_residue_index == len(residue_pointers):
        stop = len(atom_names)
    else:
        stop = residue_pointers[spec.extract_residue_index] - 1

    rows: list[TableRow] = []
    for idx in range(start, stop):
        rows.append(TableRow(
            terminal_state=spec.terminal_state,
            residue_name=spec.residue_name,
            atom_name=atom_names[idx],
            partial_charge=charges[idx],
            pb_radius=radii[idx],
        ))
    return rows


def build_specs(amberhome: pathlib.Path) -> list[TemplateSpec]:
    specs: list[TemplateSpec] = []
    internal_units, nterm_units, cterm_units = amber_unit_sets(amberhome)

    supported_internal_variants = [
        unit for unit in MODEL_VARIANT_UNITS if unit in internal_units
    ]
    supported_terminal_variants = [
        unit for unit in MODEL_VARIANT_UNITS
        if f"N{unit}" in nterm_units and f"C{unit}" in cterm_units
    ]

    for unit in STANDARD_UNITS + supported_internal_variants:
        specs.append(TemplateSpec(
            terminal_state="INTERNAL",
            residue_name=unit,
            sequence_units=("ALA", unit, "ALA"),
            extract_residue_index=2,
        ))

    for unit in STANDARD_UNITS + supported_terminal_variants:
        specs.append(TemplateSpec(
            terminal_state="NTERM",
            residue_name=unit,
            sequence_units=(f"N{unit}", "ALA"),
            extract_residue_index=1,
        ))
        specs.append(TemplateSpec(
            terminal_state="CTERM",
            residue_name=unit,
            sequence_units=("ALA", f"C{unit}"),
            extract_residue_index=2,
        ))

    return specs


def write_table(path: pathlib.Path, rows: list[TableRow], amberhome: pathlib.Path) -> None:
    internal_units, nterm_units, cterm_units = amber_unit_sets(amberhome)
    supported_internal_variants = [
        unit for unit in MODEL_VARIANT_UNITS if unit in internal_units
    ]
    unsupported_internal_variants = [
        unit for unit in MODEL_VARIANT_UNITS if unit not in internal_units
    ]
    supported_terminal_variants = [
        unit for unit in MODEL_VARIANT_UNITS
        if f"N{unit}" in nterm_units and f"C{unit}" in cterm_units
    ]
    unsupported_terminal_variants = [
        unit for unit in MODEL_VARIANT_UNITS
        if unit not in supported_terminal_variants
    ]

    with path.open("w") as out:
        out.write("# ff14SB partial charges and AMBER Poisson-Boltzmann radii\n")
        out.write("# Source: AmberTools tleap, leaprc.protein.ff14SB\n")
        out.write("# Template source: amino12.lib, aminont12.lib, aminoct12.lib\n")
        out.write("# PBRadii: mbondi2 (set default PBRadii mbondi2)\n")
        out.write("# CHARGE: AMBER prmtop CHARGE converted to elementary charge by /18.2223\n")
        out.write("# Format: TERMINAL_STATE RESNAME ATOMNAME CHARGE(e) PB_RADIUS(Angstrom)\n")
        out.write("# Terminal states: INTERNAL, NTERM, CTERM\n")
        out.write("# Internal variants present: " + " ".join(supported_internal_variants) + "\n")
        out.write("# Internal variants absent: " + " ".join(unsupported_internal_variants) + "\n")
        out.write("# Terminal variants present: " + " ".join(supported_terminal_variants) + "\n")
        out.write("# Terminal variants absent: " + " ".join(unsupported_terminal_variants) + "\n")
        out.write("# Absent variants must fail explicitly or be prepared by a separate AMBER template/capping workflow.\n")
        for row in rows:
            out.write(
                f"{row.terminal_state} {row.residue_name} {row.atom_name} "
                f"{row.partial_charge:.6f} {row.pb_radius:.4f}\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tleap",
        type=pathlib.Path,
        default=pathlib.Path("/home/jessica/micromamba/envs/mm/bin/tleap"),
        help="AmberTools tleap executable",
    )
    parser.add_argument(
        "--amberhome",
        type=pathlib.Path,
        default=None,
        help="AmberTools root; defaults to AMBERHOME or tleap/..",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("data/ff14sb_params.dat"),
        help="output table path",
    )
    args = parser.parse_args()
    amberhome = infer_amberhome(args.tleap, args.amberhome)

    rows: list[TableRow] = []
    with tempfile.TemporaryDirectory(prefix="nmr_ff14sb_pb_") as tmp:
        workdir = pathlib.Path(tmp)
        for spec in build_specs(amberhome):
            rows.extend(rows_for_template(args.tleap, spec, workdir))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_table(args.output, rows, amberhome)
    print(f"wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
