#!/usr/bin/env python3
"""Axis-1 output regression guard.

Usage:
    tools/axis1_keep_regression_guard.py OUT_DIR

The guard is intentionally output-facing: it checks the sidecars and manifest
that the analysis/poster/thesis consume, not just manifest declarations.
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import sys


def fail(msg: str) -> None:
    raise SystemExit(f"axis1_keep_regression_guard: {msg}")


def rows(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        fail(f"missing {path.name}")
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def finite(value: str | None) -> bool:
    if value in (None, ""):
        return False
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def require(condition: bool, msg: str) -> None:
    if not condition:
        fail(msg)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        fail("usage: axis1_keep_regression_guard.py OUT_DIR")
    out = pathlib.Path(argv[1])
    manifest_path = out / "manifest.json"
    if not manifest_path.exists():
        fail("missing manifest.json")
    manifest = json.loads(manifest_path.read_text())

    bounded = rows(out / "bounded_sigma.csv")
    classical = rows(out / "classical_source_terms.csv")
    leaf = rows(out / "classical_source_terms_atom_leaf.csv")
    matrices = rows(out / "source_family_matrices.csv")
    subspace = rows(out / "subspace_overlaps.csv")
    eta = rows(out / "eta2_by_well.csv")

    counts = manifest.get("run", {}).get("object_counts", {})
    require(counts.get("atom") == 846, "atom count changed from 846")
    require(counts.get("ring") == 16, "ring count changed from 16")
    require(counts.get("residue") == 54, "residue count changed from 54")
    require(manifest.get("run", {}).get("sigma_step_count") == 751,
            "sigma step count changed from 751")

    max_roundtrip = max(
        (float(r["sigma_roundtrip_abs"]) for r in bounded if finite(r.get("sigma_roundtrip_abs"))),
        default=math.nan,
    )
    require(math.isfinite(max_roundtrip) and max_roundtrip <= 1.0e-12,
            f"bounded sigma roundtrip too large: {max_roundtrip}")

    for col in ("support_class", "finite_frac", "singleton_flag", "missing_n", "missing_reason"):
        require(col in bounded[0], f"bounded_sigma missing support column {col}")
        require(col in matrices[0], f"source_family_matrices missing support column {col}")
    for col in ("support_class", "finite_frac", "singleton_flag", "missing_n",
                "missing_reason_count", "independence_verdict"):
        require(col in subspace[0], f"subspace_overlaps missing support column {col}")

    require({r.get("granularity") for r in classical} == {"residue_iupac"},
            "classical primary is not residue_iupac")
    require({r.get("granularity") for r in leaf} == {"atom_leaf"},
            "classical leaf audit is not atom_leaf")
    for suffix in ("p05", "p25", "median", "p75", "p95"):
        require(f"contrib_larsen_{suffix}" in classical[0],
                f"classical descriptor missing contrib_larsen_{suffix}")
    require(any(
        finite(r.get("scale_factor")) and finite(r.get("slope_cl_qm"))
        and abs(float(r.get("r_cl_qm") or "nan")) < 0.95
        and abs(float(r["scale_factor"]) - float(r["slope_cl_qm"])) > 1.0e-9
        for r in classical
    ), "no classical row demonstrates scale_factor distinct from slope")
    require(any(r.get("larsen_term_present") == "true" for r in classical),
            "no classical row has Larsen contribution present")

    channel_payload = "\n".join(r.get("channel_names", "") for r in matrices)
    require("shielding_mopac_coulomb" not in channel_payload,
            "forbidden shielding_mopac_coulomb output family remains")
    require("hbond.nearest_dir.x" not in channel_payload
            and "hbond.nearest_dir.y" not in channel_payload
            and "hbond.nearest_dir.z" not in channel_payload,
            "raw hbond nearest_dir xyz channels remain unlabeled")
    require("hbond.nearest_dir.mol_" in channel_payload,
            "molecular-frame hbond nearest_dir channels missing")

    read_ids = {r.get("read_id") for r in subspace}
    for required in ("G2_field_sources", "G6_field_vs_efg", "G4_bs_vs_hm_tensor_components"):
        require(required in read_ids, f"missing subspace read {required}")
    require(any(r.get("read_id") == "G4_bs_vs_hm_tensor_components"
                and (r.get("independence_verdict") or "").startswith("independent_forms_checked")
                for r in subspace),
            "G4 BS/HM independence verdict missing")

    eta_families = {r.get("well_family") for r in eta}
    require("dihedral_rotamer" in eta_families, "eta2 dihedral well missing")
    require(len(eta_families - {"dihedral_rotamer"}) >= 3,
            "eta2 registry is still effectively DIHEDRAL-only")

    lit = manifest.get("literature_constants", {})
    require(lit.get("header") == "src/rediscover/LiteratureConstants.h",
            "LiteratureConstants ledger missing from manifest")

    print(json.dumps({
        "out": str(out),
        "bounded_rows": len(bounded),
        "classical_rows": len(classical),
        "classical_leaf_rows": len(leaf),
        "matrix_rows": len(matrices),
        "subspace_rows": len(subspace),
        "eta_rows": len(eta),
        "max_roundtrip": max_roundtrip,
        "eta_families": sorted(eta_families),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
