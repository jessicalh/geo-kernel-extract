#!/usr/bin/env python3
"""Verify all Stage 1 numbers against the JSON outputs.

Checks that the numbers cited in stage1-mutations/notes/ match the
actual analysis output.  Run after run_all.sh to confirm reproducibility.

Usage:
    cd learn/src
    python3 ../stage1-mutations/analysis/verify_numbers.py
"""

import json
import sys
from pathlib import Path

OUT = Path("output/actual_physics")
PASS = 0
FAIL = 0


def check(name: str, actual: float, expected: float, tol: float = 0.002):
    global PASS, FAIL
    if abs(actual - expected) <= tol:
        PASS += 1
    else:
        FAIL += 1
        print(f"  FAIL: {name}: got {actual:.4f}, expected {expected:.4f} "
              f"(diff {actual - expected:+.4f})")


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


print("Verifying Stage 1 numbers against analysis output...\n")

# ── Full space analysis ──────────────────────────────────────────
print("Full space analysis:")
fsa = load_json(OUT / "full_space" / "full_space_analysis.json")

# Per-group R² (normalised) — from notes/normalisation_physics.md
check("H ring_current norm", fsa["H"]["per_group_r2"]["norm"]["ring_current"], 0.782, 0.005)
check("H ALL norm", fsa["H"]["per_group_r2"]["norm"]["ALL"], 0.856, 0.005)
check("O dispersion norm", fsa["O"]["per_group_r2"]["norm"]["dispersion"], 0.234, 0.005)
check("O dispersion raw", fsa["O"]["per_group_r2"]["raw"]["dispersion"], 0.058, 0.005)
check("N ALL norm", fsa["N"]["per_group_r2"]["norm"]["ALL"], 0.267, 0.005)
check("C mopac_efg raw", fsa["C"]["per_group_r2"]["raw"]["mopac_efg"], 0.392, 0.005)
check("C ff14sb_efg raw", fsa["C"]["per_group_r2"]["raw"]["ff14sb_efg"], 0.223, 0.005)

# Atom counts
check("H atoms", fsa["H"]["n_atoms"], 127460, 500)
check("C atoms", fsa["C"]["n_atoms"], 74759, 500)

print(f"  {PASS} passed\n")

# ── Dimensionality ───────────────────────────────────────────────
print("Dimensionality:")
dim = load_json(OUT / "dimensionality" / "dimensionality.json")

check("H norm plateau", dim["H"]["pca_ridge"]["plateau_k"], 20, 3)
check("C norm plateau", dim["C"]["pca_ridge"]["plateau_k"], 6, 2)
check("N norm plateau", dim["N"]["pca_ridge"]["plateau_k"], 3, 1)
check("O norm plateau", dim["O"]["pca_ridge"]["plateau_k"], 12, 3)

# All raw plateaus should be 3
for en in ["H", "C", "N", "O"]:
    check(f"{en} raw plateau", dim[en]["pca_ridge_raw_plateau"], 3, 0)

# Near-field dims
check("N near dims", dim["N"]["distance_bands"]["0-4A"]["predictive_dims"], 4, 1)
check("O near dims", dim["O"]["distance_bands"]["0-4A"]["predictive_dims"], 8, 2)

print(f"  {PASS} passed\n")

# ── Geometry-only basis ──────────────────────────────────────────
print("Geometry-only basis:")
geo = load_json(OUT / "geometry_only_basis" / "geometry_only_basis.json")

check("H fair full", geo["H"]["progressive"]["fair"]["full"], 0.928, 0.005)
check("H fair geo", geo["H"]["progressive"]["fair"]["geo"], 0.921, 0.005)
check("C fair full", geo["C"]["progressive"]["fair"]["full"], 0.562, 0.010)
check("C fair geo-full gap",
      geo["C"]["progressive"]["fair"]["full"] - geo["C"]["progressive"]["fair"]["geo"],
      0.197, 0.010)

print(f"  {PASS} passed\n")

# ── Summary ──────────────────────────────────────────────────────
print(f"{'=' * 50}")
print(f"  TOTAL: {PASS} passed, {FAIL} failed")
if FAIL > 0:
    print(f"  *** {FAIL} numbers do not match cited values ***")
    sys.exit(1)
else:
    print(f"  All numbers verified.")
