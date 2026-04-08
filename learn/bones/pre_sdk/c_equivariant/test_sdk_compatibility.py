#!/usr/bin/env python3
"""
Validate that nmr_extract SDK can load calibration data and that all
fields needed by the calibration pipeline are accessible.

Tests against CalibrationExtractionTest (existing data). Run with:
    python learn/c_equivariant/test_sdk_compatibility.py
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
FEATURES = REPO / "calibration" / "features" / "GatedCalibration"

# SDK import
try:
    from nmr_extract import (
        load, Protein, RingType, BondCategory,
        ShieldingTensor, EFGTensor, VectorField,
        PerRingTypeT2, PerBondCategoryT2,
        RingContributions, RingGeometry,
    )
    print("OK  nmr_extract imported")
except ImportError as e:
    print(f"FAIL  nmr_extract import: {e}")
    print("      pip install -e /shared/2026Thesis/nmr-shielding/python")
    sys.exit(1)

# Pick a protein that has delta (mutant extraction)
test_proteins = sorted(d.name for d in FEATURES.iterdir() if d.is_dir())
if not test_proteins:
    print(f"FAIL  no proteins in {FEATURES}")
    sys.exit(1)

# Find one with delta
pid = None
for candidate in test_proteins[:20]:
    if (FEATURES / candidate / "delta_shielding.npy").exists():
        pid = candidate
        break
if pid is None:
    pid = test_proteins[0]
    print(f"WARN  no mutant protein found in first 20, using {pid}")

print(f"Loading {pid}...")
p = load(str(FEATURES / pid))

errors = 0

def check(label, condition, detail=""):
    global errors
    if condition:
        print(f"OK  {label}")
    else:
        print(f"FAIL  {label}  {detail}")
        errors += 1

# ── Identity ──
check("protein loaded", isinstance(p, Protein))
check("n_atoms > 0", p.n_atoms > 0, f"n_atoms={p.n_atoms}")
check("element shape", p.element.shape == (p.n_atoms,))
N = p.n_atoms

# ── Ring current kernels (needed for kernel layout slots 0-23) ──
check("biot_savart.shielding", isinstance(p.biot_savart.shielding, ShieldingTensor))
check("biot_savart.per_type_T2", isinstance(p.biot_savart.per_type_T2, PerRingTypeT2))
check("BS per_type_T2 shape", p.biot_savart.per_type_T2.data.shape == (N, 40))
check("BS as_block shape", p.biot_savart.per_type_T2.as_block().shape == (N, 8, 5))
check("BS T2 for PHE", p.biot_savart.per_type_T2.for_type(RingType.PHE).shape == (N, 5))
check("BS shielding.T2", p.biot_savart.shielding.T2.shape == (N, 5))
check("BS shielding.T0", p.biot_savart.shielding.T0.shape == (N, 1))
check("BS shielding.isotropic", p.biot_savart.shielding.isotropic.shape == (N,))

check("haigh_mallion.per_type_T2", p.haigh_mallion.per_type_T2.data.shape == (N, 40))
check("dispersion.per_type_T2", p.dispersion.per_type_T2.data.shape == (N, 40))

# ── Per-type T0 (needed for scalar features) ──
check("BS per_type_T0", p.biot_savart.per_type_T0.data.shape == (N, 8))
check("HM per_type_T0", p.haigh_mallion.per_type_T0.data.shape == (N, 8))
check("Disp per_type_T0", p.dispersion.per_type_T0.data.shape == (N, 8))

# ── Bond calculators (kernel layout slots 24-39) ──
check("mcconnell.shielding", p.mcconnell.shielding.T2.shape == (N, 5))
check("mcconnell.category_T2", p.mcconnell.category_T2.data.shape == (N, 25))
check("MC category as_block", p.mcconnell.category_T2.as_block().shape == (N, 5, 5))
check("mcconnell.scalars", p.mcconnell.scalars.data.shape == (N, 6))

check("coulomb.shielding", p.coulomb.shielding.T2.shape == (N, 5))
check("coulomb.E", isinstance(p.coulomb.E, VectorField))
check("coulomb.efg_backbone", isinstance(p.coulomb.efg_backbone, EFGTensor))
check("coulomb.efg_aromatic", isinstance(p.coulomb.efg_aromatic, EFGTensor))
check("coulomb.scalars", p.coulomb.scalars.data.shape == (N, 4))

check("hbond.shielding", p.hbond.shielding.T2.shape == (N, 5))
check("hbond.scalars", p.hbond.scalars.data.shape == (N, 3))

check("ring_susceptibility", isinstance(p.ring_susceptibility, ShieldingTensor))

# ── EFG kernels (layout slots 40-45) ──
check("EFG backbone T2", p.coulomb.efg_backbone.T2.shape == (N, 5))
check("EFG aromatic T2", p.coulomb.efg_aromatic.T2.shape == (N, 5))

# ── MOPAC (optional but needed for MopacEFG_aro — the dominant kernel) ──
if p.mopac is not None:
    check("mopac.core.charges", p.mopac.core.charges.shape == (N,))
    check("mopac.coulomb.shielding", p.mopac.coulomb.shielding.T2.shape == (N, 5))
    check("mopac.coulomb.efg_backbone", isinstance(p.mopac.coulomb.efg_backbone, EFGTensor))
    check("mopac.coulomb.efg_aromatic", isinstance(p.mopac.coulomb.efg_aromatic, EFGTensor))
    check("mopac.mcconnell.shielding", p.mopac.mcconnell.shielding.T2.shape == (N, 5))
    check("mopac.mcconnell.category_T2", p.mopac.mcconnell.category_T2.data.shape == (N, 25))
else:
    print("WARN  mopac is None — MopacEFG_aro (dominant kernel) unavailable")

# ── APBS (optional) ──
if p.apbs is not None:
    check("apbs.efg", isinstance(p.apbs.efg, EFGTensor))

# ── Delta (mutation comparison) ──
if p.delta is not None:
    check("delta.shielding", isinstance(p.delta.shielding, ShieldingTensor))
    check("delta.shielding.T2", p.delta.shielding.T2.shape == (N, 5))
    check("delta.scalars", p.delta.scalars.data.shape == (N, 6))
    check("delta matched_mask", hasattr(p.delta.scalars, 'matched_mask'))
    n_matched = p.delta.scalars.matched_mask.sum()
    check(f"delta matched atoms ({n_matched})", n_matched > 10)

    if p.delta.ring_proximity is not None:
        nr = p.delta.ring_proximity.n_removed_rings
        check(f"ring_proximity n_removed={nr}", nr > 0)
        check("ring_proximity distance(0)", p.delta.ring_proximity.distance(0).shape == (N,))
else:
    print("WARN  delta is None — not a mutant extraction")

# ── Ring contributions (sparse per-ring data) ──
check("ring_contributions", isinstance(p.ring_contributions, RingContributions))
rc = p.ring_contributions
check("rc.n_pairs > 0", rc.n_pairs > 0, f"n_pairs={rc.n_pairs}")
check("rc.bs.T2", rc.bs.T2.shape == (rc.n_pairs, 5))
check("rc.hm.T2", rc.hm.T2.shape == (rc.n_pairs, 5))
check("rc.hm_H.T2 (raw integral)", rc.hm_H.T2.shape == (rc.n_pairs, 5))

# ── Ring geometry ──
check("ring_geometry", isinstance(p.ring_geometry, RingGeometry))
rg = p.ring_geometry
check("rg.n_rings > 0", rg.n_rings > 0)
check("rg.center", rg.center.shape == (rg.n_rings, 3))
check("rg.normal", rg.normal.shape == (rg.n_rings, 3))

# ── Torch conversion ──
t = p.biot_savart.shielding.torch()
check("torch() works", t.shape == (N, 9))
check("torch dtype", str(t.dtype) == "torch.float64")

# ── Irreps ──
check("shielding irreps", str(p.biot_savart.shielding.irreps) == "1x0e+1x1o+1x2e")
check("per_type_T2 irreps", str(p.biot_savart.per_type_T2.irreps) == "8x2e")
check("category_T2 irreps", str(p.mcconnell.category_T2.irreps) == "5x2e")

# ── Cross-check: sparse sums match per-type ──
for rt in RingType:
    subset = rc.for_ring_type(rt)
    if subset.n_pairs == 0:
        continue
    sparse_sum = np.zeros((N, 5))
    np.add.at(sparse_sum, subset.atom_index, subset.bs.T2)
    presummed = p.biot_savart.per_type_T2.for_type(rt)
    if np.allclose(sparse_sum, presummed, atol=1e-12):
        check(f"sparse BS_{rt.name} matches per_type", True)
    else:
        max_err = np.abs(sparse_sum - presummed).max()
        check(f"sparse BS_{rt.name} matches per_type", False, f"max_err={max_err:.2e}")

print(f"\n{'='*50}")
if errors == 0:
    print(f"ALL CHECKS PASSED — SDK loads calibration data correctly")
else:
    print(f"{errors} CHECKS FAILED")
sys.exit(errors)
