#!/usr/bin/env python3
"""
T2 analysis: what can we see before learning any parameters?

The calculators store raw geometric kernels in native units:
  BS, HM:      dimensionless (unit current × PPM_FACTOR)
  McConnell:   Å⁻³ (no Δχ applied)
  Coulomb:     V/Å² (no Buckingham γ applied)
  etc.

We CANNOT compute a meaningful residual by summing these — the units
don't match. But we CAN measure:
  1. T2 angular independence between calculators
  2. T2 alignment of each calculator with DFT T2
  3. Scale hierarchy (which kernels are "big" in their own units)
  4. Spatial structure of DFT T2 vs ring proximity

The residual (DFT - parameterised classical) requires regression.
That is a_linear/fit.py.
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from load import (load_conformation, list_proteins, T0, T1, T2)


def t2_norm(t2_vec):
    """L2 norm of T2 vector (5 components per atom). Returns (N,)."""
    return np.sqrt(np.sum(t2_vec ** 2, axis=1))


def cos_similarity_5d(a, b):
    """Mean |cos angle| between two sets of 5D vectors.
    Random vectors in 5D have E[|cos|] ≈ 0.36."""
    norm_a = np.sqrt(np.sum(a**2, axis=1))
    norm_b = np.sqrt(np.sum(b**2, axis=1))
    mask = (norm_a > 1e-12) & (norm_b > 1e-12)
    if mask.sum() == 0:
        return 0.0, 0
    dots = np.sum(a[mask] * b[mask], axis=1) / (norm_a[mask] * norm_b[mask])
    return np.mean(np.abs(dots)), int(mask.sum())


def analyze_protein(protein_id: str):
    """Analyze T2 structure for one protein."""
    wt = load_conformation(protein_id, "wt")
    if wt is None or wt.orca is None:
        return None

    # T2 from each source
    calcs = {
        "BS":           T2(wt.bs),
        "HM":           T2(wt.hm),
        "McConnell":    T2(wt.mc),
        "RingSusc":     T2(wt.ringchi),
        "PiQuad":       T2(wt.pq),
        "Dispersion":   T2(wt.disp),
        "Coulomb":      T2(wt.coulomb),
        "HBond":        T2(wt.hbond),
        "MopacCoulomb": T2(wt.mopac_coulomb),
        "MopacMC":      T2(wt.mopac_mc),
    }
    dft_t2 = T2(wt.orca)

    result = {
        "protein_id": protein_id,
        "n_atoms": wt.n_atoms,
    }

    # Scale: mean |T2| in native units (NOT ppm)
    for name, t2 in calcs.items():
        result[f"scale_{name}"] = t2_norm(t2).mean()
    result["scale_DFT"] = t2_norm(dft_t2).mean()

    # Angular alignment of each calculator T2 with DFT T2
    for name, t2 in calcs.items():
        cos, n = cos_similarity_5d(t2, dft_t2)
        result[f"align_{name}"] = cos

    # T2 independence: pairwise cos between calculators
    names = list(calcs.keys())
    t2s = list(calcs.values())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            cos, n = cos_similarity_5d(t2s[i], t2s[j])
            result[f"indep_{names[i]}_vs_{names[j]}"] = cos

    # DFT T2 by ring proximity (in DFT ppm units — this IS meaningful)
    near = wt.ring_counts[:, 0] > 0   # within 3Å
    mid = (wt.ring_counts[:, 1] > 0) & ~near  # 3-5Å
    far = ~near & ~mid  # >5Å

    for zone, mask in [("near3A", near), ("mid5A", mid), ("far", far)]:
        if mask.any():
            result[f"dft_t2_{zone}"] = t2_norm(dft_t2)[mask].mean()
            result[f"n_{zone}"] = int(mask.sum())
        else:
            result[f"dft_t2_{zone}"] = 0.0
            result[f"n_{zone}"] = 0

    # T0 (isotropic) — for context
    dft_t0 = T0(wt.orca).ravel()
    result["dft_t0_mean"] = np.abs(dft_t0).mean()
    result["dft_t0_std"] = dft_t0.std()

    # T1 (antisymmetric) — this is rarely examined
    dft_t1 = T1(wt.orca)
    t1_norm = np.sqrt(np.sum(dft_t1**2, axis=1))
    result["dft_t1_mean"] = t1_norm.mean()

    # Which calculators produce T1? (magnetic interactions → asymmetric tensors)
    for name, arr in [("BS", wt.bs), ("McConnell", wt.mc),
                       ("RingSusc", wt.ringchi), ("HBond", wt.hbond),
                       ("MopacMC", wt.mopac_mc)]:
        t1 = T1(arr)
        t1n = np.sqrt(np.sum(t1**2, axis=1))
        result[f"t1_{name}"] = t1n.mean()

    return result


def main():
    proteins = list_proteins()
    if not proteins:
        print("No complete pairs found yet.")
        return

    # Also accept proteins with just WT (no pair requirement for this analysis)
    features_dir = Path(__file__).parent / "features" / "FirstExtraction"
    wt_only = []
    for d in sorted(features_dir.iterdir()):
        if d.is_dir() and (d / "wt" / "pos.npy").exists():
            wt_only.append(d.name)

    print(f"Found {len(wt_only)} WT conformations ({len(proteins)} complete pairs)")
    print(f"Analyzing T2 structure...\n")

    results = []
    skipped = 0
    for pid in wt_only:
        r = analyze_protein(pid)
        if r is not None:
            results.append(r)
        else:
            skipped += 1

    print(f"Analyzed: {len(results)}, Skipped (no DFT): {skipped}\n")

    if not results:
        return

    # ========================================================
    # 1. SCALE HIERARCHY
    # ========================================================
    print("=" * 70)
    print("1. SCALE HIERARCHY (mean |T2| in native units)")
    print("=" * 70)
    print("   These are NOT comparable across calculators (different units).")
    print("   They show the raw magnitude of each geometric kernel.\n")

    calc_names = ["BS", "HM", "McConnell", "RingSusc", "PiQuad",
                  "Dispersion", "Coulomb", "HBond",
                  "MopacCoulomb", "MopacMC", "DFT"]
    units = {"BS": "dimless", "HM": "dimless", "McConnell": "Å⁻³",
             "RingSusc": "Å⁻³", "PiQuad": "Å⁻⁵", "Dispersion": "Å⁻⁶",
             "Coulomb": "V/Å²", "HBond": "Å⁻³",
             "MopacCoulomb": "V/Å²", "MopacMC": "Å⁻³", "DFT": "ppm"}
    for name in calc_names:
        vals = [r[f"scale_{name}"] for r in results]
        print(f"  {name:14s}: {np.mean(vals):10.6f}  ±{np.std(vals):.6f}  ({units[name]})")
    print()

    # ========================================================
    # 2. ANGULAR ALIGNMENT WITH DFT T2
    # ========================================================
    print("=" * 70)
    print("2. ANGULAR ALIGNMENT WITH DFT T2")
    print("   |cos| in 5D: random = 0.36, parallel = 1.0")
    print("   This IS unit-independent (direction, not magnitude).")
    print("=" * 70 + "\n")

    for name in ["BS", "HM", "McConnell", "RingSusc", "PiQuad",
                 "Dispersion", "Coulomb", "HBond",
                 "MopacCoulomb", "MopacMC"]:
        vals = [r[f"align_{name}"] for r in results]
        bar = "#" * int(np.mean(vals) * 50)
        print(f"  {name:14s}: {np.mean(vals):.4f}  ±{np.std(vals):.4f}  {bar}")
    print()

    # ========================================================
    # 3. T2 INDEPENDENCE
    # ========================================================
    print("=" * 70)
    print("3. T2 INDEPENDENCE (pairwise |cos|, random = 0.36)")
    print("=" * 70 + "\n")

    pairs = []
    for r in results:
        for key in r:
            if key.startswith("indep_") and key not in [p[0] for p in pairs]:
                pairs.append((key, []))
    for key, vals in pairs:
        for r in results:
            if key in r:
                vals.append(r[key])

    for key, vals in sorted(pairs, key=lambda x: np.mean(x[1]), reverse=True):
        pair_name = key.replace("indep_", "")
        mean_cos = np.mean(vals)
        marker = "***" if mean_cos > 0.7 else "** " if mean_cos > 0.5 else "   "
        print(f"  {marker} {pair_name:35s}: {mean_cos:.4f}")
    print()

    # ========================================================
    # 4. DFT T2 BY RING PROXIMITY
    # ========================================================
    print("=" * 70)
    print("4. DFT T2 MAGNITUDE BY RING PROXIMITY (ppm)")
    print("   Does shielding anisotropy increase near rings?")
    print("=" * 70 + "\n")

    for zone in ["near3A", "mid5A", "far"]:
        vals = [r[f"dft_t2_{zone}"] for r in results if r[f"n_{zone}"] > 0]
        counts = [r[f"n_{zone}"] for r in results if r[f"n_{zone}"] > 0]
        if vals:
            print(f"  {zone:10s}: {np.mean(vals):8.3f} ppm  "
                  f"(mean {np.mean(counts):.0f} atoms/protein)")
    print()

    # ========================================================
    # 5. T1 (ANTISYMMETRIC) — RARELY EXAMINED
    # ========================================================
    print("=" * 70)
    print("5. T1 (ANTISYMMETRIC TENSOR COMPONENT)")
    print("   Physically real for magnetic interactions.")
    print("=" * 70 + "\n")

    dft_t1_vals = [r["dft_t1_mean"] for r in results]
    print(f"  DFT |T1|:       {np.mean(dft_t1_vals):.4f} ppm")
    for name in ["BS", "McConnell", "RingSusc", "HBond", "MopacMC"]:
        vals = [r[f"t1_{name}"] for r in results]
        print(f"  {name:14s} |T1|: {np.mean(vals):.6f}")
    print()

    # ========================================================
    # 6. IRREP BUDGET: T0 vs T1 vs T2
    # ========================================================
    print("=" * 70)
    print("6. DFT IRREP BUDGET (ppm)")
    print("=" * 70 + "\n")

    t0_vals = [r["dft_t0_mean"] for r in results]
    t1_vals = [r["dft_t1_mean"] for r in results]
    t2_vals = [r["scale_DFT"] for r in results]
    print(f"  |T0| (isotropic):     {np.mean(t0_vals):8.3f} ppm")
    print(f"  |T1| (antisymmetric): {np.mean(t1_vals):8.3f} ppm")
    print(f"  |T2| (anisotropic):   {np.mean(t2_vals):8.3f} ppm")
    print(f"  T2/T0 ratio:          {np.mean(t2_vals)/np.mean(t0_vals):.3f}")
    print()


if __name__ == "__main__":
    main()
