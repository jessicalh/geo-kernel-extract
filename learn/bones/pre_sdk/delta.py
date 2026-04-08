#!/usr/bin/env python3
"""
WT-ALA delta T2 analysis.

Match atoms between WT and ALA by position, compute delta shielding.
The delta isolates the effect of the removed ring — backbone contributions
cancel because backbone geometry is nearly identical.

Atom matching: heavy atoms matched by nearest position (< 0.5 Å tolerance).
Hydrogens are excluded (names differ between WT and ALA builds).
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from load import (load_conformation, load_pair, list_proteins,
                  T0, T1, T2, RING_TYPES, _load_optional)


def match_atoms(wt_pos, wt_elem, ala_pos, ala_elem, tol=0.5):
    """Match heavy atoms between WT and ALA by position.

    Returns:
        wt_idx: indices into WT arrays
        ala_idx: corresponding indices into ALA arrays
        distances: match distances in Angstroms
    """
    # Heavy atoms only (element > 1)
    wt_heavy = np.where(wt_elem > 1)[0]
    ala_heavy = np.where(ala_elem > 1)[0]

    wt_idx = []
    ala_idx = []
    distances = []

    # For each ALA heavy atom, find nearest WT heavy atom
    ala_matched = set()
    wt_matched = set()

    for wi in wt_heavy:
        dists = np.sqrt(np.sum((ala_pos[ala_heavy] - wt_pos[wi]) ** 2, axis=1))
        best = np.argmin(dists)
        ai = ala_heavy[best]

        if dists[best] < tol and ai not in ala_matched and wi not in wt_matched:
            # Check element match
            if wt_elem[wi] == ala_elem[ai]:
                wt_idx.append(wi)
                ala_idx.append(ai)
                distances.append(dists[best])
                wt_matched.add(wi)
                ala_matched.add(ai)

    return np.array(wt_idx), np.array(ala_idx), np.array(distances)


def analyze_delta(protein_id: str):
    """Compute and analyze WT-ALA delta T2 for one protein."""
    wt, ala = load_pair(protein_id)
    if wt is None or ala is None:
        return None
    if wt.orca is None or ala.orca is None:
        return None

    wt_idx, ala_idx, match_dist = match_atoms(
        wt.pos, wt.element, ala.pos, ala.element
    )

    if len(wt_idx) < 10:
        return None

    # Delta = WT - ALA (effect of the ring being present)
    delta_dft = wt.orca[wt_idx] - ala.orca[ala_idx]       # (M, 9) ST

    # Diamagnetic and paramagnetic decomposition (if available)
    wt_dia = _load_optional(wt.path / "orca_diamagnetic.npy")
    ala_dia = _load_optional(ala.path / "orca_diamagnetic.npy")
    wt_para = _load_optional(wt.path / "orca_paramagnetic.npy")
    ala_para = _load_optional(ala.path / "orca_paramagnetic.npy")
    has_decomp = all(x is not None for x in [wt_dia, ala_dia, wt_para, ala_para])
    if has_decomp:
        delta_dia = wt_dia[wt_idx] - ala_dia[ala_idx]
        delta_para = wt_para[wt_idx] - ala_para[ala_idx]
    delta_bs = wt.bs[wt_idx] - ala.bs[ala_idx]
    delta_hm = wt.hm[wt_idx] - ala.hm[ala_idx]
    delta_mc = wt.mc[wt_idx] - ala.mc[ala_idx]
    delta_coulomb = wt.coulomb[wt_idx] - ala.coulomb[ala_idx]
    delta_ringchi = wt.ringchi[wt_idx] - ala.ringchi[ala_idx]
    delta_pq = wt.pq[wt_idx] - ala.pq[ala_idx]
    delta_disp = wt.disp[wt_idx] - ala.disp[ala_idx]
    delta_hbond = wt.hbond[wt_idx] - ala.hbond[ala_idx]

    delta_mopac_coulomb = wt.mopac_coulomb[wt_idx] - ala.mopac_coulomb[ala_idx]
    delta_mopac_mc = wt.mopac_mc[wt_idx] - ala.mopac_mc[ala_idx]

    # T2 of deltas
    dft_t2 = T2(delta_dft)
    calcs = {
        "BS":           T2(delta_bs),
        "HM":           T2(delta_hm),
        "McConnell":    T2(delta_mc),
        "Coulomb":      T2(delta_coulomb),
        "RingSusc":     T2(delta_ringchi),
        "PiQuad":       T2(delta_pq),
        "Dispersion":   T2(delta_disp),
        "HBond":        T2(delta_hbond),
        "MopacCoulomb": T2(delta_mopac_coulomb),
        "MopacMC":      T2(delta_mopac_mc),
    }

    def t2_norm(v):
        return np.sqrt(np.sum(v ** 2, axis=1))

    def cos5d(a, b):
        na = np.sqrt(np.sum(a**2, axis=1))
        nb = np.sqrt(np.sum(b**2, axis=1))
        mask = (na > 1e-12) & (nb > 1e-12)
        if mask.sum() == 0:
            return 0.0, 0.0, 0
        dots = np.sum(a[mask] * b[mask], axis=1) / (na[mask] * nb[mask])
        return np.mean(np.abs(dots)), np.mean(dots), int(mask.sum())

    result = {
        "protein_id": protein_id,
        "n_matched": len(wt_idx),
        "n_wt": wt.n_atoms,
        "n_ala": ala.n_atoms,
        "mean_match_dist": match_dist.mean(),
        "delta_dft_t2_mean": t2_norm(dft_t2).mean(),
        "delta_dft_t0_mean": np.abs(T0(delta_dft)).mean(),
    }

    # Ring proximity of matched atoms (from WT)
    ring_counts = wt.ring_counts[wt_idx]
    near = ring_counts[:, 0] > 0  # within 3Å
    mid = (ring_counts[:, 1] > 0) & ~near  # 3-5Å
    far = ~near & ~mid

    for zone, mask in [("near3A", near), ("mid5A", mid), ("far", far)]:
        if mask.any():
            result[f"delta_dft_t2_{zone}"] = t2_norm(dft_t2)[mask].mean()
            result[f"n_{zone}"] = int(mask.sum())
        else:
            result[f"delta_dft_t2_{zone}"] = 0.0
            result[f"n_{zone}"] = 0

    # Alignment of each calculator's delta T2 with DFT delta T2
    for name, ct2 in calcs.items():
        abs_cos, signed_cos, n = cos5d(ct2, dft_t2)
        result[f"align_{name}"] = abs_cos
        result[f"sign_{name}"] = signed_cos
        result[f"scale_{name}"] = t2_norm(ct2).mean()

    # Diamagnetic vs paramagnetic T2
    if has_decomp:
        result["delta_dia_t2_mean"] = t2_norm(T2(delta_dia)).mean()
        result["delta_para_t2_mean"] = t2_norm(T2(delta_para)).mean()
        result["delta_dia_t0_mean"] = np.abs(T0(delta_dia)).mean()
        result["delta_para_t0_mean"] = np.abs(T0(delta_para)).mean()
        # Which part of DFT delta T2 is diamagnetic vs paramagnetic?
        abs_cos_dia, _, _ = cos5d(T2(delta_dia), dft_t2)
        abs_cos_para, _, _ = cos5d(T2(delta_para), dft_t2)
        result["align_dia_with_total"] = abs_cos_dia
        result["align_para_with_total"] = abs_cos_para
        # BS alignment with each component
        abs_cos_bs_dia, _, _ = cos5d(T2(delta_bs), T2(delta_dia))
        abs_cos_bs_para, _, _ = cos5d(T2(delta_bs), T2(delta_para))
        result["align_BS_with_dia"] = abs_cos_bs_dia
        result["align_BS_with_para"] = abs_cos_bs_para

    # Element breakdown
    for elem_z, elem_name in [(6, "C"), (7, "N"), (8, "O")]:
        emask = wt.element[wt_idx] == elem_z
        if emask.sum() > 5:
            result[f"delta_dft_t2_{elem_name}"] = t2_norm(dft_t2)[emask].mean()
            for name, ct2 in calcs.items():
                ac, sc, _ = cos5d(ct2[emask], dft_t2[emask])
                result[f"align_{name}_{elem_name}"] = ac

    return result


def main():
    proteins = list_proteins()
    if not proteins:
        print("No complete pairs found.")
        return

    print(f"Computing WT-ALA delta T2 for {len(proteins)} proteins...\n")

    results = []
    for pid in proteins:
        r = analyze_delta(pid)
        if r is not None:
            results.append(r)

    print(f"Analyzed: {len(results)} proteins\n")
    if not results:
        return

    # ========================================================
    print("=" * 70)
    print("WT-ALA DELTA T2 ANALYSIS")
    print("  Delta = WT - ALA = effect of the ring being present")
    print("=" * 70 + "\n")

    # Matching stats
    n_matched = [r["n_matched"] for r in results]
    print(f"Atom matching: {np.mean(n_matched):.0f} heavy atoms matched/protein "
          f"(mean dist {np.mean([r['mean_match_dist'] for r in results]):.4f} Å)\n")

    # Delta DFT scale
    dt2 = [r["delta_dft_t2_mean"] for r in results]
    dt0 = [r["delta_dft_t0_mean"] for r in results]
    print(f"Delta DFT magnitude (ppm):")
    print(f"  |ΔT0| mean: {np.mean(dt0):.4f}")
    print(f"  |ΔT2| mean: {np.mean(dt2):.4f}")
    print()

    # Delta T2 by ring proximity
    print("Delta |T2| by ring proximity:")
    for zone in ["near3A", "mid5A", "far"]:
        vals = [r[f"delta_dft_t2_{zone}"] for r in results if r.get(f"n_{zone}", 0) > 0]
        counts = [r[f"n_{zone}"] for r in results if r.get(f"n_{zone}", 0) > 0]
        if vals:
            print(f"  {zone:10s}: {np.mean(vals):.4f} ppm  "
                  f"({np.mean(counts):.0f} atoms/protein)")
    print()

    # THE KEY RESULT: alignment of calculator delta T2 with DFT delta T2
    print("Calculator delta T2 alignment with DFT delta T2:")
    print(f"  {'Calculator':12s} {'|cos|':>8s} {'signed cos':>12s} {'|ΔT2|':>10s}")
    print("  " + "-" * 46)
    calc_names = ["BS", "HM", "McConnell", "RingSusc", "PiQuad",
                  "Dispersion", "Coulomb", "HBond",
                  "MopacCoulomb", "MopacMC"]
    for name in calc_names:
        abs_cos = np.mean([r[f"align_{name}"] for r in results])
        sign_cos = np.mean([r[f"sign_{name}"] for r in results])
        scale = np.mean([r[f"scale_{name}"] for r in results])
        bar = "#" * int(abs_cos * 50)
        print(f"  {name:12s} {abs_cos:8.4f} {sign_cos:+12.4f} {scale:10.6f}  {bar}")
    print(f"  {'random 5D':12s} {'0.3600':>8s}")
    print()

    # Diamagnetic vs paramagnetic
    dia_keys = ["delta_dia_t2_mean", "delta_para_t2_mean"]
    if all(dia_keys[0] in r for r in results):
        print("Diamagnetic vs paramagnetic decomposition of delta T2:")
        dia_t2 = [r["delta_dia_t2_mean"] for r in results]
        para_t2 = [r["delta_para_t2_mean"] for r in results]
        print(f"  |Δ diamagnetic T2|:  {np.mean(dia_t2):.4f} ppm")
        print(f"  |Δ paramagnetic T2|: {np.mean(para_t2):.4f} ppm")
        dia_t0 = [r["delta_dia_t0_mean"] for r in results]
        para_t0 = [r["delta_para_t0_mean"] for r in results]
        print(f"  |Δ diamagnetic T0|:  {np.mean(dia_t0):.4f} ppm")
        print(f"  |Δ paramagnetic T0|: {np.mean(para_t0):.4f} ppm")
        print()
        align_dia = [r["align_dia_with_total"] for r in results]
        align_para = [r["align_para_with_total"] for r in results]
        print(f"  Diamagnetic T2 alignment with total delta T2:  {np.mean(align_dia):.4f}")
        print(f"  Paramagnetic T2 alignment with total delta T2: {np.mean(align_para):.4f}")
        bs_dia = [r["align_BS_with_dia"] for r in results]
        bs_para = [r["align_BS_with_para"] for r in results]
        print(f"  BS alignment with diamagnetic delta T2:  {np.mean(bs_dia):.4f}")
        print(f"  BS alignment with paramagnetic delta T2: {np.mean(bs_para):.4f}")
        print()

    # Per-element breakdown
    for elem_name in ["C", "N", "O"]:
        key = f"delta_dft_t2_{elem_name}"
        vals = [r[key] for r in results if key in r]
        if vals:
            print(f"Delta |T2| for {elem_name}: {np.mean(vals):.4f} ppm")
            for name in ["BS", "HM", "McConnell", "Coulomb"]:
                akey = f"align_{name}_{elem_name}"
                avals = [r[akey] for r in results if akey in r]
                if avals:
                    print(f"  {name:12s} |cos|={np.mean(avals):.4f}")
        print()


if __name__ == "__main__":
    main()
