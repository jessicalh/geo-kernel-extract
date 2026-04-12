#!/usr/bin/env python3
"""Expected ensemble variance from first principles.

For each (atom, ring) pair, the kernel sensitivity to displacement is:
    δ(kernel) / kernel ≈ n × δr / r

where n is the power-law exponent (3 for BS/HM/MC, 5 for PQ, 6 for Disp).
Given a displacement amplitude δr (from MD RMSF), the fractional
variance is (n × δr / r)².

This tells us, for each kernel at each atom:
  - What is the current magnitude (the "mean" we'd accumulate)?
  - What is the expected variance relative to that mean?
  - At what displacement does the variance exceed the mean (signal lost)?
  - Which atoms are in the high-sensitivity regime?

No magic distances — the kernel magnitudes and distances sort themselves.
No model — just dimensional analysis + the measured power laws.

Usage:
    cd learn/src
    python3 -c "
    from mutation_set.config import load_config
    import sys; cfg = load_config('calibration.toml')
    sys.path.insert(0, str(cfg.paths.sdk))
    from actual_physics.expected_ensemble_variance import run
    run(cfg)
    "
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from secondary.loader import iter_proteins, ELEMENT_NAMES, RING_TYPE_NAMES
from mutation_set.config import Config
from mutation_set.kernels import KernelLayout


# Power-law exponents from tensor_realities.py (measured, not assumed):
# BS: -3.04, HM: -3.06, PQ: -5.05
# MC: -3 (dipolar), Coulomb EFG: -3 (1/r³ for point charges)
# Disp: -6 (1/r⁶ for scalar, 1/r⁸ for tensor, but we use scalar × T2)
# RingSusc: -3 (point magnetic dipole, same as MC)
# HBond: -3 (dipolar)
KERNEL_EXPONENTS = {
    "BS": 3.04,
    "HM": 3.06,
    "PQ": 5.05,
    "MC": 3.0,
    "Coulomb": 3.0,
    "EFG": 3.0,
    "Disp": 6.0,
    "RingSusc": 3.0,
    "HBond": 3.0,
    "Chi": 3.0,  # ring susceptibility, same dipolar
}

def _exponent_for_kernel(name):
    """Map kernel name to its distance power-law exponent."""
    for prefix, exp in KERNEL_EXPONENTS.items():
        if name.startswith(prefix):
            return exp
    return 3.0  # default


# Typical RMSF values for different structural contexts (Angstroms)
# From well-tempered metadynamics literature + general MD experience
RMSF_SCENARIOS = {
    "rigid_helix":     0.5,   # backbone heavy atoms in alpha helix
    "stable_sheet":    0.7,   # backbone in beta sheet
    "moderate_loop":   1.5,   # backbone in structured loop
    "flexible_loop":   2.5,   # backbone in disordered loop
    "sidechain_core":  1.0,   # buried sidechain
    "sidechain_surface": 2.0, # exposed sidechain
    "ring_chi1_only":  1.5,   # aromatic ring, only chi1 rotates
    "ring_chi1chi2":   2.5,   # aromatic ring, chi1+chi2 both rotate
}


def run(cfg: Config, max_proteins: int = 0):
    out = Path("output/actual_physics/ensemble_variance")
    out.mkdir(parents=True, exist_ok=True)

    layout = KernelLayout.from_config(cfg)

    ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O"}

    # Per-ring calculator names and their exponents
    per_ring_calcs = {
        "bs": ("Biot-Savart", 3.04),
        "hm": ("Haigh-Mallion", 3.06),
        "pq": ("Pi-Quadrupole", 5.05),
        "chi": ("Ring Susceptibility", 3.0),
    }

    # Accumulators: per element, per calculator
    # For each (atom, ring) pair: distance, kernel_mag, theta
    pair_data = {e: {calc: {"dist": [], "mag": [], "theta": [], "ring_type": []}
                     for calc in per_ring_calcs}
                 for e in ELEMENTS}

    # Also accumulate per-atom total kernel magnitudes for the summed kernels
    atom_data = {e: {"mc_total_mag": [], "coulomb_mag": [], "efg_aro_mag": [],
                      "nearest_dist": [], "n_rings_in_range": []}
                 for e in ELEMENTS}

    print("Loading proteins with per-ring data...")
    n_loaded = 0
    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_loaded += 1
        p = rec.protein
        idx = rec.matched_idx
        rc = p.ring_contributions
        M = len(idx)

        if rc is None or rc.n_pairs == 0:
            continue

        for j in range(M):
            e = int(rec.element[j])
            if e not in ELEMENTS:
                continue

            atom_rc = rc.for_atom(idx[j])
            if atom_rc.n_pairs == 0:
                continue

            # Per-ring pair data
            for ci, (calc, (calc_name, _)) in enumerate(per_ring_calcs.items()):
                t2 = getattr(atom_rc, calc).T2  # (n_pairs, 5)
                mags = np.linalg.norm(t2, axis=1)
                pair_data[e][calc]["dist"].extend(atom_rc.distance.tolist())
                pair_data[e][calc]["mag"].extend(mags.tolist())
                pair_data[e][calc]["theta"].extend(atom_rc.theta.tolist())
                pair_data[e][calc]["ring_type"].extend(atom_rc.ring_type.astype(int).tolist())

            # Per-atom totals
            nearest = np.argmin(atom_rc.distance)
            atom_data[e]["nearest_dist"].append(float(atom_rc.distance[nearest]))
            atom_data[e]["n_rings_in_range"].append(atom_rc.n_pairs)

            # McConnell, Coulomb, EFG from the protein-level shielding objects
            mc_t2 = p.mcconnell.shielding.T2[idx[j]]
            atom_data[e]["mc_total_mag"].append(float(np.linalg.norm(mc_t2)))
            c_t2 = p.coulomb.shielding.T2[idx[j]]
            atom_data[e]["coulomb_mag"].append(float(np.linalg.norm(c_t2)))
            efg_t2 = p.coulomb.efg_aromatic.T2[idx[j]]
            atom_data[e]["efg_aro_mag"].append(float(np.linalg.norm(efg_t2)))

        if n_loaded % 50 == 0:
            print(f"  {n_loaded} proteins")

    print(f"Loaded {n_loaded} proteins\n")

    results = {}

    for e, ename in ELEMENTS.items():
        print("=" * 70)
        print(f"  {ename}")
        print("=" * 70)
        er = {}

        # ── Per-ring-calculator sensitivity analysis ──────────────

        for calc, (calc_name, exponent) in per_ring_calcs.items():
            d = pair_data[e][calc]
            dist = np.array(d["dist"])
            mag = np.array(d["mag"])
            theta = np.array(d["theta"])

            valid = (dist > 1.0) & (mag > 1e-15)
            if valid.sum() < 100:
                continue

            dist_v = dist[valid]
            mag_v = mag[valid]

            n_pairs = valid.sum()

            print(f"\n  {calc_name} ({n_pairs:,d} atom-ring pairs)")
            print(f"  Power law exponent: n = {exponent:.2f}")

            # Fractional sensitivity: δ(kernel)/kernel = n × δr / r
            # For each RMSF scenario, compute the fractional variance
            print(f"\n    Expected fractional std (|δK|/K = n × δr / r):")
            print(f"    {'scenario':<22s}  {'δr':>5s}  "
                  f"{'p10':>8s}  {'p25':>8s}  {'median':>8s}  {'p75':>8s}  {'p90':>8s}")

            scenario_data = {}
            for scenario, rmsf in RMSF_SCENARIOS.items():
                frac_std = exponent * rmsf / dist_v
                pcts = np.percentile(frac_std, [10, 25, 50, 75, 90])
                print(f"    {scenario:<22s}  {rmsf:5.1f}  "
                      f"{pcts[0]:8.3f}  {pcts[1]:8.3f}  {pcts[2]:8.3f}  "
                      f"{pcts[3]:8.3f}  {pcts[4]:8.3f}")
                scenario_data[scenario] = {
                    "rmsf": rmsf,
                    "percentiles": {
                        "p10": round(pcts[0], 4),
                        "p25": round(pcts[1], 4),
                        "median": round(pcts[2], 4),
                        "p75": round(pcts[3], 4),
                        "p90": round(pcts[4], 4),
                    }
                }

            # What fraction of pairs have fractional std > 1
            # (i.e., ensemble variance exceeds mean — signal scrambled)
            print(f"\n    Fraction of pairs where |δK|/K > 1 (signal scrambled):")
            for scenario, rmsf in [("rigid_helix", 0.5),
                                     ("moderate_loop", 1.5),
                                     ("sidechain_surface", 2.0)]:
                frac_std = exponent * rmsf / dist_v
                scrambled = (frac_std > 1.0).mean()
                print(f"      {scenario}: {scrambled:.1%}")
                scenario_data[scenario]["frac_scrambled"] = round(scrambled, 4)

            # Distance distribution of the pairs
            pcts_d = np.percentile(dist_v, [10, 25, 50, 75, 90])
            print(f"\n    Distance distribution (Angstroms):")
            print(f"      p10={pcts_d[0]:.1f}  p25={pcts_d[1]:.1f}  "
                  f"median={pcts_d[2]:.1f}  p75={pcts_d[3]:.1f}  p90={pcts_d[4]:.1f}")

            # Magnitude distribution
            pcts_m = np.percentile(mag_v, [10, 25, 50, 75, 90])
            print(f"    Magnitude distribution:")
            print(f"      p10={pcts_m[0]:.6f}  p25={pcts_m[1]:.6f}  "
                  f"median={pcts_m[2]:.6f}  p75={pcts_m[3]:.6f}  p90={pcts_m[4]:.6f}")

            # The critical distance: where fractional std = 1 for moderate dynamics
            # n × δr / r_crit = 1 → r_crit = n × δr
            for rmsf_val, label in [(0.5, "rigid"), (1.5, "moderate"), (2.5, "flexible")]:
                r_crit = exponent * rmsf_val
                frac_below = (dist_v < r_crit).mean()
                print(f"    Critical distance (δK/K=1 at {label} δr={rmsf_val}A): "
                      f"r_crit = {r_crit:.1f}A, {frac_below:.1%} of pairs below")

            er[calc] = {
                "n_pairs": int(n_pairs),
                "exponent": exponent,
                "scenarios": scenario_data,
                "dist_percentiles": {f"p{p}": round(float(v), 2)
                                     for p, v in zip([10,25,50,75,90], pcts_d)},
                "mag_percentiles": {f"p{p}": round(float(v), 6)
                                    for p, v in zip([10,25,50,75,90], pcts_m)},
            }

        # ── Summed kernel sensitivity ─────────────────────────────

        ad = atom_data[e]
        if len(ad["nearest_dist"]) < 100:
            continue

        nearest_dist = np.array(ad["nearest_dist"])
        mc_mag = np.array(ad["mc_total_mag"])
        coulomb_mag = np.array(ad["coulomb_mag"])
        efg_mag = np.array(ad["efg_aro_mag"])
        n_rings = np.array(ad["n_rings_in_range"])

        N = len(nearest_dist)

        print(f"\n  Summed kernel magnitudes ({N:,d} atoms):")
        for name, mag, exp in [("MC_total", mc_mag, 3.0),
                                 ("Coulomb_total", coulomb_mag, 3.0),
                                 ("EFG_aro", efg_mag, 3.0)]:
            valid = mag > 1e-15
            if valid.sum() < 50:
                continue
            pcts = np.percentile(mag[valid], [10, 25, 50, 75, 90])
            # For summed kernels, sensitivity is more complex because
            # multiple sources contribute.  But dominant source sets the scale.
            # Use nearest ring distance as proxy for dominant source distance.
            frac_std_mod = exp * 1.5 / nearest_dist[valid]
            med_frac = np.median(frac_std_mod)
            print(f"    {name:15s}: median_mag={pcts[2]:.4f}  "
                  f"δK/K at δr=1.5A ≈ {med_frac:.2f}")

        # ── Summary: what the ensemble observer should expect ─────

        print(f"\n  EXPECTED ENSEMBLE REGIME:")
        # For the key ring current kernels:
        bs_pairs = pair_data[e]["bs"]
        if len(bs_pairs["dist"]) > 100:
            bs_dist = np.array(bs_pairs["dist"])
            bs_mag = np.array(bs_pairs["mag"])
            valid = (bs_dist > 1) & (bs_mag > 1e-15)

            # At moderate dynamics (1.5A RMSF):
            frac_std = 3.04 * 1.5 / bs_dist[valid]

            # Categorize atoms by regime
            stable = frac_std < 0.3   # variance < 30% of mean
            dynamic = (frac_std >= 0.3) & (frac_std < 1.0)
            scrambled = frac_std >= 1.0

            print(f"    BS ring current at δr=1.5A RMSF:")
            print(f"      Stable  (δK/K < 0.3):  {stable.mean():6.1%}  "
                  f"(mean dist {bs_dist[valid][stable].mean():.1f}A)" if stable.any() else
                  f"      Stable: 0%")
            print(f"      Dynamic (0.3-1.0):      {dynamic.mean():6.1%}  "
                  f"(mean dist {bs_dist[valid][dynamic].mean():.1f}A)" if dynamic.any() else
                  f"      Dynamic: 0%")
            print(f"      Scrambled (>1.0):        {scrambled.mean():6.1%}  "
                  f"(mean dist {bs_dist[valid][scrambled].mean():.1f}A)" if scrambled.any() else
                  f"      Scrambled: 0%")

        results[ename] = er

    with open(out / "expected_ensemble_variance.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults in {out}/")
