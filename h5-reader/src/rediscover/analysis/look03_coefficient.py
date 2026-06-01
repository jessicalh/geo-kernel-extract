#!/usr/bin/env python3
"""look03 — does the classical ring-current FORM drop out, with a sane coefficient?

Reads the EMITTED aggregated reducer outputs (no Python re-sum of the kernel):
  K_valid = sum_dipolar_producer_valid     C++ through-space sum, self/bonded
                                           excluded by identity (= Sigma dipolar)
The per-ring-type emitted sums (sum_dipolar_ringtype_*) carry the geometry split.

Per feedback_no_python_physics_except_labeled_integrity_test (2026-06-01) and the
lead's decision, Python does NOT re-sum `ring_intensity * dipolar` here: the
intensity-weighted aggregate is a C++ reducer output to ADD when wanted (an
intensity-weighted per-type aggregate), not a Python recompute. Until that reducer
lands, the within-atom analysis below uses the emitted unweighted through-space
sum, which is what the producer's own kernel excludes self/bonded from.

Then, within-atom (the modulation the kernel claims):
  - how well does the emitted through-space sum track the producer's BS kernel B
    (bare_T0)?
  - how well does it track DFT sigma_iso D?
  - recover the slopes (the embedded / empirical ring-current constants).
"""
import sys
import numpy as np
import pandas as pd

AGG = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rediscover-out/ring_current_aggregated.csv"

a = pd.read_csv(AGG)
a = a[(a.dft_present == 1) & (a.bare_kernel_present == 1)].copy()
# Emitted reducer output — the through-space sum (self/bonded excluded by identity).
a["K_valid"] = a["sum_dipolar_producer_valid"]
print(f"merged rows: {len(a)}   atoms: {a.atom_index.nunique()}")
print(f"K_valid (emitted sum_dipolar_producer_valid): "
      f"min={a.K_valid.min():+.4g} median={a.K_valid.median():+.4g} max={a.K_valid.max():+.4g}")


def demean(frame, col):
    return frame[col] - frame.groupby("atom_index")[col].transform("mean")


def wcorr_slope(frame, xcol, ycol):
    x = demean(frame, xcol).to_numpy()
    y = demean(frame, ycol).to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or x[ok].std() == 0:
        return np.nan, np.nan
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    slope = np.polyfit(x[ok], y[ok], 1)[0]
    return r, slope


print("\n== within-atom: emitted through-space sum vs the producer kernel B (bare_T0) ==")
r, s = wcorr_slope(a, "K_valid", "bare_T0")
print(f"  B ~ K_valid  r={r:+.4f}  slope={s:+.5g}   (producer's embedded constant)")

print("\n== within-atom: emitted through-space sum vs DFT sigma_iso D ==")
r, s = wcorr_slope(a, "K_valid", "dft_sigma_iso")
print(f"  D ~ K_valid  r={r:+.4f}  slope={s:+.5g} ppm/(kernel unit)")

# Cross-check: producer kernel B vs DFT D directly (their physics, ceiling for ours).
r, s = wcorr_slope(a, "bare_T0", "dft_sigma_iso")
print(f"\n  ceiling  D ~ B(producer)  r={r:+.4f}  slope={s:+.5g} ppm/ppm")

# ---- per-residue: where the signal lives ----
print("\n== per-residue through-space coupling (atoms with real ring motion) ==")
rows = []
for aidx, g in a.groupby("atom_index"):
    if g.K_valid.std() == 0 or len(g) < 10:
        continue
    x = (g.K_valid - g.K_valid.mean()).to_numpy()
    y = (g.dft_sigma_iso - g.dft_sigma_iso.mean()).to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10 or x[ok].std() == 0:
        continue
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    rows.append((aidx, g.atom_name.iloc[0], int(g.residue_number.iloc[0]),
                 g.K_valid.std(), r, g.n_sources_valid.mean()))
pr = pd.DataFrame(rows, columns=["atom_index", "atom_name", "resnum", "sd_Kvalid",
                                 "r_KvalidD", "mean_n_valid"])
pr = pr.sort_values("sd_Kvalid", ascending=False)
with pd.option_context("display.float_format", lambda v: f"{v:+.3g}"):
    print(pr.head(18).to_string(index=False))
print(f"\n  atoms with through-space rings: {len(pr)};  "
      f"median |r_KvalidD|={pr.r_KvalidD.abs().median():.3f};  "
      f"frac |r|>0.5: {(pr.r_KvalidD.abs()>0.5).mean():.2f}")
print("\n  NOTE: intensity-weighted aggregate intentionally omitted (it is a C++ "
      "reducer to add, not a Python re-sum of ring_intensity*dipolar).")
