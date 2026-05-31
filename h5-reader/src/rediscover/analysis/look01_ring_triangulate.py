#!/usr/bin/env python3
"""look01 — ring current, first honest look at the substrate.

Runs on the CURRENT emitted CSVs (self-ring NOT yet excluded), because the
within-atom de-meaned test is largely robust to a near-constant self-ring
offset: the aromatic H sits at a roughly fixed spot in its own ring, so that
term de-means away. This is a FIRST read to see whether there is any structure
worth the re-run, not the final analysis.

Three scalars per (atom, frame):
  S = sum_dipolar       our recomputed Pople sum  Sigma (3cos^2-1)/r^3
  B = bare_T0           the producer's own BS ring-current kernel (ppm)
  D = dft_sigma_iso     the DFT truth, T0 = trace/3 (ppm), rotation-safe

Three questions, kept apart:
  (1) S vs B  -> do WE reproduce the producer's kernel?      (extraction check)
  (2) B vs D  -> does the producer's kernel track DFT?       (their physics)
  (3) S vs D  -> does our sum track DFT?                      (our physics)
Each asked pooled AND within-atom de-meaned (the ring-current claim is a
modulation on a per-atom chemical baseline; pooling conflates between-atom
chemistry with within-atom geometry).
"""
import sys
import numpy as np
import pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rediscover-out/ring_current_aggregated.csv"

df = pd.read_csv(CSV)
print(f"loaded {CSV}")
print(f"  rows={len(df)}  unique atoms={df.atom_index.nunique()}  unique frames={df.h5_row.nunique()}")

# Keep only rows where both the DFT target and the producer kernel are present.
m = (df.dft_present == 1) & (df.bare_kernel_present == 1)
d = df[m].copy()
print(f"  rows with dft & bare present: {len(d)}  ({d.atom_index.nunique()} atoms)")

S = d["sum_dipolar"].to_numpy()
B = d["bare_T0"].to_numpy()
D = d["dft_sigma_iso"].to_numpy()

def stats(name, x):
    print(f"  {name:14s} n={len(x):6d}  mean={np.nanmean(x):+.4g}  sd={np.nanstd(x):.4g}"
          f"  min={np.nanmin(x):+.4g}  max={np.nanmax(x):+.4g}  nan={np.isnan(x).sum()}")

print("\n== marginal distributions ==")
stats("S sum_dipolar", S)
stats("B bare_T0", B)
stats("D dft_sigma_iso", D)

def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.nanstd(a[ok]) == 0 or np.nanstd(b[ok]) == 0:
        return np.nan, ok.sum()
    return np.corrcoef(a[ok], b[ok])[0, 1], ok.sum()

print("\n== pooled correlations (all atoms together) ==")
for label, a, b in [("S vs B (extraction)", S, B),
                    ("B vs D (their physics)", B, D),
                    ("S vs D (our physics)", S, D)]:
    r, n = corr(a, b)
    print(f"  {label:26s} r={r:+.4f}  (n={n})")

# ---- within-atom de-meaning ----
# Subtract each atom's own across-frame mean from S, B, D, then pool. This is
# the part the ring-current kernel actually claims: geometric modulation.
print("\n== within-atom de-meaned (the ring-current modulation claim) ==")
d["S_w"] = S - d.groupby("atom_index")["sum_dipolar"].transform("mean")
d["B_w"] = B - d.groupby("atom_index")["bare_T0"].transform("mean")
d["D_w"] = D - d.groupby("atom_index")["dft_sigma_iso"].transform("mean")
for label, a, b in [("S vs B", d.S_w.to_numpy(), d.B_w.to_numpy()),
                    ("B vs D", d.B_w.to_numpy(), d.D_w.to_numpy()),
                    ("S vs D", d.S_w.to_numpy(), d.D_w.to_numpy())]:
    r, n = corr(a, b)
    # OLS slope b on a (the coefficient k), via least squares on finite rows
    ok = np.isfinite(a) & np.isfinite(b)
    slope = np.polyfit(a[ok], b[ok], 1)[0] if ok.sum() > 2 and np.std(a[ok]) > 0 else np.nan
    print(f"  {label:8s} r={r:+.4f}  slope(b~a)={slope:+.4g}  (n={n})")

# ---- per-atom coupling distribution ----
# Within each atom, correlate S_t and D_t across frames. A strong coupling for
# atoms with a real moving ring nearby is the signal; flat atoms are the floor.
print("\n== per-atom S-vs-D coupling across frames ==")
rows = []
for aidx, g in d.groupby("atom_index"):
    if len(g) < 10:
        continue
    s, dd = g["sum_dipolar"].to_numpy(), g["dft_sigma_iso"].to_numpy()
    r, n = corr(s, dd)
    rows.append((aidx, g["atom_name"].iloc[0], int(g["residue_number"].iloc[0]),
                 len(g), np.std(s), np.std(dd), r))
pa = pd.DataFrame(rows, columns=["atom_index", "atom_name", "resnum", "n_frames",
                                 "sd_S", "sd_D", "r_SD"])
print(f"  atoms with >=10 frames: {len(pa)}")
print(f"  median |r_SD| = {pa.r_SD.abs().median():.3f}   "
      f"frac |r|>0.3: {(pa.r_SD.abs()>0.3).mean():.2f}   "
      f"frac |r|>0.5: {(pa.r_SD.abs()>0.5).mean():.2f}")
print("\n  most strongly coupled atoms (by |r_SD|, needs ring motion: high sd_S):")
top = pa.reindex(pa.r_SD.abs().sort_values(ascending=False).index).head(15)
with pd.option_context("display.float_format", lambda v: f"{v:+.3g}"):
    print(top.to_string(index=False))

print("\n  atoms with the most ring-current modulation (highest sd_S):")
topvar = pa.sort_values("sd_S", ascending=False).head(15)
with pd.option_context("display.float_format", lambda v: f"{v:+.3g}"):
    print(topvar.to_string(index=False))
