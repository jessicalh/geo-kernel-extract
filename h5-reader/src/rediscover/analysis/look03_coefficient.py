#!/usr/bin/env python3
"""look03 — does the classical ring-current FORM drop out, with a sane coefficient?

Builds, from the per-source table (excluding the self ring by geometry, pending
the identity-based C++ fix), two through-space kernels:
  K_raw = Sigma (3cos^2-1)/r^3                       unweighted geometry
  K_int = Sigma intensity * (3cos^2-1)/r^3           intensity-weighted (Pople form)
intensity = QtRing.LiteratureIntensity (nA/T), emitted per source.

Then, within-atom (the modulation the kernel claims):
  - which kernel best tracks the producer's own BS kernel B (bare_T0)?
  - which best tracks DFT sigma_iso D?
  - recover slope(B ~ K_int): the producer's embedded ring-current constant.
  - recover slope(D ~ K_int): the empirical susceptibility from DFT.
The producer's kernel embeds the textbook intensities; if K_int reproduces B up
to a single scale, we have recovered the classical constant from geometry alone.
"""
import sys
import numpy as np
import pandas as pd

SRC = "/tmp/rediscover-out/ring_current_sources.csv"
AGG = "/tmp/rediscover-out/ring_current_aggregated.csv"
SELF_R = 3.0  # A; self/bonded ring proxy (self median 2.48, 95% 2.52; through-space >3.5)

src = pd.read_csv(SRC).rename(columns={"dipolar_3cos2m1_over_r3": "dipolar"})
key = ["atom_index", "h5_row"]

ts = src[src.r >= SELF_R].copy()  # through-space only
n_self = len(src) - len(ts)
print(f"sources: {len(src)} total, excluded {n_self} self/bonded (r<{SELF_R}A), "
      f"{len(ts)} through-space")
print(f"through-space r: min={ts.r.min():.2f} median={ts.r.median():.2f} A   "
      f"intensity range: {ts.ring_intensity.min():.3g}..{ts.ring_intensity.max():.3g} nA/T")

ts["w"] = ts.ring_intensity * ts.dipolar
k = ts.groupby(key).agg(K_raw=("dipolar", "sum"), K_int=("w", "sum"),
                        n_ts=("dipolar", "size")).reset_index()

a = pd.read_csv(AGG)
a = a[(a.dft_present == 1) & (a.bare_kernel_present == 1)]
a = a[key + ["bare_T0", "dft_sigma_iso", "atom_name", "residue_number", "amino_acid_ord"]]
df = a.merge(k, on=key, how="left").fillna({"K_raw": 0.0, "K_int": 0.0, "n_ts": 0})
print(f"merged rows: {len(df)}   (atoms with >=1 through-space ring: "
      f"{(df.groupby('atom_index').n_ts.max() > 0).sum()}/{df.atom_index.nunique()})")

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

print("\n== within-atom: which kernel best tracks the producer kernel B (bare_T0)? ==")
for kc in ["K_raw", "K_int"]:
    r, s = wcorr_slope(df, kc, "bare_T0")
    print(f"  B ~ {kc:6s}  r={r:+.4f}  slope={s:+.5g}   (producer's embedded constant)")

print("\n== within-atom: which kernel best tracks DFT sigma_iso D? ==")
for kc in ["K_raw", "K_int"]:
    r, s = wcorr_slope(df, kc, "dft_sigma_iso")
    print(f"  D ~ {kc:6s}  r={r:+.4f}  slope={s:+.5g} ppm/(kernel unit)")

# Cross-check: producer kernel B vs DFT D directly (their physics, ceiling for ours).
r, s = wcorr_slope(df, "bare_T0", "dft_sigma_iso")
print(f"\n  ceiling  D ~ B(producer)  r={r:+.4f}  slope={s:+.5g} ppm/ppm")

# ---- per-residue: where the signal lives, and the neighbour ring type ----
print("\n== per-residue through-space coupling (atoms with real ring motion) ==")
rows = []
for aidx, g in df.groupby("atom_index"):
    if g.n_ts.max() == 0 or len(g) < 10:
        continue
    x = (g.K_int - g.K_int.mean()).to_numpy()
    y = (g.dft_sigma_iso - g.dft_sigma_iso.mean()).to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10 or x[ok].std() == 0:
        continue
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    # dominant neighbour ring type for this atom's through-space sources
    sg = ts[ts.atom_index == aidx]
    dom = sg.groupby("ring_type_index").dipolar.apply(lambda v: v.abs().sum())
    dom_type = int(dom.idxmax()) if len(dom) else -1
    rows.append((aidx, g.atom_name.iloc[0], int(g.residue_number.iloc[0]),
                 g.K_int.std(), r, g.n_ts.mean(), dom_type))
pr = pd.DataFrame(rows, columns=["atom_index", "atom_name", "resnum", "sd_Kint",
                                 "r_KintD", "mean_n_ts", "dom_ring_type"])
pr = pr.sort_values("sd_Kint", ascending=False)
with pd.option_context("display.float_format", lambda v: f"{v:+.3g}"):
    print(pr.head(18).to_string(index=False))
print(f"\n  atoms with through-space rings: {len(pr)};  "
      f"median |r_KintD|={pr.r_KintD.abs().median():.3f};  "
      f"frac |r|>0.5: {(pr.r_KintD.abs()>0.5).mean():.2f}")
