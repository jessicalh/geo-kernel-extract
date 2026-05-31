#!/usr/bin/env python3
"""diag_differencing — does rate-of-change across frame windows aid identifiability?

The identifiability killer is collinearity: when two competing sources move
together across frames, a levels regression can't separate their coefficients.
The hypothesis (instantaneous model, chain rule): the DIFFERENCES decorrelate
even when the levels don't, because different motions have different rates.

Measured on what we have (ring-current cell, v2 substrate). Two gates:
  1. SMOOTHNESS — is the signal smooth across the DFT stride (~20-40 ps) so a
     finite difference is a real derivative, or is it already near-white frame-
     to-frame (then differencing is just noise)? lag-1 autocorrelation answers.
  2. DECORRELATION — for pairs of through-space rings on the same atom, does
     |corr| drop from levels to differences (windows w=1,3,5)?
Plus the noise cost: Var(diff)/(2 Var(level)) ~ 1 - autocorr (>~1 => differencing
roughly doubles relative noise). Report numbers; let the data vote.
"""
import numpy as np
import pandas as pd

SRC = "/tmp/rediscover-out-v2/ring_current_sources.csv"
AGG = "/tmp/rediscover-out-v2/ring_current_aggregated.csv"
df = pd.read_csv(SRC).rename(columns={"dipolar_3cos2m1_over_r3": "dipolar"})
df = df[(df.dft_present == 1) & (df.is_self_or_bonded == 0)].copy()

def autocorr1(x):
    x = x - x.mean()
    d = (x**2).sum()
    return float((x[:-1]*x[1:]).sum()/d) if d > 0 else np.nan

# Build per (atom, ring) dipolar series across DFT frames (sorted by frame).
series = {}  # (atom, ring) -> (frames array, dipolar array)
for (a, rg), g in df.groupby(["atom_index", "ring_index"]):
    g = g.sort_values("h5_row")
    if len(g) >= 50 and g.dipolar.std() > 1e-4:
        series[(a, rg)] = (g.h5_row.to_numpy(), g.dipolar.to_numpy())

# Gate 1: smoothness (lag-1 autocorr of each ring series, and of sigma per atom).
ac_ring = np.array([autocorr1(v[1]) for v in series.values()])
agg = pd.read_csv(AGG); agg = agg[agg.dft_present == 1]
ac_sig = []
for a, g in agg.groupby("atom_index"):
    g = g.sort_values("h5_row")
    if len(g) >= 50:
        ac_sig.append(autocorr1(g.dft_sigma_iso.to_numpy()))
ac_sig = np.array(ac_sig)
print("== Gate 1: smoothness across the DFT stride (lag-1 autocorr) ==")
print(f"  ring dipolar series: median={np.nanmedian(ac_ring):+.3f}  "
      f"(>~0.3 => smooth enough that a difference is a derivative, not noise)")
print(f"  sigma_iso per atom : median={np.nanmedian(ac_sig):+.3f}")
print(f"  noise cost Var(diff)/2Var(level) ~ 1-autocorr: ring median={1-np.nanmedian(ac_ring):.2f}")

# Gate 2: decorrelation of collinear ring pairs (levels vs differences).
def diff_w(frames, vals, w):
    # finite difference over a window of w STRIDED frames where both ends exist
    fr2idx = {f: i for i, f in enumerate(frames)}
    out_a, out_b = [], []
    return frames, vals  # placeholder, replaced below

def aligned_pair(s1, s2):
    f1, v1 = s1; f2, v2 = s2
    common = np.intersect1d(f1, f2)
    if len(common) < 30: return None
    i1 = {f: i for i, f in enumerate(f1)}; i2 = {f: i for i, f in enumerate(f2)}
    return common, np.array([v1[i1[f]] for f in common]), np.array([v2[i2[f]] for f in common])

def corr(a, b):
    if a.std() == 0 or b.std() == 0: return np.nan
    return np.corrcoef(a, b)[0, 1]

def windowed_diff(common, x, w):
    # difference between frames w apart in the (already strided) common index
    return x[w:] - x[:-w]

rows = []
keys_by_atom = {}
for (a, rg) in series: keys_by_atom.setdefault(a, []).append((a, rg))
for a, ks in keys_by_atom.items():
    for i in range(len(ks)):
        for j in range(i+1, len(ks)):
            ap = aligned_pair(series[ks[i]], series[ks[j]])
            if ap is None: continue
            common, x, y = ap
            lev = abs(corr(x, y))
            d1 = abs(corr(windowed_diff(common, x, 1), windowed_diff(common, y, 1)))
            d3 = abs(corr(windowed_diff(common, x, 3), windowed_diff(common, y, 3)))
            d5 = abs(corr(windowed_diff(common, x, 5), windowed_diff(common, y, 5)))
            rows.append((a, lev, d1, d3, d5))
pp = pd.DataFrame(rows, columns=["atom", "lev", "d1", "d3", "d5"]).dropna()
print(f"\n== Gate 2: collinear ring-pair |corr|, levels vs differences ({len(pp)} pairs) ==")
print(f"  median |corr|  levels={pp.lev.median():.3f}  diff w=1={pp.d1.median():.3f}  "
      f"w=3={pp.d3.median():.3f}  w=5={pp.d5.median():.3f}")
hi = pp[pp.lev > 0.5]
print(f"  among strongly-collinear level pairs (|corr|>0.5, n={len(hi)}):")
print(f"    median |corr| levels={hi.lev.median():.3f} -> diff w=1={hi.d1.median():.3f} "
      f"w=3={hi.d3.median():.3f} w=5={hi.d5.median():.3f}")
print(f"    fraction of those decorrelated below 0.3 by w=1: {(hi.d1<0.3).mean():.2f}")

# Conditioning of the per-ring-type feature matrix: levels vs differences.
tcols = [f"sum_dipolar_ringtype_{t}" for t in range(8)]
A = agg.sort_values(["atom_index","h5_row"])
def maxoffdiag(M):
    C = np.corrcoef(M, rowvar=False); np.fill_diagonal(C, 0)
    return np.nanmax(np.abs(C))
lev_offd, dif_offd = [], []
for a, g in A.groupby("atom_index"):
    X = g[tcols].to_numpy()
    keep = X.std(0) > 1e-5
    if keep.sum() < 2: continue
    Xk = X[:, keep]
    lev_offd.append(maxoffdiag(Xk))
    dif_offd.append(maxoffdiag(np.diff(Xk, axis=0)))
print(f"\n== per-ring-type feature matrix max off-diagonal |corr| (atoms with >=2 active types) ==")
print(f"  levels median={np.nanmedian(lev_offd):.3f}   differences median={np.nanmedian(dif_offd):.3f}  "
      f"(n_atoms={len(lev_offd)})")
print("\n(read: differencing helps identifiability iff Gate-1 smooth AND Gate-2 |corr| drops)")
