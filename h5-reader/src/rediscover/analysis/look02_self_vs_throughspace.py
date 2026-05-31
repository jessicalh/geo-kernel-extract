#!/usr/bin/env python3
"""look02 — is the signal self-ring or genuine through-space ring current?

Codex flagged that our sum includes the aromatic H's OWN ring (the producer
excludes self/bonded). The current CSV has no ring_index, so we separate by
geometry: for an aromatic C-H, its own ring centre is ~2-2.6 A away and roughly
in-plane (small |z|, large rho, dipolar ~ -1/r^3, large and negative). Through-
space rings are farther and vary across frames.

Per atom, per frame: split sources into NEAREST (putative self) and REST
(through-space). Then ask which part carries the within-atom correlation with
DFT sigma_iso. If REST carries it and NEAREST is near-constant, the signal is
genuine through-space ring current AND the self-ring is dead weight to remove.
"""
import sys
import numpy as np
import pandas as pd

SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rediscover-out/ring_current_sources.csv"
AGG = sys.argv[2] if len(sys.argv) > 2 else "/tmp/rediscover-out/ring_current_aggregated.csv"

src = pd.read_csv(SRC)
src = src.rename(columns={"dipolar_3cos2m1_over_r3": "dipolar"})
print(f"loaded {SRC}: {len(src)} source rows, {src.atom_index.nunique()} atoms")

# Nearest source per (atom, frame) = putative self ring.
src = src.sort_values(["atom_index", "h5_row", "r"])
key = ["atom_index", "h5_row"]
src["rank_r"] = src.groupby(key)["r"].rank(method="first")
nearest = src[src.rank_r == 1].copy()

print("\n== nearest-source geometry (putative self ring) ==")
print(f"  r:        median={nearest.r.median():.3f} A   "
      f"5-95%: {nearest.r.quantile(.05):.3f} - {nearest.r.quantile(.95):.3f}")
print(f"  |ring_z|: median={nearest.ring_z.abs().median():.3f} A   "
      f"ring_rho median={nearest.ring_rho.median():.3f} A")
print(f"  dipolar:  median={nearest.dipolar.median():+.4f}   "
      f"(in-plane self ring -> negative, ~ -1/r^3)")
# Per-atom: is the nearest term near-constant across frames? (self-ring signature)
g = nearest.groupby("atom_index")["dipolar"]
near_cv = (g.std() / g.mean().abs()).abs()
print(f"  per-atom nearest-dipolar coeff-of-variation: median={near_cv.median():.3f} "
      f"(small => near-constant => de-means away)")

# Split each (atom,frame) sum into nearest vs rest.
agg = src.groupby(key).agg(
    sum_all=("dipolar", "sum"),
    n_src=("dipolar", "size"),
).reset_index()
near = nearest[key + ["dipolar"]].rename(columns={"dipolar": "sum_near"})
agg = agg.merge(near, on=key)
agg["sum_rest"] = agg["sum_all"] - agg["sum_near"]

# Bring in the DFT target and identity from the aggregated table.
a = pd.read_csv(AGG)
a = a[(a.dft_present == 1)][key + ["dft_sigma_iso", "atom_name", "residue_number"]]
agg = agg.merge(a, on=key, how="inner")
print(f"\n  merged rows: {len(agg)}")

def within_corr(frame, col, target="dft_sigma_iso"):
    f = frame.copy()
    f["x_w"] = f[col] - f.groupby("atom_index")[col].transform("mean")
    f["d_w"] = f[target] - f.groupby("atom_index")[target].transform("mean")
    ok = np.isfinite(f.x_w) & np.isfinite(f.d_w)
    if ok.sum() < 3 or f.x_w[ok].std() == 0:
        return np.nan
    return np.corrcoef(f.x_w[ok], f.d_w[ok])[0, 1]

print("\n== within-atom correlation with DFT sigma_iso, by sum component ==")
print(f"  sum_all   (self + through-space): r={within_corr(agg, 'sum_all'):+.4f}")
print(f"  sum_near  (self ring only):       r={within_corr(agg, 'sum_near'):+.4f}")
print(f"  sum_rest  (through-space only):   r={within_corr(agg, 'sum_rest'):+.4f}")

# Per-atom, how much variance does the nearest (self) carry vs rest?
print("\n== variance contribution across frames (per atom, mean over atoms) ==")
v = agg.groupby("atom_index").agg(
    var_all=("sum_all", "var"), var_near=("sum_near", "var"), var_rest=("sum_rest", "var"))
print(f"  mean var(sum_all)={v.var_all.mean():.5g}  var(sum_near)={v.var_near.mean():.5g}  "
      f"var(sum_rest)={v.var_rest.mean():.5g}")
print(f"  -> the self ring contributes {100*v.var_near.mean()/v.var_all.mean():.0f}% of "
      f"per-atom sum variance on average")
