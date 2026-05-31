#!/usr/bin/env python3
"""credibility_check — is the ring-current result indicative, or a setup?

Separates the circular part from the real part, with numbers:

  CIRCULAR (sanity only): recovering the Pople form by fitting the producer
  kernel bare_T0 reverse-engineers the producer's own formula. Not tested here
  (we already know it works); flagged in the report.

  REAL TEST: does the geometry predict the INDEPENDENT DFT shielding
  OUT-OF-SAMPLE? Fit a coefficient on TRAIN frames, predict HELD-OUT frames.
  - shared single coefficient across all atoms (a universal ring-current
    constant predicting unseen frames is the strong claim);
  - per-atom out-of-sample R^2 for the strongly-coupled atoms.
  Plus effective-N via the MD autocorrelation time, so the per-atom numbers
  are not over-read.
"""
import numpy as np
import pandas as pd

AGG = "/tmp/rediscover-out/ring_current_aggregated.csv"
df = pd.read_csv(AGG)
df = df[(df.dft_present == 1) & (df.bare_kernel_present == 1)].copy()
# within-atom modulation is through-space dominated (self ring ~constant, de-means away)
PRED = "sum_dipolar"   # our geometry sum; non-circular predictor of DFT
TGT = "dft_sigma_iso"

# frame split (by frame, no leakage)
frames = np.sort(df.h5_row.unique())
rng = np.random.default_rng(0)
test_f = set(rng.choice(frames, int(0.3*len(frames)), replace=False))
df["is_test"] = df.h5_row.isin(test_f)

# de-mean by atom using TRAIN means only (no leakage from test into the centering)
tr = df[~df.is_test]
mu_x = tr.groupby("atom_index")[PRED].mean()
mu_y = tr.groupby("atom_index")[TGT].mean()
df["xw"] = df[PRED] - df.atom_index.map(mu_x)
df["yw"] = df[TGT] - df.atom_index.map(mu_y)
df = df.dropna(subset=["xw", "yw"])

def r2(y, yhat):
    return 1 - ((y-yhat)**2).sum()/((y-y.mean())**2).sum()

# ---- shared single coefficient, fit on train, evaluate out-of-sample ----
a = df[~df.is_test]; b = df[df.is_test]
k, c = np.polyfit(a.xw, a.yw, 1)
oos = r2(b.yw.to_numpy(), k*b.xw.to_numpy()+c)
ins = r2(a.yw.to_numpy(), k*a.xw.to_numpy()+c)
print("== shared single ring-current coefficient (geometry -> DFT) ==")
print(f"  fitted on train frames:  k={k:+.4g} ppm/(A^-3)   intercept={c:+.4g}")
print(f"  in-sample within-atom R2:     {ins:+.4f}")
print(f"  OUT-OF-SAMPLE within-atom R2: {oos:+.4f}   <-- the non-circular number")
print(f"  (a universal coefficient predicting unseen frames; n_atoms={df.atom_index.nunique()})")

# ---- per-atom out-of-sample, for the strongly-coupled atoms ----
print("\n== per-atom OUT-OF-SAMPLE R2 (fit k_atom on train, predict test) ==")
rows = []
for aid, g in df.groupby("atom_index"):
    gtr, gte = g[~g.is_test], g[g.is_test]
    if len(gtr) < 20 or len(gte) < 10 or gtr.xw.std() == 0:
        continue
    ka, ca = np.polyfit(gtr.xw, gtr.yw, 1)
    oos_a = r2(gte.yw.to_numpy(), ka*gte.xw.to_numpy()+ca)
    # effective N via lag-1 autocorrelation of the DFT series (train+test)
    y = g.sort_values("h5_row")[TGT].to_numpy()
    yd = y - y.mean()
    rho = (yd[:-1]*yd[1:]).sum()/(yd**2).sum() if (yd**2).sum() > 0 else 0
    neff = len(y)*(1-rho)/(1+rho) if rho < 1 else len(y)
    rows.append((aid, g.atom_name.iloc[0], int(g.residue_number.iloc[0]),
                 g[PRED].std(), oos_a, ka, rho, neff))
pa = pd.DataFrame(rows, columns=["atom_index","atom_name","resnum","sd_pred",
                                 "oos_R2","k_atom","ac_lag1","n_eff"])
pa = pa.sort_values("sd_pred", ascending=False)
with pd.option_context("display.float_format", lambda v: f"{v:+.3g}"):
    print(pa.head(14).to_string(index=False))
strong = pa.head(10)
print(f"\n  top-10 most-modulated atoms: median OUT-OF-SAMPLE R2 = {strong.oos_R2.median():+.3f}")
print(f"  their k_atom: median={strong.k_atom.median():+.3g}  spread="
      f"[{strong.k_atom.min():+.3g}, {strong.k_atom.max():+.3g}] ppm/(A^-3)")
print(f"  MD autocorrelation lag-1: median={pa.ac_lag1.median():.3f}  "
      f"=> median n_eff={pa.n_eff.median():.0f} (vs {df.groupby('atom_index').size().median():.0f} frames)")
