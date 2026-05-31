#!/usr/bin/env python3
"""credibility2 — the instantaneous framing, tested across ATOMS (no time).

Ring current is instantaneous: sigma = baseline(atom) + k * kernel(geometry).
The trajectory is a geometry SAMPLER, not a process. So the strong, autocorr-
free generalization test is leave-atoms-out: fit the universal coefficient k on
a set of atoms (all their snapshots), then predict the within-atom ring-current
modulation of atoms HELD OUT from the fit. Different atoms are independent
samples, so time/autocorrelation is irrelevant here.

baseline(atom) is a per-atom nuisance removed by de-meaning each atom with its
OWN mean (no leakage of k). The transferable claim is the single shared k.
"""
import numpy as np
import pandas as pd
import os
import sys

out_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("REDISCOVER_OUT", "/tmp/rediscover-out")
df = pd.read_csv(f"{out_dir}/ring_current_aggregated.csv")
df = df[(df.dft_present == 1)].copy()
PRED = "sum_dipolar_producer_valid" if "sum_dipolar_producer_valid" in df.columns else "sum_dipolar"
TGT = "dft_sigma_iso"

# per-atom de-mean (baseline removal, per atom, uses only that atom's own data)
df["xw"] = df[PRED] - df.groupby("atom_index")[PRED].transform("mean")
df["yw"] = df[TGT] - df.groupby("atom_index")[TGT].transform("mean")

atoms = df.atom_index.unique()
sd = df.groupby("atom_index")[PRED].std()
coupled = set(sd[sd > 0.012].index)   # atoms with real ring-current modulation
print(f"atoms={len(atoms)}  coupled (sd_pred>0.012)={len(coupled)}")

def r2(y, yhat):
    return 1 - ((y-yhat)**2).sum()/((y-y.mean())**2).sum()

rng = np.random.default_rng(0)
all_R2, coup_R2, ks = [], [], []
for trial in range(40):
    perm = rng.permutation(atoms)
    train_a, test_a = set(perm[:len(perm)//2]), set(perm[len(perm)//2:])
    tr = df[df.atom_index.isin(train_a)]
    te = df[df.atom_index.isin(test_a)]
    k, c = np.polyfit(tr.xw, tr.yw, 1)           # universal coefficient on TRAIN atoms
    ks.append(k)
    all_R2.append(r2(te.yw.to_numpy(), k*te.xw.to_numpy()))      # held-out atoms
    tec = te[te.atom_index.isin(coupled)]
    if len(tec) > 10:
        coup_R2.append(r2(tec.yw.to_numpy(), k*tec.xw.to_numpy()))

print("\n== leave-ATOMS-out (universal k fit on train atoms, predict held-out atoms) ==")
print(f"  shared coefficient k: median={np.median(ks):+.3g}  spread="
      f"[{np.percentile(ks,10):+.3g}, {np.percentile(ks,90):+.3g}] ppm/(A^-3)")
print(f"  held-out-atom within-atom R2 (all test atoms):     median={np.median(all_R2):+.3f}")
print(f"  held-out-atom within-atom R2 (coupled test atoms): median={np.median(coup_R2):+.3f}")
print(f"  (positive => one universal coefficient transfers to atoms it never saw;")
print(f"   autocorrelation is irrelevant here — atoms are independent units)")

# For contrast: a single global fit on ALL atoms (the descriptive coefficient).
kg, cg = np.polyfit(df.xw, df.yw, 1)
print(f"\n  global descriptive k (all atoms): {kg:+.3g} ppm/(A^-3)  "
      f"in-sample within-atom R2={r2(df.yw.to_numpy(), kg*df.xw.to_numpy()):+.3f}")
