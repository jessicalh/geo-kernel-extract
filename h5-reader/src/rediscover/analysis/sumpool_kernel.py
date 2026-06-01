#!/usr/bin/env python3
"""sumpool_kernel — recover the ring-current kernel as a FIT (no Python physics).

Reconstruct the producer's PURE ring-current shielding B = bare_T0 (a clean,
noise-free analytic per-atom kernel, EMITTED by the C++ producer) with a
sum-pooling model fed RAW per-source geometry:

    B(atom, frame) = SUM_sources g(r, cos_theta, intensity)        (no baseline)

This is a FIT, not a recompute: the inputs (r, cos_theta, ring_intensity) and the
target (bare_T0) are all C++-emitted columns. The learned per-source function g
is written out for symbolic distillation by pysr_distill.py; whether g matches
the Pople form is decided THERE (PySR), not by recomputing (3cos^2-1)/r^3 in
numpy. Per feedback_no_python_physics_except_labeled_integrity_test (2026-06-01
tightening) and the lead's decision: the (3cos^2-1)/r^3 recompute arrays are
DELETED; any kernel comparison reads the C++-emitted `dipolar` column or relies
on PySR. The only surviving recompute in this workspace is the pinned
change-of-basis fixture test.
"""
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

torch.manual_seed(0); np.random.seed(0)
dev = "cuda" if torch.cuda.is_available() else "cpu"
SRC = "/tmp/rediscover-out/ring_current_sources.csv"

df = pd.read_csv(SRC).rename(columns={"dipolar_3cos2m1_over_r3": "dipolar"})
df = df[df.bare_kernel_present == 1].copy()           # B is always present; keep all rings
print(f"source rows: {len(df)}  atoms={df.atom_index.nunique()}  (all rings, self included)")

df["gid"], _ = pd.factorize(df.atom_index.astype(str) + ":" + df.h5_row.astype(str))
n_groups = df.gid.max() + 1
grp = df.groupby("gid").agg(B=("bare_T0", "first"), frame=("h5_row", "first")).sort_index()
target = torch.tensor(grp.B.to_numpy(), dtype=torch.float32, device=dev)

frames = np.sort(grp.frame.unique())
test_frames = set(np.random.default_rng(0).choice(frames, int(0.2*len(frames)), replace=False))
is_test = grp.frame.isin(test_frames).to_numpy()
tr = torch.tensor(~is_test, device=dev); te = torch.tensor(is_test, device=dev)

feat_raw = df[["r", "cos_theta", "ring_intensity"]].to_numpy(np.float32)
fmean, fstd = feat_raw.mean(0), feat_raw.std(0)
feat = torch.tensor((feat_raw - fmean)/fstd, dtype=torch.float32, device=dev)
src_group = torch.tensor(df.gid.to_numpy(), dtype=torch.long, device=dev)


class SumPool(nn.Module):
    def __init__(self):
        super().__init__()
        self.g = nn.Sequential(nn.Linear(3, 64), nn.SiLU(),
                               nn.Linear(64, 64), nn.SiLU(), nn.Linear(64, 1))

    def forward(self, feat, src_group, n_groups):
        c = self.g(feat).squeeze(-1)
        pooled = torch.zeros(n_groups, device=feat.device)
        return pooled.scatter_add_(0, src_group, c)


def r2(pred, tgt):
    return float(1 - ((tgt-pred)**2).sum()/((tgt-tgt.mean())**2).sum())


model = SumPool().to(dev)
opt = torch.optim.Adam(model.parameters(), lr=3e-3, weight_decay=1e-7)
for ep in range(5000):
    model.train(); opt.zero_grad()
    pred = model(feat, src_group, n_groups)
    loss = ((pred[tr]-target[tr])**2).mean()
    loss.backward(); opt.step()
    if ep % 1000 == 0 or ep == 4999:
        model.eval()
        with torch.no_grad():
            pr = model(feat, src_group, n_groups)
        print(f"  ep {ep:4d}  train MSE={loss.item():.5f}  "
              f"R2 train={r2(pr[tr],target[tr]):+.4f}  test={r2(pr[te],target[te]):+.4f}")

# ---- read out the learned per-source g and SAVE it for PySR distillation ----
# No Python kernel recompute: we emit the learned g together with the C++-emitted
# geometry columns (r, cos_theta, ring_intensity) and the emitted `dipolar`
# column. pysr_distill.py decides symbolically whether g is the Pople form;
# refine_kernel.py draws g against the emitted `dipolar`. We do NOT build
# (3cos^2-1)/r^3 here.
print("\n== reading out learned g(r, cos_theta, intensity) for PySR distillation ==")
model.eval()
with torch.no_grad():
    gsrc = model.g(feat).squeeze(-1).cpu().numpy()
r = df.r.to_numpy(); ct = df.cos_theta.to_numpy(); inten = df.ring_intensity.to_numpy()
dipolar = df.dipolar.to_numpy()  # C++-emitted (3cos^2-1)/r^3 per source

# correlation of the learned g with the EMITTED kernel column (read, not recompute)
ok = np.isfinite(gsrc) & np.isfinite(dipolar)
r_emit = np.corrcoef(gsrc[ok], dipolar[ok])[0, 1]
r_emit_int = np.corrcoef(gsrc[ok], (inten * dipolar)[ok])[0, 1]
print(f"  corr(learned g, emitted dipolar)            r={r_emit:+.4f}")
print(f"  corr(learned g, emitted intensity*dipolar)  r={r_emit_int:+.4f}")
print("  (symbolic form recovery is done by pysr_distill.py, not by recomputing the kernel)")
np.savez("/tmp/rediscover-out/sumpool_kernel_readout.npz", gsrc=gsrc,
         r=r, ct=ct, inten=inten, dipolar=dipolar)
