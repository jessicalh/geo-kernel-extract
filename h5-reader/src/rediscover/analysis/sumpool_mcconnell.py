#!/usr/bin/env python3
"""sumpool_mcconnell — does the bond-anisotropy (McConnell) kernel also fall out?

Same pooling recovery as the ring case, on backbone amide HN. Per anisotropic
bond source, the contribution is modelled as

    chi[category] * h(r, cos_theta_bond_axis)

  - h(r, cos): a shared geometric MLP fed RAW geometry (never the precomputed
    dipolar). Whether it recovers the bond-anisotropy form is read against the
    C++-emitted `dipolar` column / PySR, not by recomputing (3cos^2-1)/r^3.
  - chi[category]: a learned per-bond-category scalar (peptide C=O, C-N,
    sidechain C=O, aromatic) standing in for the anisotropy Delta-chi -- the
    parameter the McConnell form puts the chemistry in.

Target is the producer's pure McConnell kernel B = bare_T0 (clean, analytic,
EMITTED). This is a FIT on emitted columns. Per
feedback_no_python_physics_except_labeled_integrity_test (2026-06-01) + the
lead's decision, the (3cos^2-1)/r^3 readout-comparison arrays are DELETED; the
learned h is correlated with the C++-emitted `dipolar` column instead.
"""
import numpy as np
import os
import pandas as pd
import sys
import torch
import torch.nn as nn

torch.manual_seed(0); np.random.seed(0)
dev = "cuda" if torch.cuda.is_available() else "cpu"
out_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("REDISCOVER_OUT", "/tmp/rediscover-out")
SRC = f"{out_dir}/mcconnell_sources.csv"

df = pd.read_csv(SRC).rename(columns={"dipolar_3cos2m1_over_r3": "dipolar",
                                      "cos_theta_bond_axis": "ct"})
df = df[df.bare_kernel_present == 1].copy()
df["cat"], cats = pd.factorize(df.bond_category)
n_cat = df.cat.max() + 1
print(f"source rows: {len(df)}  HN atoms={df.atom_index.nunique()}  bond categories={list(cats)}")

df["gid"], _ = pd.factorize(df.atom_index.astype(str) + ":" + df.h5_row.astype(str))
n_groups = df.gid.max() + 1
grp = df.groupby("gid").agg(B=("bare_T0", "first"), frame=("h5_row", "first")).sort_index()
target = torch.tensor(grp.B.to_numpy(), dtype=torch.float32, device=dev)

frames = np.sort(grp.frame.unique())
test_frames = set(np.random.default_rng(0).choice(frames, int(0.2*len(frames)), replace=False))
is_test = grp.frame.isin(test_frames).to_numpy()
tr = torch.tensor(~is_test, device=dev); te = torch.tensor(is_test, device=dev)

feat_raw = df[["r", "ct"]].to_numpy(np.float32)
fmean, fstd = feat_raw.mean(0), feat_raw.std(0)
feat = torch.tensor((feat_raw-fmean)/fstd, dtype=torch.float32, device=dev)
cat = torch.tensor(df.cat.to_numpy(), dtype=torch.long, device=dev)
src_group = torch.tensor(df.gid.to_numpy(), dtype=torch.long, device=dev)


class McPool(nn.Module):
    def __init__(self, n_cat):
        super().__init__()
        self.chi = nn.Embedding(n_cat, 1)
        nn.init.ones_(self.chi.weight)
        self.h = nn.Sequential(nn.Linear(2, 64), nn.SiLU(),
                               nn.Linear(64, 64), nn.SiLU(), nn.Linear(64, 1))

    def forward(self, feat, cat, src_group, n_groups):
        c = (self.chi(cat).squeeze(-1) * self.h(feat).squeeze(-1))
        pooled = torch.zeros(n_groups, device=feat.device)
        return pooled.scatter_add_(0, src_group, c)


def r2(p, t): return float(1 - ((t-p)**2).sum()/((t-t.mean())**2).sum())


model = McPool(n_cat).to(dev)
opt = torch.optim.Adam(model.parameters(), lr=3e-3, weight_decay=1e-7)
for ep in range(5000):
    model.train(); opt.zero_grad()
    pred = model(feat, cat, src_group, n_groups)
    loss = ((pred[tr]-target[tr])**2).mean()
    loss.backward(); opt.step()
    if ep % 1000 == 0 or ep == 4999:
        model.eval()
        with torch.no_grad():
            pr = model(feat, cat, src_group, n_groups)
        print(f"  ep {ep:4d}  MSE={loss.item():.5f}  R2 train={r2(pr[tr],target[tr]):+.4f}  "
              f"test={r2(pr[te],target[te]):+.4f}")

print("\n== learned per-category chi (relative anisotropy weights) ==")
with torch.no_grad():
    chi = model.chi.weight.squeeze(-1).cpu().numpy()
for c, name in zip(chi, cats):
    print(f"  category {name}: chi={c:+.4g}")

print("\n== read out learned h vs the EMITTED dipolar column ==")
with torch.no_grad():
    hsrc = model.h(feat).squeeze(-1).cpu().numpy()
r = df.r.to_numpy()
dipolar = df.dipolar.to_numpy()          # C++-emitted (3cos^2-1)/r^3 per source
far = r >= 3.0
ok = far & np.isfinite(hsrc) & np.isfinite(dipolar)
s, b = np.polyfit(dipolar[ok], hsrc[ok], 1)
rr = np.corrcoef(dipolar[ok], hsrc[ok])[0, 1]
R2 = 1 - ((hsrc[ok]-(s*dipolar[ok]+b))**2).sum()/((hsrc[ok]-hsrc[ok].mean())**2).sum()
print(f"  h vs emitted dipolar (r>=3): r={rr:+.4f}  R2={R2:+.4f}  slope={s:+.4g}")
print("  (symbolic form recovery is pysr_distill.py's job; no kernel recompute here)")
