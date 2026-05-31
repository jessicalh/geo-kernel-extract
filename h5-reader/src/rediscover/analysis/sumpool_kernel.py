#!/usr/bin/env python3
"""sumpool_kernel — recover the ring-current kernel as a function, two honest layers.

LAYER 1 (well-posed): reconstruct the producer's PURE ring-current shielding
B = bare_T0 (a clean, noise-free analytic sum over rings) with a sum-pooling
model fed RAW per-source geometry:

    B(atom, frame) = SUM_sources g(r, cos_theta, intensity)        (no baseline)

If g, fed only raw geometry, reconstructs B and its read-out matches the Pople
kernel (3cos^2-1)/r^3, the FORM has fallen out of the data. All sources are
included (no self-ring guess) -- the read-out shows what g does at the self
geometry, which answers the producer-convention question empirically.

LAYER 2 (physics): how much of the DFT shielding does that kernel actually
explain? Reported by look01 (within-atom B vs DFT). Here we just print the
ceiling for context.
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

# ---- read out g vs the analytic Pople kernel (the equation check) ----
print("\n== g (from raw geometry) vs Pople kernel, on the real source distribution ==")
model.eval()
with torch.no_grad():
    gsrc = model.g(feat).squeeze(-1).cpu().numpy()
r = df.r.to_numpy(); ct = df.cos_theta.to_numpy(); inten = df.ring_intensity.to_numpy()
kern = {"(3cos^2-1)/r^3": (3*ct**2-1)/r**3,
        "intensity*(3cos^2-1)/r^3": inten*(3*ct**2-1)/r**3}

def fit(x, y):
    ok = np.isfinite(x)&np.isfinite(y)
    s, b = np.polyfit(x[ok], y[ok], 1)
    yh = s*x[ok]+b
    return s, np.corrcoef(x[ok], y[ok])[0,1], 1-((y[ok]-yh)**2).sum()/((y[ok]-y[ok].mean())**2).sum()

for name, k in kern.items():
    s, rr, R2 = fit(k, gsrc)
    print(f"  g vs {name:26s}  r={rr:+.4f}  R2={R2:+.4f}  slope={s:+.5g}")

# angular law at fixed r, median intensity
print("\n== angular law of g (r=4 A, median intensity): is it (3cos^2-1)? ==")
cts = np.linspace(-1, 1, 41)
grid = np.stack([np.full_like(cts,(4.0-fmean[0])/fstd[0]), (cts-fmean[1])/fstd[1],
                 np.full_like(cts,(np.median(inten)-fmean[2])/fstd[2])], 1).astype(np.float32)
with torch.no_grad():
    gv = model.g(torch.tensor(grid, device=dev)).squeeze(-1).cpu().numpy()
s, rr, R2 = fit(3*cts**2-1, gv)
print(f"  g(cos) vs (3cos^2-1):  R2={R2:+.4f}  slope={s:+.4g}")

# radial law at fixed cos (theta~0, in-plane) and fixed cos~1 (axial)
print("\n== radial law of g vs 1/r^3 (median intensity) ==")
for label, cval in [("in-plane cos=0", 0.0), ("axial cos=0.9", 0.9)]:
    rs = np.linspace(3.0, 12.0, 40)
    grid = np.stack([(rs-fmean[0])/fstd[0], np.full_like(rs,(cval-fmean[1])/fstd[1]),
                     np.full_like(rs,(np.median(inten)-fmean[2])/fstd[2])], 1).astype(np.float32)
    with torch.no_grad():
        gv = model.g(torch.tensor(grid, device=dev)).squeeze(-1).cpu().numpy()
    s, rr, R2 = fit(1.0/rs**3, gv)
    print(f"  {label:16s}  g vs 1/r^3:  R2={R2:+.4f}  slope={s:+.4g}")
np.savez("/tmp/rediscover-out/sumpool_kernel_readout.npz", cts=cts, gv=gv, gsrc=gsrc,
         r=r, ct=ct, inten=inten)
