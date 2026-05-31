#!/usr/bin/env python3
"""sumpool_t0 — does the ring-current kernel fall out of a flexible sum-pooling model?

Model (permutation-invariant Deep Sets, the scalar precursor to the e3nn T2 model):

    sigma_iso(atom, frame) = baseline[atom] + SUM_sources g(r, cos_theta, intensity)

  - g is a small MLP fed RAW per-source geometry only -- never the precomputed
    dipolar. If g reconstructs (3cos^2-1)/r^3 it rediscovered the Pople kernel
    from geometry; we do not hand it the answer.
  - baseline[atom] is a learned per-atom intercept: it absorbs the chemical
    baseline and the near-constant self ring, so g is left to fit the
    geometric MODULATION (the only thing the kernel claims).

Honesty guards:
  - frame-wise train/test split: g must predict the shielding modulation on
    HELD-OUT frames, not memorise.
  - report WITHIN-ATOM R^2 (the modulation), not the easy between-atom R^2.
  - read out g and regress it on the analytic Pople kernel: that R^2 + slope is
    the literal "did the equation fall out" number.
"""
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

torch.manual_seed(0)
np.random.seed(0)
dev = "cuda" if torch.cuda.is_available() else "cpu"

SRC = "/tmp/rediscover-out/ring_current_sources.csv"
SELF_R = 3.0  # self/bonded ring proxy (identity-based exclusion comes with the C++ re-run)

df = pd.read_csv(SRC).rename(columns={"dipolar_3cos2m1_over_r3": "dipolar"})
df = df[(df.dft_present == 1) & (df.r >= SELF_R)].copy()
print(f"through-space source rows with DFT: {len(df)}  atoms={df.atom_index.nunique()}")

# Contiguous group (atom,frame) and atom ids.
df["gid"], _ = pd.factorize(df.atom_index.astype(str) + ":" + df.h5_row.astype(str))
df["aid"], _ = pd.factorize(df.atom_index)
n_groups, n_atoms = df.gid.max() + 1, df.aid.max() + 1

# group -> atom, group -> target sigma_iso, group -> frame (for the split)
grp = df.groupby("gid").agg(aid=("aid", "first"), sig=("dft_sigma_iso", "first"),
                            frame=("h5_row", "first")).sort_index()
group_atom = torch.tensor(grp.aid.to_numpy(), dtype=torch.long, device=dev)
target = torch.tensor(grp.sig.to_numpy(), dtype=torch.float32, device=dev)

# Frame-wise split: 80% of frames train, 20% test.
frames = np.sort(grp.frame.unique())
rng = np.random.default_rng(0)
test_frames = set(rng.choice(frames, size=int(0.2 * len(frames)), replace=False))
grp_is_test = grp.frame.isin(test_frames).to_numpy()
g_train = torch.tensor(~grp_is_test, device=dev)
g_test = torch.tensor(grp_is_test, device=dev)

# Per-source raw features (standardised) + the source->group index.
feat_raw = df[["r", "cos_theta", "ring_intensity"]].to_numpy(np.float32)
fmean, fstd = feat_raw.mean(0), feat_raw.std(0)
feat = torch.tensor((feat_raw - fmean) / fstd, dtype=torch.float32, device=dev)
src_group = torch.tensor(df.gid.to_numpy(), dtype=torch.long, device=dev)
print(f"groups(atom,frame)={n_groups}  atoms={n_atoms}  test frames={len(test_frames)}/{len(frames)}")


class SumPool(nn.Module):
    def __init__(self, n_atoms):
        super().__init__()
        self.baseline = nn.Embedding(n_atoms, 1)
        nn.init.zeros_(self.baseline.weight)
        self.g = nn.Sequential(nn.Linear(3, 64), nn.SiLU(),
                               nn.Linear(64, 64), nn.SiLU(), nn.Linear(64, 1))

    def forward(self, feat, src_group, group_atom, n_groups):
        contrib = self.g(feat).squeeze(-1)                       # per source
        pooled = torch.zeros(n_groups, device=feat.device)
        pooled.scatter_add_(0, src_group, contrib)               # SUM over sources
        return pooled + self.baseline(group_atom).squeeze(-1)


model = SumPool(n_atoms).to(dev)
opt = torch.optim.Adam(model.parameters(), lr=5e-3, weight_decay=1e-6)


def within_atom_r2(pred, tgt, atom):
    # de-mean both by atom, on the given subset, then 1 - SS_res/SS_tot.
    pred, tgt, atom = pred.detach(), tgt.detach(), atom.detach()
    out = {}
    for a in atom.unique():
        m = atom == a
        if m.sum() < 2:
            continue
        out[int(a)] = (pred[m] - pred[m].mean(), tgt[m] - tgt[m].mean())
    p = torch.cat([v[0] for v in out.values()])
    t = torch.cat([v[1] for v in out.values()])
    ss_res = ((p - t) ** 2).sum() if False else ((t - p) ** 2).sum()
    ss_tot = (t ** 2).sum()
    return float(1 - ss_res / ss_tot)


for ep in range(4000):
    model.train(); opt.zero_grad()
    pred = model(feat, src_group, group_atom, n_groups)
    loss = ((pred[g_train] - target[g_train]) ** 2).mean()
    loss.backward(); opt.step()
    if ep % 800 == 0 or ep == 3999:
        model.eval()
        with torch.no_grad():
            pr = model(feat, src_group, group_atom, n_groups)
            r2_tr = within_atom_r2(pr[g_train], target[g_train], group_atom[g_train])
            r2_te = within_atom_r2(pr[g_test], target[g_test], group_atom[g_test])
        print(f"  ep {ep:4d}  train MSE={loss.item():.4f}  "
              f"within-atom R2 train={r2_tr:+.3f}  test={r2_te:+.3f}")

# ---- read out g and compare to the analytic Pople kernel ----
print("\n== read out g(r, cos_theta, intensity) vs the Pople kernel ==")
model.eval()
with torch.no_grad():
    g_on_src = model.g(feat).squeeze(-1).cpu().numpy()
r = df.r.to_numpy(); ct = df.cos_theta.to_numpy(); inten = df.ring_intensity.to_numpy()
pople_geo = (3 * ct**2 - 1) / r**3            # pure geometry
pople_int = inten * pople_geo                 # intensity-weighted (physical form)

def fit_r2(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    b = np.polyfit(x[ok], y[ok], 1)
    yhat = np.polyval(b, x[ok])
    r2 = 1 - ((y[ok]-yhat)**2).sum()/((y[ok]-y[ok].mean())**2).sum()
    return b[0], r2, np.corrcoef(x[ok], y[ok])[0, 1]

for name, kern in [("(3cos^2-1)/r^3", pople_geo), ("intensity*(3cos^2-1)/r^3", pople_int)]:
    slope, r2, rr = fit_r2(kern, g_on_src)
    print(f"  g  vs  {name:26s}  r={rr:+.4f}  R2={r2:+.4f}  slope={slope:+.4g}")

# cos_theta shape at fixed r, median intensity: is it (3cos^2-1)?
print("\n== angular shape of g at r=4.0 A, median intensity ==")
med_i = (np.median(inten) - fmean[2]) / fstd[2]
r_fix = (4.0 - fmean[0]) / fstd[0]
cts = np.linspace(-1, 1, 41)
grid = np.stack([np.full_like(cts, r_fix), (cts - fmean[1]) / fstd[1],
                 np.full_like(cts, med_i)], 1).astype(np.float32)
with torch.no_grad():
    gv = model.g(torch.tensor(grid, device=dev)).squeeze(-1).cpu().numpy()
slope, r2, rr = fit_r2(3 * cts**2 - 1, gv)
print(f"  g(cos) vs (3cos^2-1):  R2={r2:+.4f}  slope={slope:+.4g}  "
      f"(R2~1 => the angular law is the Pople (3cos^2-1))")
np.savez("/tmp/rediscover-out/sumpool_readout.npz", cts=cts, gv=gv,
         g_on_src=g_on_src, pople_geo=pople_geo, pople_int=pople_int)
print("  saved readout -> /tmp/rediscover-out/sumpool_readout.npz")
