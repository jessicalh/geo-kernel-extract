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

This is a FIT (DeepSets pooling on C++-emitted raw geometry against C++-emitted
DFT sigma_iso) — no Python physics recompute. The previous "read out g and
regress on (3cos^2-1)/r^3" block is DELETED per
feedback_no_python_physics_except_labeled_integrity_test (2026-06-01) + the
lead's decision: whether g is the Pople form is decided by pysr_distill.py
(symbolic), and any kernel comparison reads the C++-emitted `dipolar` column.
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

# ---- read out the learned g and correlate with the EMITTED kernel column ----
# (no recompute: `dipolar` is the C++-emitted (3cos^2-1)/r^3 per source; the
# symbolic form check is pysr_distill.py's job.)
print("\n== learned g(r, cos_theta, intensity) vs the EMITTED dipolar column ==")
model.eval()
with torch.no_grad():
    g_on_src = model.g(feat).squeeze(-1).cpu().numpy()
dipolar = df.dipolar.to_numpy()              # C++-emitted (3cos^2-1)/r^3
inten = df.ring_intensity.to_numpy()
ok = np.isfinite(g_on_src) & np.isfinite(dipolar)
print(f"  corr(g, emitted dipolar)            r={np.corrcoef(g_on_src[ok], dipolar[ok])[0,1]:+.4f}")
ok2 = ok & np.isfinite(inten)
print(f"  corr(g, emitted intensity*dipolar)  r={np.corrcoef(g_on_src[ok2], (inten*dipolar)[ok2])[0,1]:+.4f}")
print("  (symbolic 'is it the Pople form' is decided by pysr_distill.py, not here)")
