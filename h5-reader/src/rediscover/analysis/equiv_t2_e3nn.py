#!/usr/bin/env python3
"""equiv_t2_e3nn — the ring-current T2 equivariant fitter, on e3nn (replaces the
hand-rolled equiv_t2.py and its numpy lib_T2 end-run).

The model is the same physics ansatz as before — a permutation-invariant sum over
through-space ring sources of an l=2 contribution oriented by the displacement
direction r_hat and the source ring normal n_hat, gated by an invariant radial
function of (r, intensity) — but EVERY equivariant feature is now computed by
e3nn (o3.spherical_harmonics for the pure-vector l=2, an o3 tensor product for
the r_hat (x) n_hat cross term), so equivariance is the library's (Wigner-D
tested) and there is NO numpy spherical-harmonic projection anywhere.

  contribution_2e(source) = w_A(r,i) * Y2(r_hat)            # o3.spherical_harmonics
                          + w_B(r,i) * Y2(n_hat)            # o3.spherical_harmonics
                          + cross_term(r_hat, n_hat)        # o3 1o (x) 1o -> 2e
  node_2e = SUM_sources contribution_2e                     # scatter-add (Deep Sets)

THE BOUNDARY (feedback_no_python_physics_except_labeled_integrity_test +
feedback_model_is_spine):
  * Inputs are READ from the C++ export: per-source disp_local (vec3),
    source_normal_local (vec3), invariant scalars [r, ring_intensity], and the
    (atom, frame) group index. No kernel, no projection, no field is recomputed.
  * The TARGET is the C++-emitted library-basis local-frame DFT T2 NPY sidecar
    `rediscover_ring_current_sources_target_local_T2.npy`. It is REQUIRED — the
    fitter fails loud if absent; there is NO Python projection fallback (we never
    recompute the library 5-vector from the emitted 3x3). The target is mapped
    into e3nn's 2e convention by the PINNED constant change_of_basis.C (a single
    matmul on emitted data, NOT a projection). |T2| is invariant under that
    orthogonal map, so the reported magnitude correlation is identical in either
    basis.
  * e3nn's internal (y,z,x) axis convention is handled once, in change_of_basis.

Honesty (unchanged from equiv_t2.py): de-mean target AND prediction per atom
(strip the near-constant local C-H baseline tensor); frame-split is the GATE;
leave-atoms-out is reported (not gated, thin at ~7 coupled aromatic H); report
|T2| (rotation invariant) and the 5-vector R2.

Cross-term mode (decision 4: test exact vs learnable, pick better, report both):
  --cross exact      : o3.TensorProduct(1o,1o->2e) with FIXED unit path weight
                       (matches the hand-roll's symmetric outer-product ansatz)
  --cross learnable  : o3.FullyConnectedTensorProduct(1o,1o->2e), weights learned
  --cross both       : train both, print both, exit with the better test R2 model

Run (system python; torch cu130 needs the nvidia cu13 libs on LD_LIBRARY_PATH or
it segfaults — see requirements-e3nn.txt / ENV.md):
    LD_LIBRARY_PATH=<cu13 libs>:<torch/lib> python3 equiv_t2_e3nn.py [out_dir] --cross both
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from e3nn import o3

import change_of_basis as cob

torch.manual_seed(0)
np.random.seed(0)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?",
                    default=os.environ.get("REDISCOVER_OUT", "/tmp/rediscover-out-v2"))
    ap.add_argument("--cross", choices=["exact", "learnable", "both"], default="both")
    ap.add_argument("--epochs", type=int, default=4000)
    ap.add_argument("--lr", type=float, default=3e-3)
    return ap.parse_args()


# ── Load the substrate. READ-ONLY: every column/array below is C++-emitted. ───
def load(out_dir):
    src_csv = f"{out_dir}/ring_current_sources.csv"
    df = pd.read_csv(src_csv)
    df = df.rename(columns={"dipolar_3cos2m1_over_r3": "dipolar"})

    # The TARGET is the C++-emitted local-frame library-basis T2 NPY sidecar,
    # row-aligned with the sources CSV (RecordSink writes one NPY row per source
    # row, in order). We FAIL LOUD if it is absent — there is NO Python
    # projection fallback. The producer already wrote this 5-vector; recomputing
    # it from the emitted 3x3 in numpy would be the projection end-run the law
    # forbids (feedback_no_python_physics_except_labeled_integrity_test). If the
    # NPY is missing, re-run the extractor for this out_dir (the sidecar is part
    # of RecordSink's output); do not synthesise it here.
    npy = f"{out_dir}/rediscover_ring_current_sources_target_local_T2.npy"
    if not os.path.exists(npy):
        sys.exit(
            f"FATAL: required emitted target NPY not found:\n  {npy}\n"
            "The equivariant fitter reads the C++-emitted local-frame T2 target; "
            "it does NOT recompute the spherical projection. Re-run the extractor "
            "so the *_sources_target_local_T2.npy sidecar is present in this dir.")
    target_lib_all = np.load(npy)  # (n_rows_all, 5) library basis, emitted
    if len(target_lib_all) != len(df):
        sys.exit(f"FATAL: NPY rows {len(target_lib_all)} != CSV rows {len(df)} "
                 "(sidecar/CSV out of sync; re-extract).")
    df["__rowid"] = np.arange(len(df))
    return df, target_lib_all


class EquivPoolE3nn(nn.Module):
    """e3nn equivariant sum-pool: per-source 2e from spherical harmonics + a
    1o(x)1o->2e cross term, gated by an invariant radial MLP, scatter-pooled to
    (atom, frame) groups, de-meaned per atom."""

    def __init__(self, n_groups, n_atoms, cross="exact"):
        super().__init__()
        self.n_groups = n_groups
        self.n_atoms = n_atoms
        self.cross = cross
        # invariant radial MLP -> 2 gates (Y2(r_hat), Y2(n_hat)) + 1 cross gate
        self.R = nn.Sequential(
            nn.Linear(2, 64), nn.SiLU(),
            nn.Linear(64, 64), nn.SiLU(),
            nn.Linear(64, 3))
        irr_1o = o3.Irreps("1o")
        irr_2e = o3.Irreps("2e")
        if cross == "exact":
            # fixed-path symmetric product 1o (x) 1o -> 2e (matches the hand-roll's
            # Y2sym(r,n)); no internal weights, gated externally by R.
            self.tp = o3.TensorProduct(
                irr_1o, irr_1o, irr_2e,
                instructions=[(0, 0, 0, "uuu", False)],
                irrep_normalization="component")
        else:  # learnable
            self.tp = o3.FullyConnectedTensorProduct(
                irr_1o, irr_1o, irr_2e, irrep_normalization="component")

    def forward(self, rhat, nhat, featn, src_group, group_atom):
        w = self.R(featn)                                   # (n_src, 3) gates
        y2_r = o3.spherical_harmonics(2, rhat, normalize=True, normalization="component")
        y2_n = o3.spherical_harmonics(2, nhat, normalize=True, normalization="component")
        cross = self.tp(rhat, nhat)                         # (n_src, 5) e3nn 2e
        contrib = (w[:, 0:1] * y2_r
                   + w[:, 1:2] * y2_n
                   + w[:, 2:3] * cross)                     # (n_src, 5)
        pooled = torch.zeros(self.n_groups, 5, device=contrib.device)
        pooled.index_add_(0, src_group, contrib)
        # de-mean prediction per atom (mirror the target's per-atom de-meaning)
        mean = torch.zeros(self.n_atoms, 5, device=contrib.device)
        cnt = torch.zeros(self.n_atoms, 1, device=contrib.device)
        mean.index_add_(0, group_atom, pooled)
        cnt.index_add_(0, group_atom, torch.ones(self.n_groups, 1, device=contrib.device))
        return pooled - (mean / cnt)[group_atom]


def r2(p, y):
    return float(1 - ((y - p) ** 2).sum() / ((y - y.mean(0)) ** 2).sum())


def build_tensors(df, target_lib_all, dev):
    # through-space sources only (self/bonded is the local baseline; de-meaned)
    keep = (df.dft_present == 1) & (df.is_self_or_bonded == 0) & (df.dft_local_frame_valid == 1)
    d = df.loc[keep].copy()
    print(f"through-space source rows: {len(d)}  atoms={d.atom_index.nunique()}")

    disp = d[["disp_local_x", "disp_local_y", "disp_local_z"]].to_numpy(float)
    rhat = disp / np.linalg.norm(disp, axis=1, keepdims=True)
    nrm = d[["source_normal_local_x", "source_normal_local_y",
             "source_normal_local_z"]].to_numpy(float)
    nhat = nrm / np.clip(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-9, None)

    feat = d[["r", "ring_intensity"]].to_numpy(np.float32)
    fmean, fstd = feat.mean(0), feat.std(0)
    featn = (feat - fmean) / fstd

    # per-source emitted library target -> e3nn 2e via the pinned constant C.
    tgt_lib = target_lib_all[d["__rowid"].to_numpy()]          # (n_src, 5) emitted
    # group by (atom, frame); the target is per-group, take groupby-first
    d = d.assign(__gid=pd.factorize(d.atom_index.astype(str) + ":" + d.h5_row.astype(str))[0])
    for k in range(5):
        d[f"__t{k}"] = tgt_lib[:, k]
    n_groups = int(d["__gid"].max()) + 1
    grp = d.groupby("__gid").agg(aid=("atom_index", "first"),
                                 frame=("h5_row", "first")).sort_index()
    gT2_lib = d.groupby("__gid")[[f"__t{k}" for k in range(5)]].first().sort_index().to_numpy()
    gT2_e3nn = cob.lib_to_e3nn(gT2_lib)                        # constant matmul

    aid_of_group = grp.aid.to_numpy()
    group_atom_idx, _ = pd.factorize(aid_of_group)
    n_atoms = int(group_atom_idx.max()) + 1

    # de-mean target per atom
    gT2_dm = gT2_e3nn.copy()
    for a in np.unique(group_atom_idx):
        m = group_atom_idx == a
        gT2_dm[m] -= gT2_e3nn[m].mean(0)

    # frame split (gate); also keep atom ids for the leave-atoms-out report
    frames = np.sort(grp.frame.unique())
    test_f = set(np.random.default_rng(0).choice(frames, int(0.2 * len(frames)), replace=False))
    g_te = grp.frame.isin(test_f).to_numpy()

    t = lambda x, dt=torch.float32: torch.tensor(x, dtype=dt, device=dev)
    pack = dict(
        rhat=t(rhat), nhat=t(nhat), featn=t(featn),
        src_group=t(d["__gid"].to_numpy(), torch.long),
        group_atom=t(group_atom_idx, torch.long),
        target=t(gT2_dm),
        g_te=t(g_te, torch.bool),
        n_groups=n_groups, n_atoms=n_atoms,
        aid_of_group=aid_of_group, group_atom_idx=group_atom_idx,
        gT2_e3nn=gT2_e3nn,
    )
    return pack


def train_one(pack, dev, cross, epochs, lr):
    model = EquivPoolE3nn(pack["n_groups"], pack["n_atoms"], cross=cross).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-7)
    target = pack["target"]
    g_te = pack["g_te"]
    g_tr = ~g_te
    for ep in range(epochs):
        model.train(); opt.zero_grad()
        pred = model(pack["rhat"], pack["nhat"], pack["featn"],
                     pack["src_group"], pack["group_atom"])
        loss = ((pred[g_tr] - target[g_tr]) ** 2).mean()
        loss.backward(); opt.step()
        if ep % 1000 == 0 or ep == epochs - 1:
            model.eval()
            with torch.no_grad():
                pr = model(pack["rhat"], pack["nhat"], pack["featn"],
                           pack["src_group"], pack["group_atom"])
            print(f"  [{cross:9s}] ep {ep:4d} loss={loss.item():.4e} "
                  f"T2 R2 train={r2(pr[g_tr], target[g_tr]):+.3f} "
                  f"test={r2(pr[g_te], target[g_te]):+.3f}")
    model.eval()
    with torch.no_grad():
        pr = model(pack["rhat"], pack["nhat"], pack["featn"],
                   pack["src_group"], pack["group_atom"]).cpu().numpy()
    tg = target.cpu().numpy()
    te = g_te.cpu().numpy()
    mag_p, mag_t = np.linalg.norm(pr, axis=1), np.linalg.norm(tg, axis=1)
    r_mag = np.corrcoef(mag_p[te], mag_t[te])[0, 1]
    r2_comp = r2(torch.tensor(pr[te]), torch.tensor(tg[te]))
    print(f"  [{cross:9s}] frame-split: T2 5-vec R2(test)={r2_comp:+.3f}  "
          f"|T2| corr(test) r={r_mag:+.3f}")

    # leave-atoms-out (reported, NOT gated): refit gates per held-out atom set.
    loo = leave_atoms_out(pack, dev, cross, epochs, lr)
    print(f"  [{cross:9s}] leave-atoms-out: T2 5-vec R2={loo:+.3f} "
          f"(reported, not gated; thin ~7 coupled aromatic H)")
    return r2_comp, r_mag, loo


def leave_atoms_out(pack, dev, cross, epochs, lr):
    """Reported honesty check: predict the T2 of atoms held entirely out of the
    fit. Thin (~7 coupled aromatic H); not the gate."""
    atoms = np.unique(pack["group_atom_idx"])
    if len(atoms) < 3:
        return float("nan")
    pred = np.full(pack["gT2_e3nn"].shape, np.nan)
    ga = pack["group_atom_idx"]
    for held in atoms:
        tr_mask = ga != held
        te_mask = ga == held
        # train on the other atoms (shorter schedule; this is a report, not gate)
        model = EquivPoolE3nn(pack["n_groups"], pack["n_atoms"], cross=cross).to(dev)
        opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-7)
        gtr = torch.tensor(tr_mask, dtype=torch.bool, device=dev)
        for ep in range(max(800, epochs // 4)):
            model.train(); opt.zero_grad()
            p = model(pack["rhat"], pack["nhat"], pack["featn"],
                      pack["src_group"], pack["group_atom"])
            loss = ((p[gtr] - pack["target"][gtr]) ** 2).mean()
            loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            p = model(pack["rhat"], pack["nhat"], pack["featn"],
                      pack["src_group"], pack["group_atom"]).cpu().numpy()
        pred[te_mask] = p[te_mask]
    # de-mean per atom for a fair within-atom modulation R2
    tg = pack["gT2_e3nn"].copy()
    for a in atoms:
        m = ga == a
        tg[m] -= pack["gT2_e3nn"][m].mean(0)
        pred[m] -= pred[m].mean(0)
    ok = np.isfinite(pred).all(1)
    return r2(torch.tensor(pred[ok]), torch.tensor(tg[ok]))


def main():
    args = parse_args()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={dev}  out_dir={args.out_dir}  cross={args.cross}")
    # pin/derive C up front so a convention mismatch fails loudly before training
    C = cob.get_C()
    print(f"change-of-basis C orthogonality |C^T C - I|max={np.abs(C.T @ C - np.eye(5)).max():.2e}")

    df, target_lib_all = load(args.out_dir)
    pack = build_tensors(df, target_lib_all, dev)
    print(f"groups={pack['n_groups']} atoms={pack['n_atoms']}")

    results = {}
    modes = ["exact", "learnable"] if args.cross == "both" else [args.cross]
    for m in modes:
        print(f"\n== cross-term: {m} ==")
        results[m] = train_one(pack, dev, m, args.epochs, args.lr)

    if len(results) > 1:
        best = max(results, key=lambda k: results[k][0])  # by frame-split R2
        print(f"\nBETTER cross-term by frame-split T2 R2: {best} "
              f"(R2={results[best][0]:+.3f}, |T2| r={results[best][1]:+.3f})")
    print("\nGATE: frame-split T2 R2 (baseline to reproduce/beat: 0.467; |T2| r: 0.756)")


if __name__ == "__main__":
    main()
