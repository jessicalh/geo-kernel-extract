#!/usr/bin/env python3
"""equiv_t2_e3nn — the ring-current T2 equivariant fitter, on e3nn (replaces the
hand-rolled equiv_t2.py and its numpy lib_T2 end-run).

================================================================================
PARITY DEFECT — STALE FITTER, DO NOT COPY THE CROSS TERM (note 2026-06-13)
================================================================================
NATURE: the `cross_term(r_hat, n_hat)` below is built as an o3 TensorProduct
  `1o (x) 1o -> 2e`, declaring BOTH inputs polar. But n_hat is the ring NORMAL,
  an AXIAL pseudovector (parity EVEN, irrep `1e`) — a cross product of two
  in-plane polar vectors. r_hat is polar (`1o`). So an axial vector is being fed
  through a polar slot.
EXTENT / why it is wrong: the parity-valid product of a polar r_hat (1o) and an
  axial n_hat (1e) is `1o (x) 1e -> 2o` (ODD). The shielding T2 target is EVEN
  (`2e`) — because B is axial and sigma maps axial->axial (B = nabla x A; first
  principles, verified 2026-06-13). So this cross term silently injects a
  parity-FORBIDDEN (2o-flavoured) contribution labelled 2e: it PASSES
  rotation-only equivariance checks and BREAKS under reflection/inversion. There
  is NO parity-valid `2e` bilinear from (polar)(x)(axial). The only even ways to
  couple r_hat and n_hat are: Y2(r_hat), Y2(n_hat) (each already here, each `2e`),
  or gating those by the even invariant (r_hat . n_hat)^2 (`0e`). Note
  (r_hat . n_hat) alone is a PSEUDOscalar (`0o`).
WHY NOT FIX/EXTEND HERE: this fitter is STALE and superseded. The parity-correct
  pattern already exists in equiv_t2_backbone_e3nn.py — ring normals labelled
  `1e`, orientation carried by Y2(axis) (= traceless axis(x)axis, the `2e`
  director), and NO normal/axis ever fed through a `1o` cross path. A fix here is
  a redesign (drop the cross term — Y2(r)+Y2(n) already supply the valid `2e` — or
  use the (r.n)^2 even-gate), to be done ONLY if this fitter is revived. It is
  not retrofitted now, and the scope of this defect is exactly this one cross
  term in this one stale file — nothing else here is implicated.
DO NOT GENERALISE: do not copy this `1o (x) 1o -> 2e` cross term into a new model
  (e.g. the equivariant side model, EQUIVARIANT_SIDE_MODEL_SPEC_2026-06-13.md), and
  do not extend the pattern to other axial inputs (B-field, shielding T1, bond
  directors). Build on equiv_t2_backbone_e3nn.py, and make a rotation + REFLECTION
  equivariance test a build gate — a reflection test catches exactly this defect;
  rotation-only checks miss it.
================================================================================

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

Honesty (EFG-path protocol): blocked/purged frame-split is the GATE; target,
prediction, and invariant-feature normalization are centered/scaled from TRAIN
rows only. Leave-atoms-out is reported (not gated, thin at ~7 coupled aromatic H);
report |T2| (rotation invariant) and the 5-vector R2.

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
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from e3nn import o3

import change_of_basis as cob
from e3nn_protocol import centred_by_train_atom, make_split_masks, normalised_by_train_rows

torch.manual_seed(0)
np.random.seed(0)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?",
                    default=os.environ.get("REDISCOVER_OUT", "/tmp/rediscover-out-v2"))
    ap.add_argument("--cross", choices=["exact", "learnable", "both"], default="both")
    ap.add_argument("--epochs", type=int, default=4000)
    ap.add_argument("--lr", type=float, default=3e-3)
    ap.add_argument("--artifact-dir", default=None,
                    help="write predictions, metrics, run audit, and advisor graph here")
    ap.add_argument("--t0-loss-weight", type=float, default=0.25,
                    help="same-model T0 companion loss weight after train-only scaling")
    ap.add_argument("--loao", action="store_true",
                    help="also run the leave-atoms-out audit. It is not the graph gate.")
    ap.add_argument("--no-graph", action="store_true",
                    help="skip rendering the PNG/PDF graph, but still write metrics/predictions")
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
    (atom, frame) groups, de-meaned per atom from the supplied centering mask.
    A scalar 0e/T0 companion head is summed over the same emitted source rows; it
    is not a separate scalar model."""

    def __init__(self, n_groups, n_atoms, cross="exact"):
        super().__init__()
        self.n_groups = n_groups
        self.n_atoms = n_atoms
        self.cross = cross
        # invariant radial MLP -> 2 gates (Y2(r_hat), Y2(n_hat)) + 1 cross gate
        self.R2 = nn.Sequential(
            nn.Linear(2, 64), nn.SiLU(),
            nn.Linear(64, 64), nn.SiLU(),
            nn.Linear(64, 3))
        self.R0 = nn.Sequential(
            nn.Linear(2, 64), nn.SiLU(),
            nn.Linear(64, 64), nn.SiLU(),
            nn.Linear(64, 1))
        # >>> PARITY DEFECT (see module docstring) — DO NOT COPY THIS CROSS TERM <<<
        # This `1o (x) 1o -> 2e` treats n_hat (the AXIAL ring normal, `1e`) as polar.
        # (polar r_hat) (x) (axial n_hat) is parity `2o`, FORBIDDEN for the even
        # (`2e`) shielding T2. Stale/superseded — use equiv_t2_backbone_e3nn.py's
        # Y2(axis) / `2e`-director pattern instead.
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

    def forward(self, rhat, nhat, featn, src_group, group_atom, center_mask=None):
        w = self.R2(featn)                                  # (n_src, 3) gates
        y2_r = o3.spherical_harmonics(2, rhat, normalize=True, normalization="component")
        y2_n = o3.spherical_harmonics(2, nhat, normalize=True, normalization="component")
        cross = self.tp(rhat, nhat)  # PARITY DEFECT (see docstring): nhat is axial 1e, not 1o
        contrib2 = (w[:, 0:1] * y2_r
                    + w[:, 1:2] * y2_n
                    + w[:, 2:3] * cross)                    # (n_src, 5)
        contrib0 = self.R0(featn)                            # (n_src, 1), 0e
        pooled2 = torch.zeros(self.n_groups, 5, device=contrib2.device)
        pooled0 = torch.zeros(self.n_groups, 1, device=contrib2.device)
        pooled2.index_add_(0, src_group, contrib2)
        pooled0.index_add_(0, src_group, contrib0)
        # De-mean prediction per atom using train rows only; this mirrors the
        # target centering and keeps held-out/purged groups out of the mean.
        if center_mask is None:
            center_mask = torch.ones(group_atom.shape, dtype=torch.bool, device=group_atom.device)
        mean2 = torch.zeros(self.n_atoms, 5, device=contrib2.device)
        mean0 = torch.zeros(self.n_atoms, 1, device=contrib2.device)
        cnt = torch.zeros(self.n_atoms, 1, device=contrib2.device)
        center_atom = group_atom[center_mask]
        mean2.index_add_(0, center_atom, pooled2[center_mask])
        mean0.index_add_(0, center_atom, pooled0[center_mask])
        cnt.index_add_(0, center_atom,
                       torch.ones(center_atom.shape[0], 1, device=contrib2.device))
        centered2 = pooled2 - (mean2 / cnt.clamp_min(1.0))[group_atom]
        centered0 = pooled0 - (mean0 / cnt.clamp_min(1.0))[group_atom]
        return centered2, centered0.squeeze(1)


def r2(p, y):
    return float(1 - ((y - p) ** 2).sum() / ((y - y.mean(0)) ** 2).sum())


def r2_np(p, y):
    p = np.asarray(p, dtype=float)
    y = np.asarray(y, dtype=float)
    den = np.sum((y - y.mean(axis=0)) ** 2)
    if den <= 0:
        return float("nan")
    return float(1.0 - np.sum((y - p) ** 2) / den)


def corr_np(a, b):
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if len(a) < 2 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def train_atom_means(x, group_atom_idx, train_mask):
    x = np.asarray(x, dtype=float)
    group_atom_idx = np.asarray(group_atom_idx)
    train_mask = np.asarray(train_mask, dtype=bool)
    means = np.full_like(x, np.nan, dtype=float)
    for atom in np.unique(group_atom_idx):
        m = group_atom_idx == atom
        tr = m & train_mask
        if tr.sum() == 0:
            continue
        means[m] = x[tr].mean(axis=0, keepdims=True)
    return means


def split_labels(train_mask, test_mask):
    out = np.full(len(train_mask), "purged", dtype=object)
    out[np.asarray(train_mask, dtype=bool)] = "train"
    out[np.asarray(test_mask, dtype=bool)] = "test"
    return out


def disk_audit_for_write(path):
    path = Path(path).resolve()
    base = path if path.exists() else path.parent
    usage = shutil.disk_usage(base)
    min_free = 20 * 1024 ** 3
    if usage.free < min_free:
        raise SystemExit(
            f"FATAL: disk gate failed before write: {base} has "
            f"{usage.free / 1024 ** 3:.1f} GiB free < 20 GiB")
    runs = Path("/tmp/rediscover-runs")
    total = 0
    if runs.exists():
        for root, dirs, files in os.walk(runs):
            for name in files:
                try:
                    total += (Path(root) / name).stat().st_size
                except OSError:
                    pass
    max_total = 15 * 1024 ** 3
    if total > max_total:
        raise SystemExit(
            f"FATAL: rediscover run budget failed before write: "
            f"{total / 1024 ** 3:.2f} GiB > 15 GiB")
    return {
        "path": str(path),
        "free_bytes": int(usage.free),
        "free_gib": usage.free / 1024 ** 3,
        "rediscover_runs_bytes": int(total),
        "rediscover_runs_gib": total / 1024 ** 3,
    }


def build_tensors(df, target_lib_all, dev):
    # through-space sources only (self/bonded is the local baseline; de-meaned)
    keep = (df.dft_present == 1) & (df.is_self_or_bonded == 0) & (df.dft_local_frame_valid == 1)
    d = df.loc[keep].copy()
    finite_cols = [
        "disp_local_x", "disp_local_y", "disp_local_z",
        "source_normal_local_x", "source_normal_local_y", "source_normal_local_z",
        "r", "ring_intensity", "dft_sigma_iso",
    ]
    finite = np.isfinite(d[finite_cols].to_numpy(float)).all(axis=1)
    if not finite.all():
        d = d.loc[finite].copy()
    print(f"through-space source rows: {len(d)}  atoms={d.atom_index.nunique()}")

    disp = d[["disp_local_x", "disp_local_y", "disp_local_z"]].to_numpy(float)
    rhat = disp / np.linalg.norm(disp, axis=1, keepdims=True)
    nrm = d[["source_normal_local_x", "source_normal_local_y",
             "source_normal_local_z"]].to_numpy(float)
    nhat = nrm / np.clip(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-9, None)

    feat = d[["r", "ring_intensity"]].to_numpy(np.float32)

    # per-source emitted library target -> e3nn 2e via the pinned constant C.
    tgt_lib = target_lib_all[d["__rowid"].to_numpy()]          # (n_src, 5) emitted
    # group by (atom, frame); the target is per-group, take groupby-first
    d = d.assign(__gid=pd.factorize(d.atom_index.astype(str) + ":" + d.h5_row.astype(str))[0])
    for k in range(5):
        d[f"__t{k}"] = tgt_lib[:, k]
    n_groups = int(d["__gid"].max()) + 1
    grp = d.groupby("__gid").agg(
        aid=("atom_index", "first"),
        residue_number=("residue_number", "first"),
        atom_name=("atom_name", "first"),
        frame=("h5_row", "first"),
        original_index=("original_index", "first"),
        time_ps=("time_ps", "first"),
        n_sources=("__gid", "size"),
        target_T0=("dft_sigma_iso", "first"),
    ).sort_index()
    gT2_lib = d.groupby("__gid")[[f"__t{k}" for k in range(5)]].first().sort_index().to_numpy()
    gT2_e3nn = cob.lib_to_e3nn(gT2_lib)                        # constant matmul
    gT0 = grp.target_T0.to_numpy(float)

    aid_of_group = grp.aid.to_numpy()
    group_atom_idx, _ = pd.factorize(aid_of_group)
    n_atoms = int(group_atom_idx.max()) + 1

    # EFG-path frame gate: contiguous held-out block with neighbouring train
    # frames purged. Target and feature normalization are train-only.
    g_tr, g_te, split_info = make_split_masks(grp.frame.to_numpy())
    gT2_dm = centred_by_train_atom(gT2_e3nn, group_atom_idx, g_tr)
    gT0_dm = centred_by_train_atom(gT0[:, None], group_atom_idx, g_tr)[:, 0]
    gT2_mean = train_atom_means(gT2_e3nn, group_atom_idx, g_tr)
    gT0_mean = train_atom_means(gT0[:, None], group_atom_idx, g_tr)[:, 0]
    source_train = g_tr[d["__gid"].to_numpy()]
    featn, fmean, fstd = normalised_by_train_rows(feat, source_train)

    t = lambda x, dt=torch.float32: torch.tensor(x, dtype=dt, device=dev)
    pack = dict(
        rhat=t(rhat), nhat=t(nhat), featn=t(featn),
        src_group=t(d["__gid"].to_numpy(), torch.long),
        group_atom=t(group_atom_idx, torch.long),
        target_t2=t(gT2_dm),
        target_t0=t(gT0_dm),
        g_te=t(g_te, torch.bool),
        g_tr=t(g_tr, torch.bool),
        n_groups=n_groups, n_atoms=n_atoms,
        aid_of_group=aid_of_group, group_atom_idx=group_atom_idx,
        gT2_e3nn=gT2_e3nn, gT0=gT0, gT2_mean=gT2_mean, gT0_mean=gT0_mean,
        group_frame=grp.frame.to_numpy(),
        group_original_index=grp.original_index.to_numpy(),
        group_time_ps=grp.time_ps.to_numpy(float),
        group_residue_number=grp.residue_number.to_numpy(),
        group_atom_name=grp.atom_name.to_numpy(),
        group_source_count=grp.n_sources.to_numpy(),
        split_label=split_labels(g_tr, g_te),
        feature_mean=fmean, feature_scale=fstd,
        source_rows_all=len(df),
        source_rows_kept=len(d),
        source_columns=len(df.columns),
        source_dense_bytes_estimate=int(len(df) * len(df.columns) * 8),
        **split_info,
    )
    return pack


def train_one(pack, dev, cross, epochs, lr, t0_loss_weight=0.25, run_loao=False):
    model = EquivPoolE3nn(pack["n_groups"], pack["n_atoms"], cross=cross).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-7)
    target_t2 = pack["target_t2"]
    target_t0 = pack["target_t0"]
    g_te = pack["g_te"]
    g_tr = pack["g_tr"]
    t2_scale = target_t2[g_tr].std().clamp_min(1e-6)
    t0_scale = target_t0[g_tr].std().clamp_min(1e-6)
    for ep in range(epochs):
        model.train(); opt.zero_grad()
        pred_t2, pred_t0 = model(pack["rhat"], pack["nhat"], pack["featn"],
                                 pack["src_group"], pack["group_atom"], center_mask=g_tr)
        loss_t2 = (((pred_t2[g_tr] - target_t2[g_tr]) / t2_scale) ** 2).mean()
        loss_t0 = (((pred_t0[g_tr] - target_t0[g_tr]) / t0_scale) ** 2).mean()
        loss = loss_t2 + t0_loss_weight * loss_t0
        loss.backward(); opt.step()
        if ep % 1000 == 0 or ep == epochs - 1:
            model.eval()
            with torch.no_grad():
                pr2, pr0 = model(pack["rhat"], pack["nhat"], pack["featn"],
                                 pack["src_group"], pack["group_atom"], center_mask=g_tr)
            print(f"  [{cross:9s}] ep {ep:4d} loss={loss.item():.4e} "
                  f"T2 R2 train={r2(pr2[g_tr], target_t2[g_tr]):+.3f} "
                  f"test={r2(pr2[g_te], target_t2[g_te]):+.3f} "
                  f"T0 test={r2(pr0[g_te], target_t0[g_te]):+.3f}")
    model.eval()
    with torch.no_grad():
        pr_t2_t, pr_t0_t = model(pack["rhat"], pack["nhat"], pack["featn"],
                                 pack["src_group"], pack["group_atom"], center_mask=g_tr)
    pr = pr_t2_t.cpu().numpy()
    pr0 = pr_t0_t.cpu().numpy()
    tg = target_t2.cpu().numpy()
    tg0 = target_t0.cpu().numpy()
    te = g_te.cpu().numpy()
    tr = g_tr.cpu().numpy()
    mag_p, mag_t = np.linalg.norm(pr, axis=1), np.linalg.norm(tg, axis=1)
    r_mag = corr_np(mag_p[te], mag_t[te])
    r2_comp = r2_np(pr[te], tg[te])
    r2_train = r2_np(pr[tr], tg[tr])
    r2_mag = r2_np(mag_p[te], mag_t[te])
    r2_t0 = r2_np(pr0[te], tg0[te])
    pred_t0_restored = pr0 + pack["gT0_mean"]
    r2_t0_restored = r2_np(pred_t0_restored[te], pack["gT0"][te])
    print(f"  [{cross:9s}] frame-split: T2 5-vec R2(test)={r2_comp:+.3f}  "
          f"|T2| corr(test) r={r_mag:+.3f}  T0 modulation R2(test)={r2_t0:+.3f}")

    # leave-atoms-out (reported, NOT gated): refit gates per held-out atom set.
    loo = float("nan")
    if run_loao:
        loo = leave_atoms_out(pack, dev, cross, epochs, lr)
        print(f"  [{cross:9s}] leave-atoms-out: T2 5-vec R2={loo:+.3f} "
              f"(reported, not gated; thin coupled aromatic H)")
    return {
        "cross": cross,
        "t2_r2_test": r2_comp,
        "t2_r2_train": r2_train,
        "t2_mag_r_test": r_mag,
        "t2_mag_r2_test": r2_mag,
        "t0_mod_r2_test": r2_t0,
        "t0_restored_r2_test": r2_t0_restored,
        "loao_t2_r2": loo,
        "pred_t2": pr,
        "pred_t0": pr0,
    }


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
            p, _ = model(pack["rhat"], pack["nhat"], pack["featn"],
                         pack["src_group"], pack["group_atom"], center_mask=gtr)
            loss = ((p[gtr] - pack["target_t2"][gtr]) ** 2).mean()
            loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            p, _ = model(pack["rhat"], pack["nhat"], pack["featn"],
                         pack["src_group"], pack["group_atom"], center_mask=gtr)
            p = p.cpu().numpy()
        pred[te_mask] = p[te_mask]
    tg = pack["target_t2"].cpu().numpy()
    ok = np.isfinite(pred).all(1) & np.isfinite(tg).all(1)
    if ok.sum() < 2:
        return float("nan")
    return r2(torch.tensor(pred[ok]), torch.tensor(tg[ok]))


def write_predictions(path, pack, result):
    pred_t2 = result["pred_t2"]
    pred_t0 = result["pred_t0"]
    target_t2 = pack["target_t2"].cpu().numpy()
    target_t0 = pack["target_t0"].cpu().numpy()
    pred_t0_restored = pred_t0 + pack["gT0_mean"]
    rows = {
        "atom_index": pack["aid_of_group"],
        "residue_number": pack["group_residue_number"],
        "atom_name": pack["group_atom_name"],
        "h5_row": pack["group_frame"],
        "original_frame_index": pack["group_original_index"],
        "time_ps": pack["group_time_ps"],
        "split": pack["split_label"],
        "source_rows_in_group": pack["group_source_count"],
        "target_T0_centered": target_t0,
        "pred_T0_centered": pred_t0,
        "target_T0_restored": pack["gT0"],
        "pred_T0_restored": pred_t0_restored,
    }
    for k in range(5):
        rows[f"target_T2_{k}_centered"] = target_t2[:, k]
        rows[f"pred_T2_{k}_centered"] = pred_t2[:, k]
    pd.DataFrame(rows).to_csv(path, index=False)


def write_json(path, payload):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=True)
        fh.write("\n")


def scatter_limits(*arrays):
    vals = np.concatenate([np.asarray(a, dtype=float).ravel() for a in arrays])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (-1.0, 1.0)
    lo, hi = np.percentile(vals, [0.5, 99.5])
    if lo == hi:
        pad = max(1.0, abs(lo) * 0.1)
    else:
        pad = 0.08 * (hi - lo)
    return float(lo - pad), float(hi + pad)


def render_graph(path_png, path_pdf, pack, result, metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    te = pack["g_te"].cpu().numpy()
    target_t2 = pack["target_t2"].cpu().numpy()
    target_t0 = pack["target_t0"].cpu().numpy()
    pred_t2 = result["pred_t2"]
    pred_t0 = result["pred_t0"]
    mag_t = np.linalg.norm(target_t2, axis=1)
    mag_p = np.linalg.norm(pred_t2, axis=1)
    t0_restored = pack["gT0"]
    p0_restored = pred_t0 + pack["gT0_mean"]

    fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.5), constrained_layout=True)
    ax = axs[0, 0]
    x = target_t2[te].ravel()
    y = pred_t2[te].ravel()
    lim = scatter_limits(x, y)
    ax.scatter(x, y, s=9, alpha=0.36, linewidths=0, color="#1f5a78")
    ax.plot(lim, lim, color="black", lw=1.0, alpha=0.75)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title("A. Held-out T2 components")
    ax.set_xlabel("DFT T2 modulation, train-atom centered")
    ax.set_ylabel("e3nn prediction")
    ax.text(0.04, 0.96, f"clean R2 = {metrics['t2_r2_test']:+.3f}",
            transform=ax.transAxes, va="top", ha="left",
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9))

    ax = axs[0, 1]
    lim = scatter_limits(mag_t[te], mag_p[te])
    ax.scatter(mag_t[te], mag_p[te], s=16, alpha=0.45, linewidths=0, color="#7d5a1f")
    ax.plot(lim, lim, color="black", lw=1.0, alpha=0.75)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title("B. |T2| modulation recovery")
    ax.set_xlabel("DFT |T2| modulation")
    ax.set_ylabel("predicted |T2| modulation")
    ax.text(0.04, 0.96,
            f"r = {metrics['t2_mag_r_test']:+.3f}\nmag R2 = {metrics['t2_mag_r2_test']:+.3f}",
            transform=ax.transAxes, va="top", ha="left",
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9))

    ax = axs[1, 0]
    lim = scatter_limits(target_t0[te], pred_t0[te])
    ax.scatter(target_t0[te], pred_t0[te], s=16, alpha=0.45, linewidths=0, color="#3a6f44")
    ax.plot(lim, lim, color="black", lw=1.0, alpha=0.75)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title("C. sigma_iso/T0 companion")
    ax.set_xlabel("DFT T0 modulation, train-atom centered")
    ax.set_ylabel("predicted T0 modulation")
    ax.text(0.04, 0.96, f"MODULATION R2 = {metrics['t0_mod_r2_test']:+.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=11,
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.95))
    iax = inset_axes(ax, width="38%", height="38%", loc="lower right", borderpad=1.2)
    ilim = scatter_limits(t0_restored[te], p0_restored[te])
    iax.scatter(t0_restored[te], p0_restored[te], s=5, alpha=0.35, linewidths=0,
                color="#6b6b6b")
    iax.plot(ilim, ilim, color="black", lw=0.7, alpha=0.6)
    iax.set_xlim(ilim); iax.set_ylim(ilim)
    iax.set_xticks([]); iax.set_yticks([])
    iax.set_title("restored inset", fontsize=8)
    ax.text(0.98, 0.03,
            f"restored R2 = {metrics['t0_restored_r2_test']:+.3f}\n"
            "baseline-dominated inset; not skill metric",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="0.85", alpha=0.9))

    ax = axs[1, 1]
    ax.axis("off")
    bucket = metrics["interpretation_bucket"]
    text = (
        "D. Direction, not destination\n\n"
        "1P9J geometry sampler -> emitted C++ source rows ->\n"
        "e3nn source-sum predictor -> held-out DFT modulation\n\n"
        f"Model: source-e3nn v0 ({result['cross']} cross), 5-component T2.\n"
        f"Protocol: clean blocked/purged split; lag-1 train/test pairs = "
        f"{metrics['cross_split_lag1_pairs']}.\n"
        f"Support: {metrics['source_rows_kept']:,} through-space source rows, "
        f"{metrics['groups']:,} atom-frame groups, {metrics['atoms']} atoms.\n"
        f"Bucket: {bucket}.\n\n"
        "Within-axis 1P9J interpolation only. The trajectory samples\n"
        "instantaneous geometries; this is not a dynamics/process model.\n"
        "The graph shows correlate-not-match recovery of held-out DFT\n"
        "modulation and is a direction signal toward Stage-3, not a\n"
        "transferability or BMRB validation claim."
    )
    ax.text(0.02, 0.98, text, va="top", ha="left", fontsize=10,
            linespacing=1.25, family="DejaVu Sans")

    fig.suptitle("1P9J Source-e3nn Interpolation Advisor Graph", fontsize=14)
    fig.savefig(path_png, dpi=220)
    fig.savefig(path_pdf)
    plt.close(fig)


def interpretation_bucket(t2_r2):
    if not np.isfinite(t2_r2) or t2_r2 < 0.05:
        return "negative/NO-GO graph"
    if t2_r2 < 0.25:
        return "weak direction signal"
    return "direction signal"


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
    print(f"split={pack['split_strategy']} test_frames={pack['test_frames']} "
          f"purged_frames={pack['purged_train_frames']} "
          f"cross_split_lag1_pairs={pack['cross_split_lag1_pairs']}")

    results = {}
    modes = ["exact", "learnable"] if args.cross == "both" else [args.cross]
    for m in modes:
        print(f"\n== cross-term: {m} ==")
        results[m] = train_one(pack, dev, m, args.epochs, args.lr,
                               t0_loss_weight=args.t0_loss_weight,
                               run_loao=args.loao)

    best = max(results, key=lambda k: results[k]["t2_r2_test"])  # by clean held-out T2 R2
    best_result = results[best]
    print(f"\nBETTER cross-term by clean held-out T2 R2: {best} "
          f"(R2={best_result['t2_r2_test']:+.3f}, "
          f"|T2| r={best_result['t2_mag_r_test']:+.3f}, "
          f"T0 mod R2={best_result['t0_mod_r2_test']:+.3f})")
    print("GATE: clean blocked/purged held-out T2 R2 is the build's baseline. "
          "Prior 0.466/0.757 is pre-protocol-fix leaky-suspect and is not used.")

    artifact_dir = Path(args.artifact_dir) if args.artifact_dir else None
    if artifact_dir:
        disk_audit = disk_audit_for_write(artifact_dir)
        artifact_dir.mkdir(parents=True, exist_ok=True)
        pred_csv = artifact_dir / "interp_1p9j_predictions.csv"
        metrics_json = artifact_dir / "interp_1p9j_metrics.json"
        audit_json = artifact_dir / "interp_1p9j_run_audit.json"
        png = artifact_dir / "interp_1p9j_advisor_graph.png"
        pdf = artifact_dir / "interp_1p9j_advisor_graph.pdf"

        te = pack["g_te"].cpu().numpy()
        tr = pack["g_tr"].cpu().numpy()
        metrics = {
            "model_mode": "source-e3nn v0",
            "input_dir": str(Path(args.out_dir).resolve()),
            "best_cross": best,
            "all_cross_metrics": {
                k: {mk: mv for mk, mv in v.items() if not mk.startswith("pred_")}
                for k, v in results.items()
            },
            "t2_r2_test": best_result["t2_r2_test"],
            "t2_r2_train": best_result["t2_r2_train"],
            "t2_mag_r_test": best_result["t2_mag_r_test"],
            "t2_mag_r2_test": best_result["t2_mag_r2_test"],
            "t0_mod_r2_test": best_result["t0_mod_r2_test"],
            "t0_restored_r2_test": best_result["t0_restored_r2_test"],
            "interpretation_bucket": interpretation_bucket(best_result["t2_r2_test"]),
            "protocol": "clean_blocked_purged_train_only_centering",
            "split_strategy": pack["split_strategy"],
            "test_frames": int(pack["test_frames"]),
            "purged_train_frames": int(pack["purged_train_frames"]),
            "cross_split_lag1_pairs": int(pack["cross_split_lag1_pairs"]),
            "groups": int(pack["n_groups"]),
            "train_groups": int(tr.sum()),
            "test_groups": int(te.sum()),
            "purged_groups": int((~tr & ~te).sum()),
            "atoms": int(pack["n_atoms"]),
            "frames": int(len(np.unique(pack["group_frame"]))),
            "source_rows_all": int(pack["source_rows_all"]),
            "source_rows_kept": int(pack["source_rows_kept"]),
            "source_columns": int(pack["source_columns"]),
            "source_dense_bytes_estimate": int(pack["source_dense_bytes_estimate"]),
            "held_out_rows_nonempty": bool(te.sum() > 0),
            "prediction_rows": int(pack["n_groups"]),
            "leaky_suspect_prior_not_used": {
                "t2_r2": 0.466,
                "t2_mag_r": 0.757,
                "label": "pre-protocol-fix leaky-suspect; clean re-run held before this build",
            },
        }

        write_predictions(pred_csv, pack, best_result)
        write_json(metrics_json, metrics)
        if not args.no_graph:
            render_graph(png, pdf, pack, best_result, metrics)
        audit = {
            "disk_before_artifact_write": disk_audit,
            "boundary": {
                "python_opened_trajectory_h5": False,
                "python_spatial_or_protein_graph": False,
                "python_physics_recompute": False,
                "basis_conversion": "change_of_basis.get_C/lib_to_e3nn only",
                "extractor_modified_by_this_script": False,
                "target": "5-component T2 plus companion T0; T2 not scalar-collapsed",
            },
            "artifacts": {
                "predictions_csv": str(pred_csv),
                "metrics_json": str(metrics_json),
                "advisor_graph_png": str(png) if not args.no_graph else None,
                "advisor_graph_pdf": str(pdf) if not args.no_graph else None,
            },
        }
        write_json(audit_json, audit)
        print(f"wrote predictions: {pred_csv}")
        print(f"wrote metrics: {metrics_json}")
        if not args.no_graph:
            print(f"wrote graph: {png}")
            print(f"wrote graph: {pdf}")


if __name__ == "__main__":
    main()
