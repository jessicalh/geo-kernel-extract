#!/usr/bin/env python3
"""Standalone APBS EFG T2 -> DFT T2 fitter for the six backbone strata.

Inputs are only the producer-emitted EFG sidecars:
  efg_aggregated.csv, efg_feature_T2.npy, efg_target_T2.npy

The model is the Schur-scalar law for a 2e -> 2e equivariant map:
  pred_T2(atom, frame) = g(|EFG|) * EFG_T2(atom, frame)

Both feature and target are mapped through the pinned library->e3nn 2e basis
constant from change_of_basis.get_C(). No tensor construction or protein data
loading happens here.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

import change_of_basis as cob


torch.manual_seed(0)
np.random.seed(0)

FV_HN = {1, 2}          # HN_Standard, HN_NTerminus
FV_N = {4, 5}           # BackboneN, BackboneN_NTerminus
FV_CA = {6}             # BackboneCA
FV_C = {7}              # BackboneCarbonylC
FV_O = {8}              # BackboneCarbonylO
FV_HA = {9}             # BackboneHA

STRATA_ORDER = ["N", "CA", "C", "O", "HN", "HA"]

TEST_FRAME_FRACTION = 0.20
FRAME_SPLIT_SEED = 0
MIN_FRAME_SPLIT_FRAMES = 5
THIN_ATOM_WARN = 10
DEFAULT_SPLIT = "blocked"
PURGE_NEIGHBOUR_FRAMES = 1


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?",
                    default=os.environ.get("REDISCOVER_OUT", "/tmp/rdc-efg"))
    ap.add_argument("--epochs", type=int, default=4000)
    ap.add_argument("--lr", type=float, default=3e-3)
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--loao", action="store_true",
                    help="also run leave-atoms-out retraining. The frame split "
                         "remains the reported gate.")
    ap.add_argument("--split", choices=["blocked", "random"], default=DEFAULT_SPLIT,
                    help="frame split for the emitted EFG within gate. Default "
                         "is a contiguous temporal block with purged neighbours.")
    ap.add_argument("--purge-frames", type=int, default=PURGE_NEIGHBOUR_FRAMES,
                    help="for --split blocked, drop this many neighbouring frame "
                         "indices from train around the held-out block.")
    return ap.parse_args()


def stratum_of(frame_variant, atom_name):
    """Six-class backbone mapping from emitted frame_variant + atom_name.

    The frame_variant sets mirror equiv_t2_backbone_e3nn.py. That exemplar splits
    glycine HA2/HA3 from HA for the broad eight-class fit; this EFG brief asks for
    the six backbone classes, so all FV_HA rows are reported as HA.
    """
    del atom_name
    fv = int(frame_variant)
    if fv in FV_N:
        return "N"
    if fv in FV_CA:
        return "CA"
    if fv in FV_C:
        return "C"
    if fv in FV_O:
        return "O"
    if fv in FV_HN:
        return "HN"
    if fv in FV_HA:
        return "HA"
    return None


def r2(pred, target):
    denom = ((target - target.mean(axis=0, keepdims=True)) ** 2).sum()
    if denom <= 0:
        return float("nan")
    return float(1.0 - ((target - pred) ** 2).sum() / denom)


def corrcoef(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 2:
        return float("nan")
    a = a[ok]
    b = b[ok]
    if np.std(a) <= 0 or np.std(b) <= 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def require(path, what):
    if not os.path.exists(path):
        sys.exit(f"FATAL: required emitted {what} not found:\n  {path}\n"
                 "Run the extractor with --case efg into this output directory.")


def load(out_dir):
    agg_csv = os.path.join(out_dir, "efg_aggregated.csv")
    feature_npy = os.path.join(out_dir, "efg_feature_T2.npy")
    target_npy = os.path.join(out_dir, "efg_target_T2.npy")
    require(agg_csv, "aggregated CSV")
    require(feature_npy, "feature T2 NPY")
    require(target_npy, "target T2 NPY")

    agg = pd.read_csv(agg_csv)
    feature_lib = np.load(feature_npy)
    target_lib = np.load(target_npy)
    if feature_lib.shape != target_lib.shape:
        sys.exit(f"FATAL: feature shape {feature_lib.shape} != target shape "
                 f"{target_lib.shape} (sidecars out of sync).")
    if len(agg) != feature_lib.shape[0]:
        sys.exit(f"FATAL: CSV rows {len(agg)} != NPY rows {feature_lib.shape[0]} "
                 "(sidecars out of sync).")
    if feature_lib.ndim != 2 or feature_lib.shape[1] != 5:
        sys.exit(f"FATAL: efg_feature_T2.npy has shape {feature_lib.shape}, "
                 "expected (rows, 5).")
    return agg, feature_lib, target_lib


class ScalarGate(nn.Module):
    def __init__(self, hidden, init_gate):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden), nn.SiLU(),
            nn.Linear(hidden, 1),
        )
        final = self.net[-1]
        nn.init.zeros_(final.weight)
        nn.init.constant_(final.bias, float(init_gate))

    def forward(self, mag_norm):
        return self.net(mag_norm)


class EfgLaw(nn.Module):
    def __init__(self, n_atoms, hidden, init_gate):
        super().__init__()
        self.n_atoms = n_atoms
        self.gate = ScalarGate(hidden, init_gate)

    def forward(self, feature, mag_norm, group_atom, center_mask=None):
        pred = self.gate(mag_norm) * feature
        if center_mask is None:
            center_mask = torch.ones(group_atom.shape, dtype=torch.bool, device=group_atom.device)
        mean = torch.zeros(self.n_atoms, 5, dtype=pred.dtype, device=pred.device)
        cnt = torch.zeros(self.n_atoms, 1, dtype=pred.dtype, device=pred.device)
        center_atom = group_atom[center_mask]
        center_pred = pred[center_mask]
        mean.index_add_(0, center_atom, center_pred)
        cnt.index_add_(0, center_atom, torch.ones_like(cnt[center_atom]))
        return pred - (mean / cnt.clamp_min(1.0))[group_atom]


def centred_by_train_atom(x, group_atom_idx, train_mask):
    out = np.full_like(np.asarray(x, dtype=float), np.nan, dtype=float)
    for atom in np.unique(group_atom_idx):
        m = group_atom_idx == atom
        train_atom = m & train_mask
        if train_atom.sum() == 0:
            continue
        out[m] = x[m] - x[train_atom].mean(axis=0, keepdims=True)
    return out


def make_split_masks(row_frames, strategy=DEFAULT_SPLIT, seed=FRAME_SPLIT_SEED,
                     purge_frames=PURGE_NEIGHBOUR_FRAMES):
    frames = np.sort(np.unique(row_frames))
    empty = np.zeros(len(row_frames), dtype=bool)
    if len(frames) < MIN_FRAME_SPLIT_FRAMES:
        return empty.copy(), empty.copy(), {
            "split_strategy": strategy,
            "test_frames": 0,
            "purged_train_frames": 0,
            "cross_split_lag1_pairs": 0,
        }
    n_test = max(1, int(TEST_FRAME_FRACTION * len(frames)))
    rng = np.random.default_rng(seed)

    if strategy == "random":
        test_frames = set(rng.choice(frames, n_test, replace=False))
        train_frames = set(frames) - test_frames
        purged = set()
    elif strategy == "blocked":
        start = int(rng.integers(0, len(frames) - n_test + 1))
        stop = start + n_test
        test_frames = set(frames[start:stop])
        purge_lo = max(0, start - max(0, purge_frames))
        purge_hi = min(len(frames), stop + max(0, purge_frames))
        purged = set(frames[purge_lo:start]) | set(frames[stop:purge_hi])
        train_frames = set(frames) - test_frames - purged
    else:
        raise ValueError(f"unknown split strategy {strategy!r}")

    train = np.array([f in train_frames for f in row_frames], dtype=bool)
    test = np.array([f in test_frames for f in row_frames], dtype=bool)

    frame_to_split = {f: ("test" if f in test_frames else "train" if f in train_frames else "purged")
                      for f in frames}
    cross = 0
    for a, b in zip(frames[:-1], frames[1:]):
        sa = frame_to_split[a]
        sb = frame_to_split[b]
        if {sa, sb} == {"train", "test"}:
            cross += 1
    return train, test, {
        "split_strategy": strategy,
        "test_frames": int(len(test_frames)),
        "purged_train_frames": int(len(purged)),
        "cross_split_lag1_pairs": int(cross),
    }


def build_pack(agg_s, feature_lib_s, target_lib_s, C, dev,
               split_strategy=DEFAULT_SPLIT, split_seed=FRAME_SPLIT_SEED,
               purge_frames=PURGE_NEIGHBOUR_FRAMES):
    present = (
        (agg_s["dft_present"].to_numpy() == 1)
        & (agg_s["apbs_efg_present"].to_numpy() == 1)
        & np.isfinite(agg_s["efg_T2_magnitude"].to_numpy(float))
        & np.isfinite(feature_lib_s).all(axis=1)
        & np.isfinite(target_lib_s).all(axis=1)
    )
    a = agg_s.loc[present].copy().reset_index(drop=True)
    if not len(a):
        return None
    feature_lib = feature_lib_s[present]
    target_lib = target_lib_s[present]

    feature = np.asarray(feature_lib, dtype=float) @ C.T
    target = np.asarray(target_lib, dtype=float) @ C.T

    group_atom_idx, atom_labels = pd.factorize(a["atom_index"].to_numpy())
    n_atoms = int(group_atom_idx.max()) + 1
    g_tr, g_te, split_info = make_split_masks(a["h5_row"].to_numpy(),
                                              strategy=split_strategy,
                                              seed=split_seed,
                                              purge_frames=purge_frames)
    target_dm = centred_by_train_atom(target, group_atom_idx, g_tr)
    feature_dm_constant_gate = centred_by_train_atom(feature, group_atom_idx, g_tr)

    mag = a["efg_T2_magnitude"].to_numpy(float).reshape(-1, 1)
    norm_source = mag[g_tr] if g_tr.any() else mag
    mag_mean = norm_source.mean(axis=0)
    mag_std = norm_source.std(axis=0)
    mag_norm = (mag - mag_mean) / np.where(mag_std > 1e-8, mag_std, 1.0)

    denom = (feature_dm_constant_gate[g_tr] ** 2).sum() if g_tr.any() else 0.0
    if denom > 0:
        init_gate = float((target_dm[g_tr] * feature_dm_constant_gate[g_tr]).sum() / denom)
    else:
        init_gate = 0.0

    t = lambda x, dt=torch.float32: torch.tensor(x, dtype=dt, device=dev)
    return {
        "feature": t(feature),
        "target": t(target_dm),
        "mag_norm": t(mag_norm),
        "mag": mag.reshape(-1),
        "mag_mean": float(mag_mean[0]),
        "mag_std": float(mag_std[0] if mag_std[0] > 1e-8 else 1.0),
        "group_atom": t(group_atom_idx, torch.long),
        "group_atom_idx": group_atom_idx,
        "atom_labels": atom_labels,
        "h5_row": a["h5_row"].to_numpy(),
        "g_te": t(g_te, torch.bool),
        "g_tr": t(g_tr, torch.bool),
        "n_groups": len(a),
        "n_atoms": n_atoms,
        "n_frames": int(a["h5_row"].nunique()),
        "test_groups": int(g_te.sum()),
        "train_groups": int(g_tr.sum()),
        **split_info,
        "init_gate": init_gate,
        "feature_scale_median": float(np.median(a["efg_T2_magnitude"].to_numpy(float))),
    }


def new_model(pack, args, dev):
    return EfgLaw(pack["n_atoms"], args.hidden, pack["init_gate"]).to(dev)


def train_model(pack, args, dev, train_mask=None, epochs=None):
    model = new_model(pack, args, dev)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-7)
    if train_mask is None:
        train_mask = pack["g_tr"]
    if epochs is None:
        epochs = args.epochs
    if int(train_mask.sum().item()) < 2:
        return None

    target = pack["target"]
    scale = target[train_mask].std().clamp_min(torch.tensor(1.0, device=dev))
    for _ in range(epochs):
        model.train()
        opt.zero_grad()
        pred = model(pack["feature"], pack["mag_norm"], pack["group_atom"],
                     center_mask=train_mask)
        loss = (((pred[train_mask] - target[train_mask]) / scale) ** 2).mean()
        loss.backward()
        opt.step()
    return model


def score_model(model, pack):
    if model is None:
        return float("nan"), float("nan")
    with torch.no_grad():
        pred = model(pack["feature"], pack["mag_norm"], pack["group_atom"],
                     center_mask=pack["g_tr"]).cpu().numpy()
    target = pack["target"].cpu().numpy()
    te = pack["g_te"].cpu().numpy()
    if te.sum() < 2:
        return float("nan"), float("nan")
    r2_t2 = r2(pred[te], target[te])
    mag_r = corrcoef(np.linalg.norm(pred[te], axis=1),
                     np.linalg.norm(target[te], axis=1))
    return r2_t2, mag_r


def leave_atoms_out(pack, args, dev):
    atoms = np.unique(pack["group_atom_idx"])
    if len(atoms) < 3:
        return float("nan")

    pred = np.full(pack["target"].shape, np.nan, dtype=float)
    ga = pack["group_atom_idx"]
    epochs = max(800, args.epochs // 4)
    for atom in atoms:
        train_mask = torch.tensor(ga != atom, dtype=torch.bool, device=dev)
        model = train_model(pack, args, dev, train_mask=train_mask, epochs=epochs)
        if model is None:
            continue
        with torch.no_grad():
            p = model(pack["feature"], pack["mag_norm"], pack["group_atom"],
                      center_mask=train_mask).cpu().numpy()
        pred[ga == atom] = p[ga == atom]

    target = pack["target"].cpu().numpy()
    ok = np.isfinite(pred).all(axis=1)
    if ok.sum() < 2:
        return float("nan")
    return r2(pred[ok], target[ok])


def fmt(x):
    return "nan" if not np.isfinite(x) else f"{x:+.3f}"


def report_stratum(name, pack, args, dev):
    thin = pack["n_atoms"] < THIN_ATOM_WARN
    thin_txt = " THIN" if thin else ""
    print(f"\n-- stratum {name}: rows={pack['n_groups']}  "
          f"effective_N={pack['n_atoms']}  frames={pack['n_frames']}  "
          f"train_rows={pack['train_groups']}  test_rows={pack['test_groups']}  "
          f"split={pack['split_strategy']}  cross_lag1_pairs={pack['cross_split_lag1_pairs']}"
          f"{thin_txt}  "
          f"median_|EFG|={pack['feature_scale_median']:.3g}",
          flush=True)
    if pack["n_atoms"] < 2 or pack["n_groups"] < 4:
        print("   insufficient effective N to fit honestly", flush=True)
        return None
    model = train_model(pack, args, dev)
    r2_t2, mag_r = score_model(model, pack)
    loao = leave_atoms_out(pack, args, dev) if args.loao else float("nan")
    loao_txt = fmt(loao) if args.loao else "skipped"
    print(f"   frame-split T2 R2={fmt(r2_t2)}  |T2| r={fmt(mag_r)}  "
          f"leave-atoms-out R2={loao_txt}",
          flush=True)
    return {
        "r2": r2_t2,
        "mag_r": mag_r,
        "n_atoms": pack["n_atoms"],
        "n_groups": pack["n_groups"],
        "thin": thin,
        "loao": loao,
    }


def main():
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    C = cob.get_C()
    print(f"device={dev}  out_dir={args.out_dir}  epochs={args.epochs}  "
          f"lr={args.lr}  hidden={args.hidden}  loao={args.loao}  "
          f"split={args.split}  purge_frames={args.purge_frames}",
          flush=True)
    print(f"change_of_basis.get_C() |C^T C - I|max="
          f"{np.abs(C.T @ C - np.eye(5)).max():.2e}", flush=True)

    agg, feature_lib, target_lib = load(args.out_dir)
    print(f"emitted rows={len(agg)}  feature_shape={feature_lib.shape}  "
          f"target_shape={target_lib.shape}", flush=True)

    agg = agg.copy()
    agg["__stratum"] = [stratum_of(fv, nm)
                        for fv, nm in zip(agg["frame_variant"], agg["atom_name"])]
    unrecognised = int(agg["__stratum"].isna().sum())
    if unrecognised:
        print(f"ignored non-backbone/unmapped rows={unrecognised}", flush=True)

    summary = {}
    for name in STRATA_ORDER:
        idx = np.flatnonzero(agg["__stratum"].to_numpy() == name)
        if len(idx) == 0:
            print(f"\n-- stratum {name}: absent", flush=True)
            summary[name] = None
            continue
        pack = build_pack(agg.iloc[idx].reset_index(drop=True),
                          feature_lib[idx], target_lib[idx], C, dev,
                          split_strategy=args.split, split_seed=args.seed,
                          purge_frames=args.purge_frames)
        if pack is None:
            print(f"\n-- stratum {name}: no rows with both DFT and APBS EFG present",
                  flush=True)
            summary[name] = None
            continue
        summary[name] = report_stratum(name, pack, args, dev)

    print("\n=== EFG scalar-gate summary ===", flush=True)
    print("stratum,T2_R2_frame_split,absT2_r_frame_split,effective_N,rows,flag",
          flush=True)
    for name in STRATA_ORDER:
        res = summary.get(name)
        if not res:
            print(f"{name},nan,nan,0,0,INSUFFICIENT_N", flush=True)
            continue
        flag = "THIN" if res["thin"] else ""
        print(f"{name},{fmt(res['r2'])},{fmt(res['mag_r'])},"
              f"{res['n_atoms']},{res['n_groups']},{flag}",
              flush=True)


if __name__ == "__main__":
    main()
