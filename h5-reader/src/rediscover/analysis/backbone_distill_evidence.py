#!/usr/bin/env python3
"""backbone_distill_evidence -- read out the broad-backbone e3nn T2 model.

This is an evidence collector for REDISCOVERY_MAP item 1. It stays on the
Python-consumer side of the boundary:

* inputs are the emitted broad_backbone CSV/NPY substrate;
* angular features are the existing e3nn features from equiv_t2_backbone_e3nn.py;
* target T2 is mapped with the frozen change_of_basis.get_C() path;
* the script reads learned radial gates from the model and fits requested
  radial-only comparators. It does not rebuild an emitted kernel, field, or T2
  target from raw coordinates.

Outputs:
  <evidence_dir>/frame_split_metrics.csv
  <evidence_dir>/readout_samples.csv.gz
  <evidence_dir>/radial_fit_summary.csv
  <evidence_dir>/atom_split_validation.csv
  <evidence_dir>/models/backbone_<stratum>.pt
  <evidence_dir>/figures/radial_<kind>.png
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

import change_of_basis as cob
from equiv_t2_backbone_e3nn import (
    KIND_FEATURES,
    STRATA_ORDER,
    build_pack,
    detect_axis_columns,
    load,
    r2,
    stratum_of,
    _new_model,
)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default="/tmp/rdc-broad-backbone-axes")
    ap.add_argument("--evidence-dir", default="/tmp/rdc-backbone-law-evidence")
    ap.add_argument("--epochs", type=int, default=4000)
    ap.add_argument("--validation-epochs", type=int, default=500)
    ap.add_argument("--validation-splits", type=int, default=2)
    ap.add_argument("--lr", type=float, default=3e-3)
    ap.add_argument("--sample-per-kind", type=int, default=50000)
    ap.add_argument("--grid-size", type=int, default=80)
    ap.add_argument("--skip-validation", action="store_true")
    ap.add_argument("--strata", nargs="*", default=STRATA_ORDER)
    return ap.parse_args()


def as_np(x):
    return x.detach().cpu().numpy()


def finite_r2(y, yhat):
    ok = np.isfinite(y) & np.isfinite(yhat)
    if ok.sum() < 3:
        return float("nan")
    yy = y[ok]
    pp = yhat[ok]
    den = ((yy - yy.mean()) ** 2).sum()
    if den <= 0:
        return float("nan")
    return float(1.0 - ((yy - pp) ** 2).sum() / den)


def ols_fit(X, y):
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if ok.sum() < X.shape[1] + 3:
        return None
    Xo = X[ok]
    yo = y[ok]
    beta, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
    pred = Xo @ beta
    den = ((yo - yo.mean()) ** 2).sum()
    r2v = float("nan") if den <= 0 else float(1.0 - ((yo - pred) ** 2).sum() / den)
    return beta, r2v, int(ok.sum())


def train_model(pack, dev, epochs, lr, train_mask=None, verbose_prefix=""):
    model = _new_model(pack, "exact", dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-7)
    target = pack["target"]
    if train_mask is None:
        train_mask_t = pack["g_tr"]
    else:
        train_mask_t = torch.tensor(train_mask, dtype=torch.bool, device=dev)
    if int(train_mask_t.sum().item()) < 3:
        return model, float("nan"), float("nan")
    for ep in range(epochs):
        model.train()
        opt.zero_grad()
        pred = model(pack["per_kind"], pack["group_atom"], center_mask=train_mask_t)
        loss = ((pred[train_mask_t] - target[train_mask_t]) ** 2).mean()
        loss.backward()
        opt.step()
        if verbose_prefix and (ep % 1000 == 0 or ep == epochs - 1):
            with torch.no_grad():
                pr = model(pack["per_kind"], pack["group_atom"], center_mask=train_mask_t)
            te = pack["g_te"]
            print(
                f"{verbose_prefix} ep={ep:4d} loss={loss.item():.4e} "
                f"train_R2={r2(pr[train_mask_t], target[train_mask_t]):+.3f} "
                f"frame_R2={r2(pr[te], target[te]):+.3f}",
                flush=True,
            )
    model.eval()
    with torch.no_grad():
        pred = model(pack["per_kind"], pack["group_atom"], center_mask=train_mask_t)
    te = pack["g_te"]
    pred_np = as_np(pred)
    tgt_np = as_np(target)
    te_np = as_np(te).astype(bool)
    r2_frame = r2(torch.tensor(pred_np[te_np]), torch.tensor(tgt_np[te_np]))
    mag_r = float(
        np.corrcoef(np.linalg.norm(pred_np[te_np], axis=1), np.linalg.norm(tgt_np[te_np], axis=1))[0, 1]
    )
    return model, r2_frame, mag_r


def prepare_strata(agg):
    out = agg.copy()
    out["__stratum"] = [stratum_of(fv, nm) for fv, nm in zip(out.frame_variant, out.atom_name)]
    return out


def collect_raw_sources(src, agg_stratum):
    keep_agg = (agg_stratum.dft_present == 1) & (agg_stratum.dft_local_frame_valid == 1)
    a = agg_stratum.loc[keep_agg].copy().reset_index(drop=True)
    a["__gid"] = np.arange(len(a))
    gid_of_rowid = dict(zip(a.row_id.to_numpy(), a.__gid.to_numpy()))
    s = src[src.row_id.isin(set(a.row_id))].copy()
    s = s[~((s.mechanism == "ring") & (s.is_self_or_bonded == 1))]
    s["__gid"] = s.row_id.map(gid_of_rowid).astype("int64")
    return {k: s[s.mechanism == k].copy() for k in KIND_FEATURES}


def kind_gate_names(n_gate):
    if n_gate == 3:
        return ["disp", "axis", "cross"]
    return ["disp"]


def chosen_gate_name(kind, n_gate):
    if kind in ("ring", "bond") and n_gate == 3:
        return "axis"
    return "disp"


def readout_model(stratum, model, pack, raw_by_kind, sample_per_kind, rng):
    rows = []
    summary = []
    for kind, raw in raw_by_kind.items():
        pk = pack["per_kind"][kind]
        n_src = int(pk["src_group"].numel())
        if n_src == 0:
            continue
        if len(raw) != n_src:
            raise RuntimeError(f"{stratum}/{kind}: raw rows {len(raw)} != tensor rows {n_src}")
        ch = model.channels[kind]
        ch.eval()
        with torch.no_grad():
            gates = as_np(ch.R(pk["featn"]))
            contrib = as_np(ch(pk))
        names = kind_gate_names(gates.shape[1])
        chosen = chosen_gate_name(kind, gates.shape[1])
        chosen_idx = names.index(chosen)
        gate_norm = np.linalg.norm(gates, axis=1)
        contrib_norm = np.linalg.norm(contrib, axis=1)

        take_n = min(sample_per_kind, n_src)
        idx = rng.choice(n_src, size=take_n, replace=False)
        out = pd.DataFrame(
            {
                "stratum": stratum,
                "kind": kind,
                "row_id": raw.row_id.to_numpy()[idx],
                "atom_index": raw.atom_index.to_numpy()[idx],
                "__gid": raw.__gid.to_numpy()[idx],
                "r": raw.r.to_numpy(float)[idx],
                "cos_theta": raw.cos_theta.to_numpy(float)[idx],
                "ring_intensity": raw.ring_intensity.to_numpy(float)[idx],
                "bond_category": raw.bond_category.to_numpy(int)[idx],
                "source_q_e": raw.source_q_e.to_numpy(float)[idx],
                "gate_chosen_name": chosen,
                "gate_chosen": gates[idx, chosen_idx],
                "gate_norm": gate_norm[idx],
                "contrib_norm": contrib_norm[idx],
            }
        )
        for gi, name in enumerate(names):
            out[f"gate_{name}"] = gates[idx, gi]
        for ci in range(5):
            out[f"contrib_{ci}"] = contrib[idx, ci]
        rows.append(out)
        summary.append(
            {
                "stratum": stratum,
                "kind": kind,
                "n_sources": n_src,
                "sampled": take_n,
                "chosen_gate": chosen,
                "gate_columns": ",".join(names),
                "gate_chosen_mean": float(np.mean(gates[:, chosen_idx])),
                "gate_chosen_std": float(np.std(gates[:, chosen_idx])),
                "gate_norm_mean": float(np.mean(gate_norm)),
                "contrib_norm_mean": float(np.mean(contrib_norm)),
            }
        )
    return rows, summary


def add_radial_fit(records, df):
    if len(df) < 20:
        return
    stratum = str(df.stratum.iloc[0])
    kind = str(df.kind.iloc[0])
    y = df.gate_chosen.to_numpy(float)
    r = df.r.to_numpy(float)
    eps = 1e-12

    if kind == "ring":
        inten = df.ring_intensity.to_numpy(float)
        ok = (r > 0) & (np.abs(inten) > 1e-8) & (np.abs(y) > eps)
        yy = np.log(np.abs(y[ok] / inten[ok]))
        X = np.column_stack([np.ones(ok.sum()), np.log(r[ok])])
        fit = ols_fit(X, yy)
        if fit:
            beta, r2_log, n = fit
            records.append(
                {
                    "stratum": stratum,
                    "kind": kind,
                    "fit": "log_abs_gate_over_intensity_vs_log_r",
                    "power": float(beta[1]),
                    "r2": r2_log,
                    "n": n,
                    "note": "chosen learned gate; radial-only log fit",
                }
            )
        x = inten * np.power(r, -3.0)
        X = np.column_stack([np.ones(len(x)), x])
        fit = ols_fit(X, y)
        if fit:
            beta, r2_lin, n = fit
            records.append(
                {
                    "stratum": stratum,
                    "kind": kind,
                    "fit": "gate_vs_intensity_r^-3",
                    "power": -3.0,
                    "r2": r2_lin,
                    "n": n,
                    "note": "requested radial-only analytic comparator",
                }
            )
    elif kind == "bond":
        cat = df.bond_category.to_numpy(int)
        ok = (r > 0) & (cat >= 0) & (np.abs(y) > eps)
        cats = sorted(set(cat[ok]))
        if cats:
            onehot = np.column_stack([(cat[ok] == c).astype(float) for c in cats])
            X = np.column_stack([np.ones(ok.sum()), np.log(r[ok]), onehot])
            fit = ols_fit(X, np.log(np.abs(y[ok])))
            if fit:
                beta, r2_log, n = fit
                records.append(
                    {
                        "stratum": stratum,
                        "kind": kind,
                        "fit": "log_abs_gate_vs_log_r_plus_category_offsets",
                        "power": float(beta[1]),
                        "r2": r2_log,
                        "n": n,
                        "note": f"categories={cats}; chosen learned gate",
                    }
                )
            cols = [np.where(cat == c, np.power(r, -3.0), 0.0) for c in cats]
            X = np.column_stack([np.ones(len(r))] + cols)
            fit = ols_fit(X, y)
            if fit:
                beta, r2_lin, n = fit
                records.append(
                    {
                        "stratum": stratum,
                        "kind": kind,
                        "fit": "gate_vs_category_r^-3",
                        "power": -3.0,
                        "r2": r2_lin,
                        "n": n,
                        "note": f"categories={cats}; radial-only analytic comparator",
                    }
                )
    elif kind == "charge":
        q = df.source_q_e.to_numpy(float)
        ok = (r > 0) & (np.abs(q) > 1e-5) & (np.abs(y) > eps)
        yy = np.log(np.abs(y[ok] / q[ok]))
        X = np.column_stack([np.ones(ok.sum()), np.log(r[ok])])
        fit = ols_fit(X, yy)
        if fit:
            beta, r2_log, n = fit
            records.append(
                {
                    "stratum": stratum,
                    "kind": kind,
                    "fit": "log_abs_gate_over_q_vs_log_r",
                    "power": float(beta[1]),
                    "r2": r2_log,
                    "n": n,
                    "note": "chosen learned gate; radial-only log fit",
                }
            )
        x = np.multiply(q, np.power(r, -3.0))
        X = np.column_stack([np.ones(len(x)), x])
        fit = ols_fit(X, y)
        if fit:
            beta, r2_lin, n = fit
            records.append(
                {
                    "stratum": stratum,
                    "kind": kind,
                    "fit": "gate_vs_q_r^-3",
                    "power": -3.0,
                    "r2": r2_lin,
                    "n": n,
                    "note": "T2/EFG radial-only analytic comparator",
                }
            )


def validation_pack_for_kind(pack, kind):
    out = dict(pack)
    out["per_kind"] = {kind: pack["per_kind"][kind]}
    out["kind_nfeat"] = {kind: pack["kind_nfeat"][kind]}
    out["kind_axis"] = {kind: pack["kind_axis"][kind]}
    return out


def atom_split_validation(stratum, kind, pack, dev, epochs, lr, splits, rng):
    atoms = np.unique(pack["group_atom_idx"])
    if len(atoms) < 3:
        return []
    out = []
    ppack = validation_pack_for_kind(pack, kind)
    ga = pack["group_atom_idx"]
    for split in range(splits):
        perm = rng.permutation(atoms)
        n_train_atoms = max(1, len(perm) // 2)
        train_atoms = set(perm[:n_train_atoms])
        test_atoms = set(perm[n_train_atoms:])
        train_mask = np.array([a in train_atoms for a in ga], dtype=bool)
        test_mask = np.array([a in test_atoms for a in ga], dtype=bool)
        model, _, _ = train_model(ppack, dev, epochs, lr, train_mask=train_mask, verbose_prefix="")
        model.eval()
        with torch.no_grad():
            train_mask_t = torch.tensor(train_mask, dtype=torch.bool, device=dev)
            pred = as_np(model(ppack["per_kind"], ppack["group_atom"],
                               center_mask=train_mask_t))
        tgt = as_np(ppack["target"])
        if test_mask.sum() < 3:
            continue
        mag_p = np.linalg.norm(pred[test_mask], axis=1)
        mag_t = np.linalg.norm(tgt[test_mask], axis=1)
        mag_r = float(np.corrcoef(mag_p, mag_t)[0, 1]) if len(mag_p) > 2 else float("nan")
        out.append(
            {
                "stratum": stratum,
                "kind": kind,
                "split": split,
                "train_atoms": len(train_atoms),
                "heldout_atoms": len(test_atoms),
                "heldout_groups": int(test_mask.sum()),
                "epochs": epochs,
                "t2_r2": r2(torch.tensor(pred[test_mask]), torch.tensor(tgt[test_mask])),
                "t2_mag_r": mag_r,
                "note": "source-kind-only e3nn model trained on atom split; target is independent DFT T2",
            }
        )
    return out


def make_plots(evidence_dir, readout):
    fig_dir = evidence_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for kind in ["ring", "bond", "charge"]:
        d = readout[readout.kind == kind]
        if d.empty:
            continue
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for stratum, g in d.groupby("stratum"):
            if len(g) > 6000:
                g = g.sample(6000, random_state=0)
            y = np.abs(g.gate_chosen.to_numpy(float))
            ok = (g.r.to_numpy(float) > 0) & (y > 1e-12)
            ax.scatter(g.r.to_numpy(float)[ok], y[ok], s=3, alpha=0.12, label=stratum)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("r (A)")
        ax.set_ylabel("|chosen learned radial gate|")
        ax.set_title(f"{kind}: learned radial readout sample")
        ax.legend(ncol=2, fontsize=7)
        fig.tight_layout()
        fig.savefig(fig_dir / f"radial_{kind}.png", dpi=130)
        plt.close(fig)


def main():
    args = parse_args()
    evidence_dir = Path(args.evidence_dir)
    evidence_dir.mkdir(parents=True, exist_ok=True)
    (evidence_dir / "models").mkdir(exist_ok=True)
    print(f"evidence_dir={evidence_dir}", flush=True)

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={dev} out_dir={args.out_dir}", flush=True)
    C = cob.get_C()
    print(f"change_of_basis max |C.T C - I|={np.abs(C.T @ C - np.eye(5)).max():.3e}", flush=True)

    src, agg, target_lib_all = load(args.out_dir)
    axis_present = detect_axis_columns(src)
    print(
        f"axis columns present: ring={axis_present['ring']} bond={axis_present['bond']}",
        flush=True,
    )
    agg = prepare_strata(agg)

    rng = np.random.default_rng(0)
    metric_rows = []
    readout_parts = []
    readout_summaries = []
    validation_rows = []

    for stratum in args.strata:
        a_orig = agg[agg["__stratum"] == stratum]
        if a_orig.empty:
            continue
        tgt_s = target_lib_all[a_orig.index.to_numpy()]
        a_s = a_orig.reset_index(drop=True)
        pack = build_pack(src, a_s, tgt_s, axis_present, True, dev)
        if pack is None or pack["n_groups"] < 4 or pack["n_atoms"] < 2:
            continue
        raw_by_kind = collect_raw_sources(src, a_s)
        print(
            f"\n[{stratum}] groups={pack['n_groups']} atoms={pack['n_atoms']} "
            f"sources ring/bond/charge="
            f"{pack['n_src'].get('ring', 0)}/{pack['n_src'].get('bond', 0)}/{pack['n_src'].get('charge', 0)}",
            flush=True,
        )
        model, r2_frame, mag_r = train_model(
            pack,
            dev,
            args.epochs,
            args.lr,
            verbose_prefix=f"[{stratum}]",
        )
        metric_rows.append(
            {
                "stratum": stratum,
                "groups": pack["n_groups"],
                "atoms": pack["n_atoms"],
                "test_groups": int(as_np(pack["g_te"]).sum()),
                "epochs": args.epochs,
                "frame_t2_r2": r2_frame,
                "frame_t2_mag_r": mag_r,
                "source_ring": pack["n_src"].get("ring", 0),
                "source_bond": pack["n_src"].get("bond", 0),
                "source_charge": pack["n_src"].get("charge", 0),
            }
        )
        torch.save(
            {
                "stratum": stratum,
                "state_dict": model.state_dict(),
                "kind_nfeat": pack["kind_nfeat"],
                "kind_axis": pack["kind_axis"],
                "epochs": args.epochs,
                "frame_t2_r2": r2_frame,
                "frame_t2_mag_r": mag_r,
            },
            evidence_dir / "models" / f"backbone_{stratum}.pt",
        )
        parts, summaries = readout_model(
            stratum, model, pack, raw_by_kind, args.sample_per_kind, rng
        )
        readout_parts.extend(parts)
        readout_summaries.extend(summaries)

        if not args.skip_validation:
            for kind in ["ring", "bond", "charge"]:
                if pack["n_src"].get(kind, 0) == 0:
                    continue
                print(
                    f"[{stratum}/{kind}] atom-split validation "
                    f"splits={args.validation_splits} epochs={args.validation_epochs}",
                    flush=True,
                )
                validation_rows.extend(
                    atom_split_validation(
                        stratum,
                        kind,
                        pack,
                        dev,
                        args.validation_epochs,
                        args.lr,
                        args.validation_splits,
                        rng,
                    )
                )

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(evidence_dir / "frame_split_metrics.csv", index=False)
    pd.DataFrame(readout_summaries).to_csv(evidence_dir / "readout_summary.csv", index=False)
    if readout_parts:
        readout = pd.concat(readout_parts, ignore_index=True)
    else:
        readout = pd.DataFrame()
    readout.to_csv(evidence_dir / "readout_samples.csv.gz", index=False)

    radial_records = []
    if not readout.empty:
        for (stratum, kind), d in readout.groupby(["stratum", "kind"]):
            add_radial_fit(radial_records, d)
        make_plots(evidence_dir, readout)
    pd.DataFrame(radial_records).to_csv(evidence_dir / "radial_fit_summary.csv", index=False)
    pd.DataFrame(validation_rows).to_csv(evidence_dir / "atom_split_validation.csv", index=False)

    manifest = {
        "out_dir": args.out_dir,
        "evidence_dir": str(evidence_dir),
        "epochs": args.epochs,
        "validation_epochs": args.validation_epochs,
        "validation_splits": args.validation_splits,
        "sample_per_kind": args.sample_per_kind,
        "axis_present": axis_present,
        "device": dev,
    }
    (evidence_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print("\nWrote evidence files:", flush=True)
    for name in [
        "frame_split_metrics.csv",
        "readout_summary.csv",
        "readout_samples.csv.gz",
        "radial_fit_summary.csv",
        "atom_split_validation.csv",
        "manifest.json",
    ]:
        print(f"  {evidence_dir / name}", flush=True)


if __name__ == "__main__":
    main()
