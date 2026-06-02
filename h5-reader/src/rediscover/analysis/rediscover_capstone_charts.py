#!/usr/bin/env python3
"""Reusable rediscover capstone chart set.

No notebooks. Reads generated CSV tables only and writes publication-oriented
matplotlib figures under a caller-provided output directory. To rerun on the
future 750-DFT substrate, regenerate the same CSVs and pass the new paths.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STRATA_ORDER = ["N", "CA", "C", "O", "HN", "HA"]
KERNEL_ORDER = ["total", "ring", "bond", "charge"]
COLORS = {
    "between": "#3b7ea1",
    "within": "#e7a246",
    "lab": "#9b4d64",
    "local": "#3f8f73",
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variance-csv", default="/tmp/rediscover-variance-decomposition-capstone/variance_decomposition.csv")
    ap.add_argument("--correlation-csv", default="/tmp/rdc-broad-backbone-capstone/broad_backbone_literature_kernel_t2_correlation.csv")
    ap.add_argument("--efg-audit-csv", default="/tmp/rdc-efg-localframe-audit-capstone/efg_localframe_audit.csv")
    ap.add_argument("--efg-distill-csv", default="/tmp/rdc-efg-arc-evidence-localframe-capstone/efg_distill_summary.csv")
    ap.add_argument("--out-dir", default="/tmp/rediscover-capstone-charts")
    return ap.parse_args()


def finite(v):
    return np.isfinite(np.asarray(v, dtype=float))


def save(fig, out_dir, stem, paths):
    png = out_dir / f"{stem}.png"
    svg = out_dir / f"{stem}.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    paths.extend([str(png), str(svg)])


def mechanism_order(values):
    preferred = [
        "ring_current",
        "bond_anisotropy_scalar",
        "ff14sb_charge_field",
        "apbs_efield_buckingham",
        "broad_total_T2",
        "ring_current_T2",
        "bond_anisotropy_T2",
        "ff14sb_charge_EFG_T2",
        "apbs_efg_T2",
        "apbs_efg_T2_random",
        "mopac_coulomb_shielding_T2",
        "mopac_coulomb_shielding_T2_random",
        "mopac_mc_shielding_T2",
        "mopac_mc_shielding_T2_random",
    ]
    seen = list(dict.fromkeys(values))
    return [m for m in preferred if m in seen] + [m for m in seen if m not in preferred]


def plot_variance_shares(var_df, out_dir, paths):
    for target in ["sigma_iso", "T2"]:
        d = var_df[var_df["target"] == target].copy()
        if d.empty:
            continue
        order = mechanism_order(d["mechanism"].tolist())
        d["mechanism"] = pd.Categorical(d["mechanism"], categories=order, ordered=True)
        d["stratum"] = pd.Categorical(d["stratum"], categories=STRATA_ORDER, ordered=True)
        d = d.sort_values(["mechanism", "stratum", "split_strategy"])
        labels = [f"{r.mechanism}  {r.stratum}" for r in d.itertuples()]
        y = np.arange(len(d))
        height = max(6.0, 0.24 * len(d))
        fig, ax = plt.subplots(figsize=(10.8, height))
        between = d["variance_share_between"].to_numpy(float)
        within = d["variance_share_within"].to_numpy(float)
        ax.barh(y, between, color=COLORS["between"], label="between atom means")
        ax.barh(y, within, left=between, color=COLORS["within"], label="within atom motion")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlim(0, 1)
        ax.set_xlabel("share of target variance")
        ax.set_title(f"{target}: between/within variance placement")
        ax.grid(axis="x", color="#d8d8d8", linewidth=0.6)
        ax.legend(loc="lower right", frameon=False)
        ax.invert_yaxis()
        fig.tight_layout()
        save(fig, out_dir, f"variance_shares_{target}", paths)


def heatmap(ax, matrix, rows, cols, title, cmap="RdBu_r", vmin=-1, vmax=1):
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=0)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows)
    ax.set_title(title)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:+.2f}", ha="center", va="center", fontsize=8)
    return im


def pivot_metric(df, row_col, col_col, value_col, rows, cols):
    matrix = np.full((len(rows), len(cols)), np.nan)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            vals = df.loc[(df[row_col] == r) & (df[col_col] == c), value_col].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals):
                matrix[i, j] = vals[0]
    return matrix


def plot_decircularising(corr_df, out_dir, paths):
    d = corr_df[corr_df["stratum"].isin(STRATA_ORDER)].copy()
    if d.empty:
        return
    rows = [k for k in KERNEL_ORDER if k in set(d["kernel"])]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)
    im = heatmap(
        axes[0],
        pivot_metric(d, "kernel", "stratum", "component_r", rows, STRATA_ORDER),
        rows,
        STRATA_ORDER,
        "component r",
    )
    heatmap(
        axes[1],
        pivot_metric(d, "kernel", "stratum", "magnitude_r", rows, STRATA_ORDER),
        rows,
        STRATA_ORDER,
        "|T2| r",
    )
    fig.colorbar(im, ax=axes, shrink=0.86, label="Pearson r")
    save(fig, out_dir, "decircularising_correlations", paths)


def plot_t2_fits(var_df, out_dir, paths):
    d = var_df[(var_df["target"] == "T2") & (var_df["stratum"].isin(STRATA_ORDER))].copy()
    if d.empty:
        return
    d = d[~((d["mechanism"] == "apbs_efg_T2_random") & (d["split_strategy"] == "random"))]
    rows = mechanism_order(d["mechanism"].tolist())
    fig, axes = plt.subplots(2, 2, figsize=(12.4, max(6.0, 0.28 * len(rows))), constrained_layout=True)
    metrics = [
        ("between_LOAO_absT2_r", "between |T2| r"),
        ("within_frame_absT2_r", "within |T2| r"),
        ("between_LOAO_component_r", "between component r"),
        ("within_frame_component_r", "within component r"),
    ]
    im = None
    for ax, (metric, title) in zip(axes.flat, metrics):
        im = heatmap(ax, pivot_metric(d, "mechanism", "stratum", metric, rows, STRATA_ORDER),
                     rows, STRATA_ORDER, title)
    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.86, label="correlation")
    save(fig, out_dir, "t2_fit_metrics", paths)


def plot_efg_fix(audit_df, distill_df, out_dir, paths):
    d = audit_df[audit_df["stratum"].isin(STRATA_ORDER)].copy()
    if d.empty:
        return
    x = np.arange(len(STRATA_ORDER))
    width = 0.36
    lab = d[d["frame"] == "lab"].set_index("stratum").reindex(STRATA_ORDER)
    local = d[d["frame"] == "local"].set_index("stratum").reindex(STRATA_ORDER)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0), constrained_layout=True)
    axes[0].bar(x - width / 2, lab["target_median_lag1_rho"], width, color=COLORS["lab"], label="lab")
    axes[0].bar(x + width / 2, local["target_median_lag1_rho"], width, color=COLORS["local"], label="local")
    axes[0].set_title("target lag-1 rho")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(STRATA_ORDER)
    axes[0].set_ylim(0, 1)
    axes[0].legend(frameon=False)
    axes[0].grid(axis="y", color="#d8d8d8", linewidth=0.6)

    axes[1].bar(x - width / 2, lab["feature_target_component_r"], width, color=COLORS["lab"])
    axes[1].bar(x + width / 2, local["feature_target_component_r"], width, color=COLORS["local"])
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    axes[1].set_title("within component r")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(STRATA_ORDER)
    axes[1].set_ylim(-0.7, 0.7)
    axes[1].grid(axis="y", color="#d8d8d8", linewidth=0.6)

    axes[2].bar(x - width / 2, lab["constant_gamma_R2"], width, color=COLORS["lab"])
    axes[2].bar(x + width / 2, local["constant_gamma_R2"], width, color=COLORS["local"])
    axes[2].axhline(0, color="#555555", linewidth=0.8)
    axes[2].set_title("constant gamma R2")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(STRATA_ORDER)
    axes[2].grid(axis="y", color="#d8d8d8", linewidth=0.6)
    fig.suptitle("APBS EFG: lab-frame rotation confound removed by local-frame emit")
    save(fig, out_dir, "efg_localframe_correction", paths)

    if distill_df is not None and not distill_df.empty:
        e = distill_df.set_index("stratum").reindex(STRATA_ORDER)
        fig, ax = plt.subplots(figsize=(7.2, 4.0))
        ax.bar(x - width / 2, e["nonlinear_test_r2"], width, color="#5b7fb9", label="nonlinear gate R2")
        ax.bar(x + width / 2, e["nonlinear_test_magnitude_r"], width, color="#c98f3c", label="|T2| r")
        ax.axhline(0, color="#555555", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(STRATA_ORDER)
        ax.set_title("corrected local-frame EFG fit")
        ax.legend(frameon=False)
        ax.grid(axis="y", color="#d8d8d8", linewidth=0.6)
        fig.tight_layout()
        save(fig, out_dir, "efg_localframe_fit", paths)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    var_df = pd.read_csv(args.variance_csv)
    corr_df = pd.read_csv(args.correlation_csv)
    audit_df = pd.read_csv(args.efg_audit_csv)
    distill_path = Path(args.efg_distill_csv)
    distill_df = pd.read_csv(distill_path) if distill_path.exists() else None

    plot_variance_shares(var_df, out_dir, paths)
    plot_decircularising(corr_df, out_dir, paths)
    plot_t2_fits(var_df, out_dir, paths)
    plot_efg_fix(audit_df, distill_df, out_dir, paths)

    manifest = {
        "variance_csv": args.variance_csv,
        "correlation_csv": args.correlation_csv,
        "efg_audit_csv": args.efg_audit_csv,
        "efg_distill_csv": args.efg_distill_csv if distill_df is not None else None,
        "figures": paths,
        "reuse_note": "Regenerate the same CSV tables for the 750-DFT substrate and pass their paths to this script; no code edits are required.",
    }
    manifest_path = out_dir / "chart_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"wrote {manifest_path}")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
