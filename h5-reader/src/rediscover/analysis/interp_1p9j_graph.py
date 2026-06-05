#!/usr/bin/env python3
"""Re-render the 1P9J equivariant shielding-tensor prediction graph.

This is a plotting-only reader for the cached run artifacts. It does not train,
fit, or import the e3nn training script.
"""

from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RUN_DIR = Path("/tmp/rediscover-runs/2026-06-04-interp-1p9j-e3nn")
T2_TARGET_COLS = [f"target_T2_{k}_centered" for k in range(5)]
T2_PRED_COLS = [f"pred_T2_{k}_centered" for k in range(5)]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Re-cut the 1P9J advisor graph from cached predictions."
    )
    ap.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    ap.add_argument("--predictions", type=Path, default=None)
    ap.add_argument("--metrics", type=Path, default=None)
    ap.add_argument("--png", type=Path, default=None)
    ap.add_argument("--pdf", type=Path, default=None)
    return ap.parse_args()


def r2_np(pred: np.ndarray, target: np.ndarray) -> float:
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    den = np.sum((target - target.mean(axis=0)) ** 2)
    if den <= 0:
        return float("nan")
    return float(1.0 - np.sum((target - pred) ** 2) / den)


def corr_np(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if len(a) < 2 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def scatter_limits(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    vals = np.concatenate([x.ravel(), y.ravel()])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return -1.0, 1.0
    lo, hi = float(vals.min()), float(vals.max())
    if lo == hi:
        pad = max(1.0, abs(lo) * 0.1)
    else:
        pad = 0.08 * (hi - lo)
    return lo - pad, hi + pad


def metric(metrics: dict, computed: dict, key: str) -> float:
    value = metrics.get(key, computed.get(key, float("nan")))
    return float(value)


def load_metrics(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def require_columns(df: pd.DataFrame, cols: list[str], source: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in {source}: {missing}")


def render_graph(
    df: pd.DataFrame,
    metrics: dict,
    path_png: Path,
    path_pdf: Path,
) -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    test = df["split"].to_numpy() == "test"
    if not np.any(test):
        raise SystemExit("Predictions CSV has no split == 'test' rows.")

    target_t2 = df[T2_TARGET_COLS].to_numpy(dtype=float)
    pred_t2 = df[T2_PRED_COLS].to_numpy(dtype=float)
    target_t0 = df["target_T0_centered"].to_numpy(dtype=float)
    pred_t0 = df["pred_T0_centered"].to_numpy(dtype=float)
    target_t0_restored = df["target_T0_restored"].to_numpy(dtype=float)
    pred_t0_restored = df["pred_T0_restored"].to_numpy(dtype=float)
    target_mag = np.linalg.norm(target_t2, axis=1)
    pred_mag = np.linalg.norm(pred_t2, axis=1)

    computed = {
        "t2_r2_test": r2_np(pred_t2[test], target_t2[test]),
        "t2_mag_r_test": corr_np(target_mag[test], pred_mag[test]),
        "t2_mag_r2_test": r2_np(pred_mag[test], target_mag[test]),
        "t0_mod_r2_test": r2_np(pred_t0[test], target_t0[test]),
        "t0_restored_r2_test": r2_np(
            pred_t0_restored[test], target_t0_restored[test]
        ),
    }
    m = {key: metric(metrics, computed, key) for key in computed}

    atoms = int(metrics.get("atoms", df["atom_index"].nunique()))
    frames = int(metrics.get("frames", df["original_frame_index"].nunique()))
    groups = int(metrics.get("groups", len(df)))
    train_groups = int(metrics.get("train_groups", int((df["split"] == "train").sum())))
    test_groups = int(metrics.get("test_groups", int(test.sum())))
    cross = str(metrics.get("best_cross", "learnable"))

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
        }
    )

    fig, axs = plt.subplots(2, 2, figsize=(11.2, 8.6))
    fig.subplots_adjust(
        left=0.075, right=0.97, bottom=0.075, top=0.90, wspace=0.30, hspace=0.42
    )

    ax = axs[0, 0]
    x = target_t2[test].ravel()
    y = pred_t2[test].ravel()
    lim = scatter_limits(x, y)
    ax.scatter(x, y, s=9, alpha=0.36, linewidths=0, color="#1f5a78")
    ax.plot(lim, lim, color="black", lw=1.0, alpha=0.75)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_title("A. Held-out T2 components")
    ax.set_xlabel("DFT T2 modulation, train-atom centered")
    ax.set_ylabel("equivariant prediction")
    ax.text(
        0.04,
        0.96,
        f"clean R2 = {m['t2_r2_test']:+.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9),
    )

    ax = axs[0, 1]
    lim = scatter_limits(target_mag[test], pred_mag[test])
    ax.scatter(
        target_mag[test],
        pred_mag[test],
        s=16,
        alpha=0.45,
        linewidths=0,
        color="#7d5a1f",
    )
    ax.plot(lim, lim, color="black", lw=1.0, alpha=0.75)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_title("B. |T2| recovery")
    ax.set_xlabel("DFT |T2| modulation")
    ax.set_ylabel("predicted |T2| modulation")
    ax.text(
        0.04,
        0.96,
        f"r = {m['t2_mag_r_test']:+.3f}\nmag R2 = {m['t2_mag_r2_test']:+.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9),
    )

    ax = axs[1, 0]
    lim = scatter_limits(target_t0[test], pred_t0[test])
    ax.scatter(
        target_t0[test],
        pred_t0[test],
        s=16,
        alpha=0.45,
        linewidths=0,
        color="#3a6f44",
    )
    ax.plot(lim, lim, color="black", lw=1.0, alpha=0.75)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_title("C. T0 modulation recovery")
    ax.set_xlabel("DFT T0 modulation, train-atom centered")
    ax.set_ylabel("predicted T0 modulation")
    ax.text(
        0.04,
        0.96,
        f"modulation R2 = {m['t0_mod_r2_test']:+.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=11,
        bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.95),
    )
    iax = inset_axes(ax, width="36%", height="36%", loc="lower right", borderpad=1.2)
    ilim = scatter_limits(target_t0_restored[test], pred_t0_restored[test])
    iax.scatter(
        target_t0_restored[test],
        pred_t0_restored[test],
        s=5,
        alpha=0.35,
        linewidths=0,
        color="#6b6b6b",
    )
    iax.plot(ilim, ilim, color="black", lw=0.7, alpha=0.6)
    iax.set_xlim(ilim)
    iax.set_ylim(ilim)
    iax.set_xticks([])
    iax.set_yticks([])
    iax.set_title("restored sigma", fontsize=8)
    ax.text(
        0.98,
        0.03,
        f"restored R2 = {m['t0_restored_r2_test']:+.3f}",
        transform=ax.transAxes,
        va="bottom",
        ha="right",
        fontsize=8,
        bbox=dict(facecolor="white", edgecolor="0.85", alpha=0.9),
    )

    ax = axs[1, 1]
    ax.axis("off")
    paragraphs = [
        (
            "Software tool: ingest a protein's local geometry and predict its "
            "NMR shielding tensor equivariantly on held-out frames. The tensor "
            "target is the 5-component T2 field, with a companion T0 channel."
        ),
        (
            "Practical reading: for 1P9J, the tool trained on a subset of "
            "frames and predicted the rest of the same protein. Held-out T2 "
            f"component R2 = {m['t2_r2_test']:+.3f}; |T2| r = "
            f"{m['t2_mag_r_test']:+.3f}; T0 modulation R2 = "
            f"{m['t0_mod_r2_test']:+.3f}."
        ),
        (
            "Per-protein bootstrap: DFT a subset of a protein's frames, run "
            "this toolchain, and predict the rest of that protein at about "
            "50% T2 variance recovery. Not great; better than nothing; real."
        ),
        (
            f"Run scale: {atoms} atoms, {frames} frames, {groups:,} "
            f"atom-frame groups ({train_groups:,} train, {test_groups:,} "
            f"held out), {cross} cross term."
        ),
        "Footnote: within-protein result; transferability/BMRB validation is separate.",
    ]
    capability_text = "D. Capability\n\n" + "\n\n".join(
        textwrap.fill(p, width=58) for p in paragraphs
    )
    ax.text(
        0.02,
        0.98,
        capability_text,
        va="top",
        ha="left",
        fontsize=10.1,
        linespacing=1.25,
    )

    fig.suptitle("Equivariant shielding-tensor prediction - 1P9J", fontsize=14)
    fig.savefig(path_png, dpi=220)
    fig.savefig(path_pdf)
    plt.close(fig)
    return m


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir
    predictions = args.predictions or run_dir / "interp_1p9j_predictions.csv"
    metrics_path = args.metrics or run_dir / "interp_1p9j_metrics.json"
    png = args.png or run_dir / "interp_1p9j_advisor_graph.png"
    pdf = args.pdf or run_dir / "interp_1p9j_advisor_graph.pdf"

    df = pd.read_csv(predictions)
    require_columns(
        df,
        [
            "atom_index",
            "original_frame_index",
            "split",
            "target_T0_centered",
            "pred_T0_centered",
            "target_T0_restored",
            "pred_T0_restored",
            *T2_TARGET_COLS,
            *T2_PRED_COLS,
        ],
        predictions,
    )
    metrics = load_metrics(metrics_path)

    png.parent.mkdir(parents=True, exist_ok=True)
    pdf.parent.mkdir(parents=True, exist_ok=True)
    rendered = render_graph(df, metrics, png, pdf)

    print(f"graph_png={png}")
    print(f"graph_pdf={pdf}")
    print(
        "summary: Equivariant shielding-tensor tool predicts held-out 1P9J frames "
        f"from a within-protein training subset; T2 R2={rendered['t2_r2_test']:+.3f}."
    )
    print(
        "summary: Panel C reports T0 modulation on the main axis; restored sigma "
        "is kept as a small inset."
    )


if __name__ == "__main__":
    main()
