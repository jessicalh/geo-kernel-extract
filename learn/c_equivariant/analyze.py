#!/usr/bin/env python3
"""
Cold analytical diagnostic — no learning, pure numpy.

For each protein, examines the raw relationship between per-protein-
normalized kernel T2 and DFT delta T2.  Produces a report on:

1. Per-kernel correlation with target (which kernels carry signal?)
2. Per-element breakdown (where is the signal strongest?)
3. Per-distance breakdown (how does signal decay with ring distance?)
4. Per-protein consistency (any outliers in kernel-target alignment?)
5. Residual structure after optimal ridge fit (what's left to learn?)

Usage:
    python learn/c_equivariant/analyze.py --run CalibrationExtractionTest
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein import load_protein, list_proteins
from features import RingType, BondCategory

REPO = Path(__file__).resolve().parent.parent.parent

KERNEL_NAMES = (
    [f"BS_{RingType(t).name}" for t in range(8)] +
    [f"HM_{RingType(t).name}" for t in range(8)] +
    [f"Disp_{RingType(t).name}" for t in range(8)] +
    [f"MC_{BondCategory(c).name}" for c in range(5)] +
    [f"MopacMC_{BondCategory(c).name}" for c in range(5)] +
    ["MC_total", "Coulomb_total", "HBond_total",
     "RingSusc_total", "MopacCoulomb_total", "MopacMC_total",
     "EFG_bb", "EFG_aro", "MopacEFG_bb", "MopacEFG_aro",
     "APBS_EFG", "DeltaAPBS_EFG"]
)

ELEMENT_NAMES = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
DIST_BINS = [(0, 4), (4, 8), (8, 12), (12, 999)]
DIST_LABELS = ["0-4A", "4-8A", "8-12A", "12+A"]


def _corr(a, b):
    """Pearson correlation between two 1D arrays."""
    if len(a) < 3:
        return float("nan")
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt((a ** 2).sum() * (b ** 2).sum())
    if denom < 1e-12:
        return 0.0
    return float((a * b).sum() / denom)


def _r2(pred, tgt):
    """R² on flattened arrays."""
    ss_res = np.sum((tgt - pred) ** 2)
    ss_tot = np.sum((tgt - tgt.mean(axis=0, keepdims=True)) ** 2)
    if ss_tot < 1e-12:
        return 0.0
    return 1.0 - ss_res / ss_tot


def load_one(protein_dir):
    """Load one protein, return per-protein normalized kernels + target + metadata."""
    from c_equivariant.dataset import _build_features
    result, err = _build_features(protein_dir)
    if result is None:
        return None
    p = load_protein(protein_dir)
    idx = np.where(p.delta.scalars.matched_mask)[0]
    return {
        "kernels": result["kernels"],    # (M, N_KERNELS, 5) per-protein normalized
        "target": result["target"],       # (M, 5)
        "ring_dist": result["ring_dist"], # (M,)
        "element": p.element[idx],
        "residue_type": p.residue_type[idx],
        "protein_id": protein_dir.name,
    }


def analyze(features_dir: Path):
    proteins = list_proteins(features_dir)
    print(f"Loading {len(proteins)} proteins...")

    all_data = []
    for pid in proteins:
        d = load_one(features_dir / pid)
        if d is not None:
            all_data.append(d)
    print(f"Loaded {len(all_data)} proteins with delta\n")

    if not all_data:
        print("No proteins with delta found.")
        return

    # Stack everything
    kernels = np.vstack([d["kernels"] for d in all_data])
    target = np.vstack([d["target"] for d in all_data])
    ring_dist = np.concatenate([d["ring_dist"] for d in all_data])
    element = np.concatenate([d["element"] for d in all_data])

    N = len(target)
    n_kernels = kernels.shape[1]
    print(f"Total atoms: {N}, kernels: {n_kernels}")
    print(f"Target T2 std: {target.std():.4f} ppm")
    print(f"Target T2 mean magnitude: {np.linalg.norm(target, axis=1).mean():.4f} ppm\n")

    # ── 1. Per-kernel correlation with target ────────────────────
    print("=" * 70)
    print("1. PER-KERNEL CORRELATION WITH DFT DELTA T2")
    print("   (Pearson r on flattened T2 components, per-protein normalized)")
    print("-" * 70)

    kernel_corrs = []
    for k in range(n_kernels):
        r = _corr(kernels[:, k, :].ravel(), target.ravel())
        kernel_corrs.append((KERNEL_NAMES[k], r))

    kernel_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for name, r in kernel_corrs:
        bar = "#" * int(abs(r) * 50)
        sign = "+" if r > 0 else "-"
        active = kernels[:, KERNEL_NAMES.index(name), :].std() > 1e-10
        flag = "" if active else " [ZERO]"
        print(f"  {name:30s}  r={r:+.4f}  {sign}{bar}{flag}")

    # Top contributors
    top5 = kernel_corrs[:5]
    print(f"\n  Top 5: {', '.join(f'{n}({r:+.3f})' for n, r in top5)}")

    # ── 2. Per-element breakdown ─────────────────────────────────
    print(f"\n{'=' * 70}")
    print("2. PER-ELEMENT TARGET MAGNITUDE AND KERNEL CORRELATION")
    print("-" * 70)

    for elem_id, elem_name in sorted(ELEMENT_NAMES.items()):
        mask = element == elem_id
        n = mask.sum()
        if n < 5:
            continue
        t = target[mask]
        mag = np.linalg.norm(t, axis=1)
        # Best kernel for this element
        best_r = 0
        best_k = ""
        for k in range(n_kernels):
            kdat = kernels[mask, k, :]
            if kdat.std() < 1e-10:
                continue
            r = _corr(kdat.ravel(), t.ravel())
            if abs(r) > abs(best_r):
                best_r = r
                best_k = KERNEL_NAMES[k]
        print(f"  {elem_name:2s}  n={n:6d}  T2_mag={mag.mean():.3f}±{mag.std():.3f} ppm  "
              f"best_kernel={best_k}(r={best_r:+.3f})")

    # ── 3. Per-distance breakdown ────────────────────────────────
    print(f"\n{'=' * 70}")
    print("3. SIGNAL BY DISTANCE TO NEAREST REMOVED RING")
    print("-" * 70)

    for (lo, hi), label in zip(DIST_BINS, DIST_LABELS):
        mask = (ring_dist >= lo) & (ring_dist < hi)
        n = mask.sum()
        if n < 5:
            print(f"  {label:6s}  n={n:6d}  (too few)")
            continue
        t = target[mask]
        mag = np.linalg.norm(t, axis=1)

        # Ridge fit on this distance band
        X = kernels[mask].reshape(n, -1)
        y = t
        try:
            lam = 1e-2
            XtX = X.T @ X + lam * np.eye(X.shape[1])
            Xty = X.T @ y
            w = np.linalg.solve(XtX, Xty)
            ridge_pred = X @ w
            r2 = _r2(ridge_pred, y)
        except np.linalg.LinAlgError:
            r2 = float("nan")

        print(f"  {label:6s}  n={n:6d}  T2_mag={mag.mean():.3f}±{mag.std():.3f} ppm  "
              f"ridge_R²={r2:.3f}")

    # ── 4. Per-protein consistency ───────────────────────────────
    print(f"\n{'=' * 70}")
    print("4. PER-PROTEIN RIDGE R² (optimal linear mixing, no MLP)")
    print("-" * 70)

    protein_r2s = []
    for d in all_data:
        k = d["kernels"]
        t = d["target"]
        M = len(t)
        if M < 20:
            continue
        X = k.reshape(M, -1)
        try:
            lam = 1e-2
            XtX = X.T @ X + lam * np.eye(X.shape[1])
            w = np.linalg.solve(XtX, X.T @ t)
            pred = X @ w
            r2 = _r2(pred, t)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        protein_r2s.append((d["protein_id"], r2, M))

    protein_r2s.sort(key=lambda x: x[1])
    vals = [r for _, r, _ in protein_r2s if not np.isnan(r)]
    print(f"  Proteins: {len(protein_r2s)}")
    print(f"  Ridge R²: median={np.median(vals):.3f}  "
          f"mean={np.mean(vals):.3f}  std={np.std(vals):.3f}")
    print(f"  Range: [{min(vals):.3f}, {max(vals):.3f}]")
    print(f"\n  Worst 5:")
    for pid, r2, m in protein_r2s[:5]:
        print(f"    {pid:20s}  R²={r2:.3f}  atoms={m}")
    print(f"  Best 5:")
    for pid, r2, m in protein_r2s[-5:]:
        print(f"    {pid:20s}  R²={r2:.3f}  atoms={m}")

    # ── 5. Greedy forward kernel selection ─────────────────────────
    print(f"\n{'=' * 70}")
    print("5. GREEDY FORWARD KERNEL SELECTION")
    print("   (add one kernel at a time, pick the one that maximises ridge R²)")
    print("-" * 70)

    n_kernels = kernels.shape[1]
    lam = 1e-2
    remaining = set(range(n_kernels))
    selected = []
    cumulative_r2 = []

    for step in range(min(n_kernels, 20)):  # stop at 20 — diminishing returns
        best_k = -1
        best_r2 = -np.inf
        for k in remaining:
            trial = selected + [k]
            X = kernels[:, trial, :].reshape(N, -1)
            try:
                XtX = X.T @ X + lam * np.eye(X.shape[1])
                w = np.linalg.solve(XtX, X.T @ target)
                pred = X @ w
                r2 = _r2(pred, target)
            except np.linalg.LinAlgError:
                r2 = -np.inf
            if r2 > best_r2:
                best_r2 = r2
                best_k = k
        if best_k < 0:
            break
        selected.append(best_k)
        remaining.remove(best_k)
        cumulative_r2.append(best_r2)
        delta = best_r2 - (cumulative_r2[-2] if len(cumulative_r2) > 1 else 0.0)
        print(f"  {step+1:2d}. +{KERNEL_NAMES[best_k]:30s}  "
              f"R²={best_r2:.4f}  (+{delta:.4f})")

    print(f"\n  Summary: {len(selected)} kernels shown, "
          f"first 5 reach R²={cumulative_r2[4]:.4f}, "
          f"all {len(selected)} reach R²={cumulative_r2[-1]:.4f}")

    # ── 6. Leave-one-out kernel ablation (from full set) ─────────
    print(f"\n{'=' * 70}")
    print("6. LEAVE-ONE-OUT KERNEL ABLATION")
    print("   (R² drop when each kernel is removed from the full ridge)")
    print("-" * 70)

    X_all = kernels.reshape(N, -1)
    try:
        XtX = X_all.T @ X_all + lam * np.eye(X_all.shape[1])
        w_full = np.linalg.solve(XtX, X_all.T @ target)
        r2_full = _r2(X_all @ w_full, target)
    except np.linalg.LinAlgError:
        r2_full = float("nan")

    print(f"  Full ridge R² ({n_kernels} kernels): {r2_full:.4f}\n")

    ablations = []
    for k in range(n_kernels):
        keep = [j for j in range(n_kernels) if j != k]
        X = kernels[:, keep, :].reshape(N, -1)
        try:
            XtX = X.T @ X + lam * np.eye(X.shape[1])
            w = np.linalg.solve(XtX, X.T @ target)
            r2 = _r2(X @ w, target)
        except np.linalg.LinAlgError:
            r2 = float("nan")
        drop = r2_full - r2
        ablations.append((KERNEL_NAMES[k], drop, r2))

    ablations.sort(key=lambda x: x[1], reverse=True)
    for name, drop, r2_without in ablations:
        if drop < 1e-5:
            bar = ""
        else:
            bar = "#" * min(int(drop * 500), 50)
        print(f"  {name:30s}  drop={drop:+.4f}  R²_without={r2_without:.4f}  {bar}")

    unique = [n for n, d, _ in ablations if d > 0.005]
    redundant = [n for n, d, _ in ablations if d < 0.0005]
    print(f"\n  Unique (drop > 0.005):    {', '.join(unique) if unique else 'none'}")
    print(f"  Redundant (drop < 0.0005): {len(redundant)} kernels")

    # ── 7. Residual subspace projection ─────────────────────────────
    print(f"\n{'=' * 70}")
    print("7. RESIDUAL SUBSPACE PROJECTION")
    print("   (is the unexplained 20% in-span or out-of-span of the kernels?)")
    print("-" * 70)

    try:
        ridge_pred = X_all @ w_full
        residual = target - ridge_pred

        print(f"  Global ridge R²: {r2_full:.4f}")
        print(f"  Residual std:    {residual.std():.4f} ppm")

        # Project residual onto each kernel direction per atom
        print(f"\n  Per-kernel alignment with residual (mean |cos angle|):")
        kernel_alignments = []
        for k in range(n_kernels):
            kv = kernels[:, k, :]  # (N, 5)
            k_norm = np.linalg.norm(kv, axis=1, keepdims=True)
            r_norm = np.linalg.norm(residual, axis=1, keepdims=True)
            valid = ((k_norm > 1e-10) & (r_norm > 1e-10)).ravel()
            if valid.sum() < 10:
                kernel_alignments.append((KERNEL_NAMES[k], 0.0, 0))
                continue
            cos = np.sum(kv[valid] * residual[valid], axis=1) / (
                k_norm[valid].ravel() * r_norm[valid].ravel())
            kernel_alignments.append((KERNEL_NAMES[k], np.mean(np.abs(cos)), int(valid.sum())))

        kernel_alignments.sort(key=lambda x: x[1], reverse=True)
        for name, mean_cos, nv in kernel_alignments[:15]:
            bar = "#" * int(mean_cos * 50)
            print(f"    {name:30s}  |cos|={mean_cos:.3f}  n={nv:6d}  {bar}")

        # Aggregate: what fraction of residual variance is in-span?
        # Compare OLS R² (no regularisation) vs ridge R². The gap is
        # in-span signal that ridge suppresses; (1 - OLS R²) is truly
        # out-of-span. Uses lstsq which handles rank-deficiency via SVD
        # on the smaller (K*5, K*5) normal equations, not the full matrix.
        ols_pred = np.zeros_like(target)
        for comp in range(5):
            w_ols, res_ols, rank_ols, sv = np.linalg.lstsq(X_all, target[:, comp], rcond=1e-8)
            ols_pred[:, comp] = X_all @ w_ols
        r2_ols = _r2(ols_pred, target)
        ols_residual = target - ols_pred
        in_span_frac = (r2_ols - r2_full) / (1 - r2_full) if r2_full < 1.0 else 0.0
        print(f"\n  OLS R² (no regularisation): {r2_ols:.4f}")
        print(f"  Ridge R² (λ=0.01):         {r2_full:.4f}")
        print(f"  Effective rank: {rank_ols} (of {X_all.shape[1]})")
        print(f"  In-span residual (OLS-ridge gap): {r2_ols - r2_full:.4f} "
              f"({in_span_frac:.1%} of unexplained)")
        print(f"  Out-of-span (1 - OLS R²):        {1 - r2_ols:.4f}")

        if in_span_frac > 0.15:
            print(f"\n  >> Ridge is leaving {in_span_frac:.0%} of the residual on the table."
                  f"\n     The MLP (environment-dependent weights) should recover this.")
        else:
            print(f"\n  >> Regularisation costs little. The missing {1-r2_ols:.0%}"
                  f" is truly out-of-span.\n     Needs new physics or nonlinear"
                  f" kernel interactions (correction head).")

        # Residual by element
        print(f"\n  Residual by element:")
        for elem_id, elem_name in sorted(ELEMENT_NAMES.items()):
            mask = element == elem_id
            n = mask.sum()
            if n < 5:
                continue
            res = residual[mask]
            print(f"    {elem_name:2s}  n={n:6d}  "
                  f"residual_mag={np.linalg.norm(res, axis=1).mean():.3f} ppm  "
                  f"bias={res.mean():.4f}")

        # Residual by distance
        print(f"\n  Residual by distance:")
        for (lo, hi), label in zip(DIST_BINS, DIST_LABELS):
            mask = (ring_dist >= lo) & (ring_dist < hi)
            n = mask.sum()
            if n < 5:
                continue
            res = residual[mask]
            print(f"    {label:6s}  n={n:6d}  "
                  f"residual_mag={np.linalg.norm(res, axis=1).mean():.3f} ppm  "
                  f"bias={res.mean():.4f}")

    except np.linalg.LinAlgError:
        print("  Ridge fit failed (singular matrix)")

    print(f"\n{'=' * 70}")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(description="Analytical diagnostic")
    parser.add_argument("--run", type=str, default="CalibrationExtractionTest",
                        help="extraction run name")
    args = parser.parse_args()
    features_dir = REPO / "calibration" / "features" / args.run
    analyze(features_dir)


if __name__ == "__main__":
    main()
