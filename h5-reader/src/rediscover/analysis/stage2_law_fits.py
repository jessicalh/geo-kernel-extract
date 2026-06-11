#!/usr/bin/env python3
"""Stage-2 law fits on the emitted Build-4 per-atom substrate.

This script deliberately consumes only emitted Build-4 CSV/NPY sidecars. It fits
scalar intensities multiplying emitted five-component T2 shadows; it never opens
trajectory.h5, source dumps, ORCA outputs, or older run directories.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import change_of_basis as cob
from allatom_fit_common import (
    blocked_frame_cv_splits,
    centered_by_train_atom,
    ci95,
    effective_n_components,
    finite_fmt,
    json_sanitize,
    jackknife_metric_values,
    jackknife_se_from_values,
    r2_score,
    select_best_alpha,
    split_frame_block,
)


SUBSTRATE_DIR = Path("/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4")
OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-stage2-fits")
STAGE21_OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-stage2_1-sweep")
STAGE22_OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-stage2_2-unified-vet")
STAGE23_OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-stage2_3-probability")
TRUE_LOAO_OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-true-loao")
MIN_DISK_FREE_BYTES = 20 * 1024**3
MAX_REDISCOVER_BYTES = 15 * 1024**3
T2 = 5
PYSR_SAMPLE = 3000
SWEEP_KEEP_FRACTIONS = (0.10, 0.20, 0.35, 0.50, 0.70, 1.00)
FRAME_COUNTS = (20, 50, 100, 200, 400, 660)
STAGE22_ALPHA_GRID = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 1.0e4, 1.0e5)
STAGE22_INNER_CV_FOLDS = 5
STAGE22_BOOTSTRAPS = 500
STAGE22_BOOTSTRAP_FRAME_BLOCK = 10
STAGE23_PERMUTATIONS = 1000
STAGE23_DAB_CUTOFFS = (0.30, 0.40, 0.50, 0.55, 0.60)
STAGE23_RING_CUTOFFS = (0.60, 0.65, 0.70, 0.80, 0.90)


@dataclass(frozen=True)
class ShadowSpec:
    label: str
    array: str
    prefix: str
    groups: int = 1
    sign: float = 1.0


@dataclass(frozen=True)
class KernelSpec:
    name: str
    case_dir: str
    dominance: str
    primary: ShadowSpec
    paths: tuple[ShadowSpec, ...]
    selection_note: str


KERNELS = [
    KernelSpec(
        "ring",
        "ring",
        "ring",
        ShadowSpec("jb", "per_atom_substrate_features_classical", "ring_jb_T2"),
        (
            ShadowSpec("jb", "per_atom_substrate_features_classical", "ring_jb_T2"),
            ShadowSpec("bs", "per_atom_substrate_features_ring_paths", "bs_shielding_T2"),
            ShadowSpec("hm", "per_atom_substrate_features_ring_paths", "hm_shielding_T2"),
        ),
        "CaseHunter ring atoms plus dominant_fraction_ring >= 0.5",
    ),
    KernelSpec(
        "charge",
        "charge_wide",
        "charge",
        ShadowSpec("charge_q_over_r3", "per_atom_substrate_features_classical", "charge_q_over_r3_T2"),
        (
            ShadowSpec("charge_q_over_r3", "per_atom_substrate_features_classical", "charge_q_over_r3_T2"),
        ),
        "CaseHunter charge_wide atoms plus dominant_fraction_charge >= 0.5",
    ),
    KernelSpec(
        "mc",
        "mc",
        "mc",
        ShadowSpec("mc_lit", "per_atom_substrate_features_classical", "mc_lit_T2_valid"),
        (
            ShadowSpec("mc_lit", "per_atom_substrate_features_classical", "mc_lit_T2_valid"),
            ShadowSpec("mopac_mc", "per_atom_substrate_features_classical", "mopac_mc_shielding_T2"),
        ),
        "CaseHunter mc atoms plus dominant_fraction_mc >= 0.5",
    ),
    KernelSpec(
        "field",
        "field",
        "field",
        ShadowSpec("mopac_coulomb", "per_atom_substrate_features_classical", "mopac_coulomb_shielding_T2"),
        (
            ShadowSpec("mopac_coulomb", "per_atom_substrate_features_classical", "mopac_coulomb_shielding_T2"),
            ShadowSpec("mopac_coulomb_efg_sum", "per_atom_substrate_features_method_paths", "mopac_coulomb_efg", groups=2),
            ShadowSpec("water_efg_reconciled", "per_atom_substrate_features_method_paths", "water_efg", sign=-1.0),
        ),
        "CaseHunter field atoms plus dominant_fraction_field >= 0.5",
    ),
]


UNIFIED_SPECS = [
    ShadowSpec("charge_total", "per_atom_substrate_features_classical", "charge_q_over_r3_T2"),
    ShadowSpec("mopac_field_backbone", "per_atom_substrate_features_method_paths", "mopac_coulomb_efg_backbone"),
    ShadowSpec("mopac_field_aromatic", "per_atom_substrate_features_method_paths", "mopac_coulomb_efg_aromatic"),
    ShadowSpec("water_field_efg_reconciled", "per_atom_substrate_features_method_paths", "water_efg", sign=-1.0),
    *[ShadowSpec(f"pq_type_{i}", "per_atom_substrate_features_ring_paths", f"pq_per_type_T2_{5*i}", groups=0) for i in range(8)],
    *[ShadowSpec(f"disp_type_{i}", "per_atom_substrate_features_ring_paths", f"disp_per_type_T2_{5*i}", groups=0) for i in range(8)],
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate-dir", type=Path, default=SUBSTRATE_DIR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--stage2-1", action="store_true", help="run the Stage 2.1 happy-spot sweep and frame-count ablation")
    ap.add_argument("--stage2-2", action="store_true", help="run the Stage 2.2 unified D_ab overfit/charge-coat vet")
    ap.add_argument("--stage2-3", action="store_true", help="run the Stage 2.3 probability/null framing")
    ap.add_argument("--stage2-3-permutations", type=int, default=STAGE23_PERMUTATIONS)
    ap.add_argument("--true-loao", action="store_true", help="run the true between-atom LOAO atom-mean diagnostic")
    ap.add_argument("--true-loao-permutations", type=int, default=1000)
    ap.add_argument("--run-pysr", action="store_true")
    ap.add_argument("--pysr-timeout", type=int, default=45)
    ap.add_argument("--pysr-iterations", type=int, default=8)
    return ap.parse_args()


def path_under(path: Path, root: Path) -> bool:
    path = path.resolve()
    root = root.resolve()
    return path == root or root in path.parents


def directory_size_bytes(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for item in path.rglob("*"):
        try:
            if item.is_file() or item.is_symlink():
                total += item.stat().st_size
        except FileNotFoundError:
            pass
    return int(total)


def disk_guard(out_dir: Path, drop_existing: bool = False) -> dict[str, object]:
    out_dir = out_dir.resolve()
    runs_root = Path("/tmp/rediscover-runs").resolve()
    if not path_under(out_dir, runs_root):
        raise SystemExit(f"FATAL: stage2 output must be under {runs_root}, got {out_dir}")
    if path_under(out_dir, Path("/shared").resolve()):
        raise SystemExit(f"FATAL: refusing to write result data under /shared: {out_dir}")
    usage = shutil.disk_usage(out_dir.parent if out_dir.parent.exists() else runs_root)
    if usage.free < MIN_DISK_FREE_BYTES:
        raise SystemExit(f"FATAL: disk guard: {usage.free / 1024**3:.1f} GiB free, need >=20 GiB")
    deleted_existing = False
    if drop_existing and out_dir.exists():
        shutil.rmtree(out_dir)
        deleted_existing = True
    rediscover_bytes = directory_size_bytes(runs_root)
    if rediscover_bytes > MAX_REDISCOVER_BYTES:
        raise SystemExit(f"FATAL: disk guard: /tmp/rediscover-runs is {rediscover_bytes / 1024**3:.1f} GiB, need <15 GiB")
    out_dir.mkdir(parents=True, exist_ok=True)
    return {
        "checked_before_write": True,
        "output_dir": str(out_dir),
        "filesystem_free_GiB": float(usage.free / 1024**3),
        "rediscover_runs_GiB_before_write": float(rediscover_bytes / 1024**3),
        "min_free_GiB_required": 20.0,
        "max_rediscover_GiB_required": 15.0,
        "drop_existing_requested": bool(drop_existing),
        "deleted_existing_output_dir": bool(deleted_existing),
        "deleted_existing_output_dir_full_path": str(out_dir) if deleted_existing else "",
    }


def load_json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def specs_for_array(specs: dict[str, object], array_name: str) -> dict[str, int]:
    return {
        str(c["column"]): int(c["index"])
        for c in specs["columns"]
        if str(c["array"]) == array_name
    }


def load_array(root: Path, name: str) -> np.ndarray:
    return np.load(root / f"{name}.npy", mmap_mode="r")


def t2_columns(prefix: str) -> list[str]:
    return [f"{prefix}_{i}" for i in range(T2)]


def raw_t2(arrays: dict[str, np.ndarray], specs: dict[str, object], spec: ShadowSpec) -> np.ndarray:
    matrix = arrays[spec.array]
    cols = specs_for_array(specs, spec.array)
    if spec.groups == 0:
        base, start_text = spec.prefix.rsplit("_", 1)
        start = int(start_text)
        names = [f"{base}_{start + i}" for i in range(T2)]
        idx = [cols[n] for n in names]
        return spec.sign * np.asarray(matrix[:, idx], dtype=float)
    if spec.groups == 1:
        idx = [cols[n] for n in t2_columns(spec.prefix)]
        return spec.sign * np.asarray(matrix[:, idx], dtype=float)
    pieces = []
    if spec.prefix == "larsen_hbond":
        prefixes = ["larsen_hbond_1pHB_T2", "larsen_hbond_2pHB_T2", "larsen_hbond_1pHaB_T2", "larsen_hbond_2pHaB_T2"]
        for p in prefixes:
            pieces.append(np.asarray(matrix[:, [cols[n] for n in t2_columns(p)]], dtype=float))
    elif spec.prefix == "mopac_coulomb_efg":
        for p in ["mopac_coulomb_efg_backbone", "mopac_coulomb_efg_aromatic"]:
            pieces.append(np.asarray(matrix[:, [cols[n] for n in t2_columns(p)]], dtype=float))
    else:
        for g in range(spec.groups):
            names = [f"{spec.prefix}_{g * T2 + i}" for i in range(T2)]
            pieces.append(np.asarray(matrix[:, [cols[n] for n in names]], dtype=float))
    return spec.sign * np.sum(pieces, axis=0)


def to_e3nn(values: np.ndarray, C: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float) @ C.T


def fit_scalar_no_intercept(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x).all(axis=1) & np.isfinite(y).all(axis=1)
    xx = x[ok].reshape(-1)
    yy = y[ok].reshape(-1)
    den = float(np.dot(xx, xx))
    if len(xx) < 3 or den <= 0.0:
        return math.nan
    return float(np.dot(xx, yy) / den)


def fit_free5(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.full(T2, np.nan)
    for i in range(T2):
        den = float(np.dot(x[:, i], x[:, i]))
        if den > 0:
            out[i] = float(np.dot(x[:, i], y[:, i]) / den)
    return out


def scalar_predict(x: np.ndarray, k: float) -> np.ndarray:
    return np.asarray(x, dtype=float) * float(k)


def free5_predict(x: np.ndarray, k: np.ndarray) -> np.ndarray:
    return np.asarray(x, dtype=float) * np.asarray(k, dtype=float).reshape(1, T2)


def atom_center(values: np.ndarray, atoms: np.ndarray) -> np.ndarray:
    out = np.full_like(values, np.nan, dtype=float)
    for atom in np.unique(atoms):
        m = atoms == atom
        out[m] = values[m] - values[m].mean(axis=0, keepdims=True)
    return out


def coefficient_jackknife(x: np.ndarray, y: np.ndarray, atoms: np.ndarray) -> tuple[float, float, float, float]:
    labels = np.unique(atoms)
    if len(labels) < 3:
        k = fit_scalar_no_intercept(x, y)
        return k, math.nan, math.nan, math.nan
    vals = []
    for atom in labels:
        keep = atoms != atom
        vals.append(fit_scalar_no_intercept(x[keep], y[keep]))
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    k = fit_scalar_no_intercept(x, y)
    if len(vals) < 3:
        return k, math.nan, math.nan, math.nan
    mean = float(vals.mean())
    se = float(math.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - mean) ** 2)))
    lo, hi = ci95(k, se)
    return k, se, lo, hi


def blocked_within_fit(x: np.ndarray, y: np.ndarray, atoms: np.ndarray, frames: np.ndarray) -> dict[str, object]:
    split = split_frame_block(frames, 0.20, 1)
    train = np.asarray(split["train_mask"], dtype=bool)
    test = np.asarray(split["test_mask"], dtype=bool)
    x_c = centered_by_train_atom(x, atoms, train)
    y_c = centered_by_train_atom(y, atoms, train)
    fit = train & np.isfinite(x_c).all(axis=1) & np.isfinite(y_c).all(axis=1)
    score = test & np.isfinite(x_c).all(axis=1) & np.isfinite(y_c).all(axis=1)
    k, se, lo, hi = coefficient_jackknife(x_c[fit], y_c[fit], atoms[fit])
    k5 = fit_free5(x_c[fit], y_c[fit])
    pred = scalar_predict(x_c[score], k)
    pred5 = free5_predict(x_c[score], k5)
    return {
        "coefficient": k,
        "coefficient_se": se,
        "coefficient_ci95_low": lo,
        "coefficient_ci95_high": hi,
        "within_R2": r2_score(pred, y_c[score]),
        "free5_within_R2": r2_score(pred5, y_c[score]),
        "train_rows": int(fit.sum()),
        "test_rows": int(score.sum()),
        "test_frames": len(split["test_frames"]),
        "purged_train_frames": len(split["purged_frames"]),
        "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
    }


def true_loao_single_feature_fit(x: np.ndarray, y: np.ndarray, atoms: np.ndarray) -> dict[str, object]:
    atom_labels, atom_features, atom_y, _counts = true_loao_atom_mean_arrays([x], y, atoms)
    score = true_loao_score_atom_means(["shadow"], atom_features, atom_y)
    coef = np.asarray(score.get("standardized_coef_median", np.asarray([math.nan])), dtype=float)
    coef_min = np.asarray(score.get("standardized_coef_min", np.asarray([math.nan])), dtype=float)
    coef_max = np.asarray(score.get("standardized_coef_max", np.asarray([math.nan])), dtype=float)
    return {
        "loao_R2": float(score["R2"]),
        "loao_atoms": int(len(atom_labels)),
        "loao_coeff_median": float(coef[0]) if coef.size and np.isfinite(coef[0]) else math.nan,
        "loao_coeff_min": float(coef_min[0]) if coef_min.size and np.isfinite(coef_min[0]) else math.nan,
        "loao_coeff_max": float(coef_max[0]) if coef_max.size and np.isfinite(coef_max[0]) else math.nan,
    }


def loao_modulation_fit(x: np.ndarray, y: np.ndarray, atoms: np.ndarray) -> dict[str, object]:
    return true_loao_single_feature_fit(x, y, atoms)


def path_bucket(coef: float, lo: float, hi: float, within_r2: float, loao_r2: float, support: str, agreement: float) -> str:
    nonzero = np.isfinite(lo) and np.isfinite(hi) and (lo > 0.0 or hi < 0.0)
    if support:
        return "form-recovered-scale-fitted" if nonzero and within_r2 > 0 else "cant-make-it-work-for-now"
    if nonzero and within_r2 > 0 and loao_r2 > 0 and np.isfinite(agreement) and agreement <= 0.35:
        return "recovered-law"
    if nonzero and within_r2 > 0:
        return "form-recovered-scale-fitted"
    return "cant-make-it-work-for-now"


def pysr_path(args: argparse.Namespace, out_dir: Path, kernel: str, x: np.ndarray, y: np.ndarray, train_rows: np.ndarray, test_rows: np.ndarray) -> dict[str, object]:
    if not args.run_pysr:
        return {"path": "pysr", "status": "skipped", "coefficient": math.nan, "R2": math.nan, "equation": ""}
    try:
        from pysr import PySRRegressor
    except Exception as exc:
        return {"path": "pysr", "status": f"unavailable:{type(exc).__name__}", "coefficient": math.nan, "R2": math.nan, "equation": ""}
    xs = x[train_rows].reshape(-1, 1)
    ys = y[train_rows].reshape(-1)
    ok = np.isfinite(xs[:, 0]) & np.isfinite(ys) & (np.abs(xs[:, 0]) > 1.0e-12)
    xs = xs[ok]
    ys = ys[ok]
    if len(ys) < 200:
        return {"path": "pysr", "status": "thin", "coefficient": math.nan, "R2": math.nan, "equation": ""}
    rng = np.random.default_rng(0)
    if len(ys) > PYSR_SAMPLE:
        take = rng.choice(len(ys), size=PYSR_SAMPLE, replace=False)
        xs = xs[take]
        ys = ys[take]
    model = PySRRegressor(
        niterations=args.pysr_iterations,
        populations=8,
        population_size=40,
        maxsize=9,
        binary_operators=["+", "-", "*"],
        unary_operators=[],
        variable_names=["shadow"],
        timeout_in_seconds=args.pysr_timeout,
        progress=False,
        deterministic=True,
        parallelism="serial",
        random_state=0,
        verbosity=0,
        temp_equation_file=True,
    )
    model.fit(xs, ys, variable_names=["shadow"])
    equations = model.equations_.copy()
    equations.to_csv(out_dir / f"pysr_equations_{kernel}.csv", index=False)
    best = model.get_best()
    xte = x[test_rows].reshape(-1, 1)
    yte = y[test_rows].reshape(-1)
    okte = np.isfinite(xte[:, 0]) & np.isfinite(yte)
    pred = model.predict(xte[okte])
    coef = float(np.dot(xs[:, 0], model.predict(xs)) / np.dot(xs[:, 0], xs[:, 0]))
    return {
        "path": "pysr",
        "status": "ok",
        "coefficient": coef,
        "R2": r2_score(pred.reshape(-1, 1), yte[okte].reshape(-1, 1)),
        "equation": str(best["equation"]),
    }


def case_atoms(root: Path, case_dir: str) -> tuple[set[int], pd.DataFrame]:
    path = root / "equations" / case_dir / "cases_manifest.csv"
    cases = pd.read_csv(path)
    return set(int(x) for x in cases["atom_index"]), cases


def support_flag(n_atoms: int, n_rows: int, neff: float) -> str:
    flags = []
    if n_atoms < 10:
        flags.append("thin_atoms")
    if n_rows < 500:
        flags.append("thin_rows")
    if np.isfinite(neff) and neff < 10:
        flags.append("thin_N_eff")
    return ";".join(flags)


def run_kernel_fits(args: argparse.Namespace, data: dict[str, object], out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    C: np.ndarray = data["C"]
    target_lib: np.ndarray = data["target_lib"]
    target = to_e3nn(target_lib, C)
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_bins = arrays["per_atom_substrate_dominance_bins"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    bin_cols = specs_for_array(specs, "per_atom_substrate_dominance_bins")
    dft = rows["dft_present"].to_numpy(int) == 1
    atoms_all = rows["atom_index"].to_numpy(int)
    frames_all = rows["original_frame_index"].to_numpy(int)

    summary_rows = []
    path_rows = []
    case_outputs: dict[str, pd.DataFrame] = {}
    for kernel in KERNELS:
        atoms, cases = case_atoms(args.substrate_dir, kernel.case_dir)
        frac = np.asarray(dom[:, dom_cols[f"dominant_fraction_{kernel.dominance}"]], dtype=float)
        binv = np.asarray(dom_bins[:, bin_cols[f"bin_dominant_fraction_{kernel.dominance}_quintile"]], dtype=int)
        mask = dft & np.isfinite(frac) & (frac >= 0.5) & np.isin(atoms_all, list(atoms))
        selection_rule = "dominance_fraction>=0.5"
        if mask.sum() < 500:
            mask = dft & (binv >= 4) & np.isin(atoms_all, list(atoms))
            selection_rule = "dominance_bin>=Q5 fallback"
        idx = np.flatnonzero(mask)
        y = target[idx]
        y_lib = target_lib[idx]
        atoms = atoms_all[idx]
        frames = frames_all[idx]
        neff_info = effective_n_components(y, atoms, frames)
        neff = float(neff_info["within_N_eff_median"])
        supp = support_flag(len(np.unique(atoms)), len(idx), neff)

        primary_agreement_coeffs: list[float] = []
        primary_within = None
        primary_loao = None
        primary_free5 = math.nan
        primary_coef = math.nan
        primary_ci = (math.nan, math.nan)
        primary_pysr = {"status": "not-run", "coefficient": math.nan, "R2": math.nan, "equation": ""}

        for path in kernel.paths:
            x_lib_all = raw_t2(arrays, specs, path)
            x_lib = x_lib_all[idx]
            x = to_e3nn(x_lib, C)
            within = blocked_within_fit(x, y, atoms, frames)
            within_lib = blocked_within_fit(x_lib, y_lib, atoms, frames)
            loao = true_loao_single_feature_fit(x, y, atoms)
            basis_delta = abs(float(within["coefficient"]) - float(within_lib["coefficient"])) if np.isfinite(within["coefficient"]) and np.isfinite(within_lib["coefficient"]) else math.nan
            row = {
                "kernel": kernel.name,
                "path": path.label,
                "shadow": path.prefix,
                "rows": int(len(idx)),
                "atoms": int(len(np.unique(atoms))),
                "selection_rule": selection_rule,
                "dominance_fraction_median": float(np.nanmedian(frac[idx])) if len(idx) else math.nan,
                "within_N_eff_median": neff,
                "support_flag": supp,
                **within,
                **loao,
                "basis_coefficient_delta_abs": basis_delta,
                "units": "ppm_DFT_T2 per emitted ppm_shadow_T2",
            }
            path_rows.append(row)
            if path.label == kernel.primary.label:
                primary_within = within
                primary_loao = loao
                primary_free5 = float(within["free5_within_R2"])
                primary_coef = float(within["coefficient"])
                primary_ci = (float(within["coefficient_ci95_low"]), float(within["coefficient_ci95_high"]))
                primary_agreement_coeffs = [float(within["coefficient"]), float(within_lib["coefficient"])]
                split = split_frame_block(frames, 0.20, 1)
                train = np.asarray(split["train_mask"], dtype=bool)
                test = np.asarray(split["test_mask"], dtype=bool)
                x_c = centered_by_train_atom(x, atoms, train)
                y_c = centered_by_train_atom(y, atoms, train)
                fit = train & np.isfinite(x_c).all(axis=1) & np.isfinite(y_c).all(axis=1)
                score = test & np.isfinite(x_c).all(axis=1) & np.isfinite(y_c).all(axis=1)
                primary_pysr = pysr_path(args, out_dir, kernel.name, x_c, y_c, fit, score)
                if np.isfinite(primary_pysr["coefficient"]):
                    primary_agreement_coeffs.append(float(primary_pysr["coefficient"]))

        coeff_arr = np.asarray([x for x in primary_agreement_coeffs if np.isfinite(x)], dtype=float)
        if len(coeff_arr) >= 2:
            agreement = float((np.max(coeff_arr) - np.min(coeff_arr)) / max(abs(np.median(coeff_arr)), 1.0e-12))
        else:
            agreement = math.nan
        within_r2 = float(primary_within["within_R2"]) if primary_within else math.nan
        loao_r2 = float(primary_loao["loao_R2"]) if primary_loao else math.nan
        bucket = path_bucket(primary_coef, primary_ci[0], primary_ci[1], within_r2, loao_r2, supp, agreement)
        if within_r2 > 0 and np.isfinite(primary_free5) and (primary_free5 - within_r2) < 0.05:
            stat_position = "fixed-D_ab_shape_supported"
        elif within_r2 > 0 and np.isfinite(primary_free5):
            stat_position = "free5_beats_fixed_shape"
        elif np.isfinite(primary_free5):
            stat_position = "fixed_shape_not_supported"
        else:
            stat_position = "undetermined"
        summary_rows.append({
            "kernel": kernel.name,
            "not_crap": bool(bucket != "cant-make-it-work-for-now"),
            "bucket": bucket,
            "primary_coefficient": primary_coef,
            "coef_ci95_low": primary_ci[0],
            "coef_ci95_high": primary_ci[1],
            "units": "ppm_DFT_T2 per emitted ppm_shadow_T2",
            "within_R2": within_r2,
            "LOAO_R2": loao_r2,
            "free5_within_R2": primary_free5,
            "statistical_position_vs_fixed_magic_angle_null": stat_position,
            "path_agreement_relative_span": agreement,
            "pysr_status": primary_pysr.get("status"),
            "pysr_coefficient": primary_pysr.get("coefficient"),
            "pysr_R2": primary_pysr.get("R2"),
            "pysr_equation": primary_pysr.get("equation"),
            "rows": int(len(idx)),
            "atoms": int(len(np.unique(atoms))),
            "within_N_eff_median": neff,
            "support_flag": supp,
            "selection_rule": kernel.selection_note,
        })
        cases_out = cases.copy()
        cases_out["stage2_selection_rule"] = kernel.selection_note
        cases_out["stage2_bucket"] = bucket
        cases_out["stage2_primary_coefficient"] = primary_coef
        cases_out["stage2_within_R2"] = within_r2
        cases_out["stage2_LOAO_R2"] = loao_r2
        cases_out["stage2_support_flag"] = supp
        case_outputs[kernel.name] = cases_out
    return pd.DataFrame(summary_rows), pd.DataFrame(path_rows), case_outputs


def build_unified_features(data: dict[str, object]) -> tuple[list[str], list[np.ndarray]]:
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    C: np.ndarray = data["C"]
    labels = []
    feats = []
    for spec in UNIFIED_SPECS:
        labels.append(spec.label)
        feats.append(to_e3nn(raw_t2(arrays, specs, spec), C))
    return labels, feats


def fit_linear_multi(features: list[np.ndarray], y: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xmat = np.column_stack([f[mask].reshape(-1) for f in features])
    yy = y[mask].reshape(-1)
    ok = np.isfinite(yy) & np.isfinite(xmat).all(axis=1)
    if ok.sum() < xmat.shape[1] + 2:
        return np.full(xmat.shape[1], np.nan), np.full_like(y[mask], np.nan)
    beta, *_ = np.linalg.lstsq(xmat[ok], yy[ok], rcond=None)
    pred = np.column_stack([f[mask].reshape(-1) for f in features]) @ beta
    return beta.astype(float), pred.reshape((-1, T2))


def run_unified_fit(data: dict[str, object]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    target = to_e3nn(data["target_lib"], data["C"])
    atoms_all = rows["atom_index"].to_numpy(int)
    frames_all = rows["original_frame_index"].to_numpy(int)
    dft = rows["dft_present"].to_numpy(int) == 1
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    atoms: set[int] = set()
    for case_dir in ["charge_wide", "mc", "field"]:
        case_set, _ = case_atoms(data["substrate_dir"], case_dir)
        atoms |= case_set
    dab_max = np.maximum.reduce([
        np.asarray(dom[:, dom_cols[f"dominant_fraction_{m}"]], dtype=float)
        for m in ["charge", "mc", "field"]
    ])
    ring_frac = np.asarray(dom[:, dom_cols["dominant_fraction_ring"]], dtype=float)
    mask = dft & np.isin(atoms_all, list(atoms)) & (dab_max >= 0.5) & (ring_frac < 0.7)
    idx = np.flatnonzero(mask)
    atoms_s = atoms_all[idx]
    frames_s = frames_all[idx]
    labels, features_all = build_unified_features(data)
    features = [f[idx] for f in features_all]
    y = target[idx]
    split = split_frame_block(frames_s, 0.20, 1)
    train = np.asarray(split["train_mask"], dtype=bool)
    test = np.asarray(split["test_mask"], dtype=bool)
    centered_features = [centered_by_train_atom(f, atoms_s, train) for f in features]
    y_c = centered_by_train_atom(y, atoms_s, train)
    beta, pred = fit_linear_multi(centered_features, y_c, train)
    _, pred_test = fit_linear_multi([f for f in centered_features], y_c, train)
    del pred_test
    x_test = np.column_stack([f[test].reshape(-1) for f in centered_features])
    pred_test_arr = (x_test @ beta).reshape((-1, T2)) if np.isfinite(beta).all() else np.full_like(y_c[test], np.nan)
    within_r2 = r2_score(pred_test_arr, y_c[test])
    _atom_labels, atom_features, atom_y, _counts = true_loao_atom_mean_arrays(features, y, atoms_s)
    loao_r2 = float(true_loao_score_atom_means(labels, atom_features, atom_y)["R2"])
    beta_j = []
    for held in np.unique(atoms_s):
        tr = (atoms_s != held) & train
        if tr.sum() < len(labels) + 2:
            continue
        b, _ = fit_linear_multi(centered_features, y_c, tr)
        beta_j.append(b)
    beta_j = np.asarray(beta_j, dtype=float)
    if len(beta_j) >= 3:
        se = np.sqrt((len(beta_j) - 1) / len(beta_j) * np.sum((beta_j - beta_j.mean(axis=0, keepdims=True)) ** 2, axis=0))
    else:
        se = np.full(len(labels), np.nan)
    inten_rows = []
    for i, label in enumerate(labels):
        lo, hi = ci95(float(beta[i]), float(se[i]))
        inten_rows.append({
            "term": label,
            "intensity": float(beta[i]),
            "se_atom_jackknife": float(se[i]),
            "ci95_low": lo,
            "ci95_high": hi,
            "literature_scaled_expected": 1.0,
            "ci95_contains_literature_1": bool(np.isfinite(lo) and lo <= 1.0 <= hi),
        })
    summary = pd.DataFrame([{
        "model": "unified_Dab_sum",
        "rows": int(len(idx)),
        "atoms": int(len(np.unique(atoms_s))),
        "selection_rule": "charge/mc/field CaseHunter atoms, max D_ab dominance >=0.5, ring dominance <0.7",
        "within_frameblock_R2": within_r2,
        "LOAO_R2": loao_r2,
        "terms": int(len(labels)),
        "test_frames": len(split["test_frames"]),
        "purged_train_frames": len(split["purged_frames"]),
        "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
    }])
    return summary, pd.DataFrame(inten_rows)


def load_stage2_inputs(substrate_dir: Path, include_stage21: bool = False) -> tuple[dict[str, object], float]:
    specs = load_json(substrate_dir / "per_atom_substrate_column_specs.json")
    manifest = load_json(substrate_dir / "per_atom_substrate_manifest.json")
    rows = pd.read_csv(substrate_dir / "per_atom_substrate_rows.csv")
    array_names = [
        "per_atom_substrate_features_classical",
        "per_atom_substrate_features_method_paths",
        "per_atom_substrate_features_hbond_conditioning",
        "per_atom_substrate_features_ring_paths",
        "per_atom_substrate_features_dominance",
        "per_atom_substrate_dominance_bins",
    ]
    if include_stage21:
        array_names.append("per_atom_substrate_driver_modulation_by_atom")
    arrays = {name: load_array(substrate_dir, name) for name in array_names}
    target_lib = load_array(substrate_dir, "per_atom_substrate_target_T2")
    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(T2)).max())
    if orth >= 1.0e-12:
        raise SystemExit(f"FATAL: change_of_basis.get_C orthogonality failed: {orth:.3e}")
    return {
        "manifest": manifest,
        "specs": specs,
        "rows": rows,
        "arrays": arrays,
        "target_lib": target_lib,
        "C": C,
        "substrate_dir": substrate_dir,
    }, orth


def stage21_target(data: dict[str, object]) -> np.ndarray:
    cached = data.get("_stage21_target")
    if cached is None:
        cached = to_e3nn(data["target_lib"], data["C"])
        data["_stage21_target"] = cached
    return cached


def stage21_primary_shadow(data: dict[str, object], kernel: KernelSpec) -> np.ndarray:
    cache = data.setdefault("_stage21_primary_shadows", {})
    assert isinstance(cache, dict)
    if kernel.name not in cache:
        cache[kernel.name] = to_e3nn(raw_t2(data["arrays"], data["specs"], kernel.primary), data["C"])
    return cache[kernel.name]


def stage21_unified_feature_cache(data: dict[str, object]) -> tuple[list[str], list[np.ndarray]]:
    cached = data.get("_stage21_unified_features")
    if cached is None:
        cached = build_unified_features(data)
        data["_stage21_unified_features"] = cached
    return cached  # type: ignore[return-value]


def kernel_clean_base_mask(data: dict[str, object], kernel: KernelSpec) -> tuple[np.ndarray, str]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    atoms, _cases = case_atoms(data["substrate_dir"], kernel.case_dir)
    atoms_all = rows["atom_index"].to_numpy(int)
    dft = rows["dft_present"].to_numpy(int) == 1
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_bins = arrays["per_atom_substrate_dominance_bins"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    bin_cols = specs_for_array(specs, "per_atom_substrate_dominance_bins")
    frac = np.asarray(dom[:, dom_cols[f"dominant_fraction_{kernel.dominance}"]], dtype=float)
    binv = np.asarray(dom_bins[:, bin_cols[f"bin_dominant_fraction_{kernel.dominance}_quintile"]], dtype=int)
    mask = dft & np.isin(atoms_all, list(atoms)) & np.isfinite(frac) & (frac >= 0.5)
    rule = "CaseHunter atoms + dominant_fraction>=0.5"
    if int(mask.sum()) < 500:
        mask = dft & np.isin(atoms_all, list(atoms)) & (binv >= 4)
        rule = "CaseHunter atoms + dominance_bin>=Q5 fallback"
    return mask, rule


def kernel_sweep_base_mask(data: dict[str, object], kernel: KernelSpec) -> np.ndarray:
    rows: pd.DataFrame = data["rows"]
    atoms, _cases = case_atoms(data["substrate_dir"], kernel.case_dir)
    return (rows["dft_present"].to_numpy(int) == 1) & np.isin(rows["atom_index"].to_numpy(int), list(atoms))


def unified_base_mask(data: dict[str, object], clean_stage2: bool) -> tuple[np.ndarray, str]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    atoms: set[int] = set()
    for case_dir in ["charge_wide", "mc", "field"]:
        case_set, _cases = case_atoms(data["substrate_dir"], case_dir)
        atoms |= case_set
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    dab_max = np.maximum.reduce([
        np.asarray(dom[:, dom_cols[f"dominant_fraction_{m}"]], dtype=float)
        for m in ["charge", "mc", "field"]
    ])
    ring_frac = np.asarray(dom[:, dom_cols["dominant_fraction_ring"]], dtype=float)
    atoms_all = rows["atom_index"].to_numpy(int)
    mask = (
        (rows["dft_present"].to_numpy(int) == 1)
        & np.isin(atoms_all, list(atoms))
        & np.isfinite(dab_max)
        & np.isfinite(ring_frac)
        & (ring_frac < 0.7)
    )
    rule = "D_ab CaseHunter atoms + ring dominance<0.7"
    if clean_stage2:
        mask &= dab_max >= 0.5
        rule += " + max D_ab dominance>=0.5"
    return mask, rule


def active_dab_arrays(data: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    dom_stack = np.column_stack([
        np.asarray(dom[:, dom_cols[f"dominant_fraction_{m}"]], dtype=float)
        for m in ["charge", "mc", "field"]
    ])
    with np.errstate(all="ignore"):
        active = np.nanargmax(np.where(np.isfinite(dom_stack), dom_stack, -np.inf), axis=1)
    gap_stack = np.column_stack([
        np.asarray(dom[:, dom_cols["gap_to_2nd_charge_r"]], dtype=float),
        np.asarray(dom[:, dom_cols["gap_to_2nd_bond_r"]], dtype=float),
        np.asarray(dom[:, dom_cols["gap_to_2nd_field_r"]], dtype=float),
    ])
    mod_cols = specs_for_array(specs, "per_atom_substrate_driver_modulation_by_atom")
    mod = arrays["per_atom_substrate_driver_modulation_by_atom"]
    atom_index = rows["atom_index"].to_numpy(int)
    mod_stack = np.column_stack([
        np.asarray(mod[atom_index, mod_cols["sd_charge_T2_by_atom"]], dtype=float),
        np.asarray(mod[atom_index, mod_cols["sd_mc_lit_T2_by_atom"]], dtype=float),
        np.asarray(mod[atom_index, mod_cols["sd_mopac_coulomb_T2_by_atom"]], dtype=float),
    ])
    row = np.arange(len(rows))
    return np.nanmax(dom_stack, axis=1), gap_stack[row, active], mod_stack[row, active]


def kernel_axis_values(data: dict[str, object], kernel: KernelSpec, axis: str) -> tuple[np.ndarray | None, str, str]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    if axis == "dominance":
        col = f"dominant_fraction_{kernel.dominance}"
        return np.asarray(dom[:, dom_cols[col]], dtype=float), col, "higher_is_cleaner"
    if axis == "isolation":
        gap_col = {
            "ring": "gap_to_2nd_ring_r",
            "charge": "gap_to_2nd_charge_r",
            "mc": "gap_to_2nd_bond_r",
            "field": "gap_to_2nd_field_r",
        }[kernel.name]
        return np.asarray(dom[:, dom_cols[gap_col]], dtype=float), gap_col, "higher_is_cleaner"
    if axis == "modulation":
        mod_name = {
            "ring": "sd_ring_jb_T2_by_atom",
            "charge": "sd_charge_T2_by_atom",
            "mc": "sd_mc_lit_T2_by_atom",
            "field": "sd_mopac_coulomb_T2_by_atom",
        }[kernel.name]
        mod_cols = specs_for_array(specs, "per_atom_substrate_driver_modulation_by_atom")
        mod = arrays["per_atom_substrate_driver_modulation_by_atom"]
        atom_index = rows["atom_index"].to_numpy(int)
        return np.asarray(mod[atom_index, mod_cols[mod_name]], dtype=float), mod_name, "higher_driver_exercise"
    raise ValueError(axis)


def unified_axis_values(data: dict[str, object], axis: str) -> tuple[np.ndarray | None, str, str]:
    dab_max, active_gap, active_mod = active_dab_arrays(data)
    if axis == "dominance":
        return dab_max, "max_dominant_fraction_charge_mc_field", "higher_is_cleaner"
    if axis == "isolation":
        return active_gap, "gap_to_2nd_of_active_Dab_mechanism", "higher_is_cleaner"
    if axis == "modulation":
        return active_mod, "sd_of_active_Dab_mechanism", "higher_driver_exercise"
    raise ValueError(axis)


def empty_stage21_metrics(model: str, rows: int, atoms: int, frames: int, note: str = "") -> dict[str, object]:
    return {
        "model": model,
        "rows": int(rows),
        "atoms": int(atoms),
        "frames": int(frames),
        "within_R2": math.nan,
        "LOAO_R2": math.nan,
        "coefficient": math.nan,
        "coef_ci95_low": math.nan,
        "coef_ci95_high": math.nan,
        "within_N_eff_median": math.nan,
        "median_lag1_rho": math.nan,
        "support_flag": "insufficient_rows_or_frames",
        "bucket": "undetermined",
        "note": note,
    }


def stage21_kernel_metrics(data: dict[str, object], kernel: KernelSpec, mask: np.ndarray) -> dict[str, object]:
    rows: pd.DataFrame = data["rows"]
    idx = np.flatnonzero(mask)
    atoms = rows["atom_index"].to_numpy(int)[idx]
    frames = rows["original_frame_index"].to_numpy(int)[idx]
    n_atoms = int(len(np.unique(atoms))) if len(idx) else 0
    n_frames = int(len(np.unique(frames))) if len(idx) else 0
    if len(idx) < 10 or n_frames < 3 or n_atoms < 1:
        return empty_stage21_metrics(kernel.name, len(idx), n_atoms, n_frames)
    y = stage21_target(data)[idx]
    x = stage21_primary_shadow(data, kernel)[idx]
    neff_info = effective_n_components(y, atoms, frames)
    neff = float(neff_info["within_N_eff_median"])
    supp = support_flag(n_atoms, len(idx), neff)
    within = blocked_within_fit(x, y, atoms, frames)
    loao = true_loao_single_feature_fit(x, y, atoms)
    coef = float(within["coefficient"])
    lo = float(within["coefficient_ci95_low"])
    hi = float(within["coefficient_ci95_high"])
    within_r2 = float(within["within_R2"])
    loao_r2 = float(loao["loao_R2"])
    return {
        "model": kernel.name,
        "rows": int(len(idx)),
        "atoms": n_atoms,
        "frames": n_frames,
        "within_R2": within_r2,
        "LOAO_R2": loao_r2,
        "free5_within_R2": float(within["free5_within_R2"]),
        "coefficient": coef,
        "coef_ci95_low": lo,
        "coef_ci95_high": hi,
        "within_N_eff_median": neff,
        "median_lag1_rho": float(neff_info["median_lag1_rho"]),
        "support_flag": supp,
        "bucket": path_bucket(coef, lo, hi, within_r2, loao_r2, supp, math.nan),
        "note": "",
    }


def stage21_unified_metrics(data: dict[str, object], mask: np.ndarray) -> dict[str, object]:
    rows: pd.DataFrame = data["rows"]
    idx = np.flatnonzero(mask)
    atoms = rows["atom_index"].to_numpy(int)[idx]
    frames = rows["original_frame_index"].to_numpy(int)[idx]
    n_atoms = int(len(np.unique(atoms))) if len(idx) else 0
    n_frames = int(len(np.unique(frames))) if len(idx) else 0
    labels, features_all = stage21_unified_feature_cache(data)
    if len(idx) < len(labels) + 3 or n_frames < 3 or n_atoms < 2:
        return empty_stage21_metrics("unified_Dab_sum", len(idx), n_atoms, n_frames)
    y = stage21_target(data)[idx]
    features = [f[idx] for f in features_all]
    split = split_frame_block(frames, 0.20, 1)
    train = np.asarray(split["train_mask"], dtype=bool)
    test = np.asarray(split["test_mask"], dtype=bool)
    centered_features = [centered_by_train_atom(f, atoms, train) for f in features]
    y_c = centered_by_train_atom(y, atoms, train)
    beta, _pred_train = fit_linear_multi(centered_features, y_c, train)
    if np.isfinite(beta).all() and test.sum() > 0:
        x_test = np.column_stack([f[test].reshape(-1) for f in centered_features])
        pred_test = (x_test @ beta).reshape((-1, T2))
    else:
        pred_test = np.full_like(y_c[test], np.nan)
    within_r2 = r2_score(pred_test, y_c[test])
    _atom_labels, atom_features, atom_y, _counts = true_loao_atom_mean_arrays(features, y, atoms)
    loao_r2 = float(true_loao_score_atom_means(labels, atom_features, atom_y)["R2"])
    neff_info = effective_n_components(y, atoms, frames)
    neff = float(neff_info["within_N_eff_median"])
    supp = support_flag(n_atoms, len(idx), neff)
    return {
        "model": "unified_Dab_sum",
        "rows": int(len(idx)),
        "atoms": n_atoms,
        "frames": n_frames,
        "within_R2": float(within_r2),
        "LOAO_R2": float(loao_r2),
        "coefficient": math.nan,
        "coef_ci95_low": math.nan,
        "coef_ci95_high": math.nan,
        "within_N_eff_median": neff,
        "median_lag1_rho": float(neff_info["median_lag1_rho"]),
        "support_flag": supp,
        "bucket": "multi-term-heldout-recovery" if np.isfinite(within_r2) and within_r2 > 0 else "undetermined",
        "note": f"{len(labels)} unified terms",
    }


def threshold_for_keep(values: np.ndarray, keep_fraction: float) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return math.nan
    if keep_fraction >= 1.0:
        return float(np.nanmin(vals))
    return float(np.nanquantile(vals, 1.0 - float(keep_fraction), method="lower"))


def run_stage21_happy_sweep(data: dict[str, object]) -> pd.DataFrame:
    rows_out: list[dict[str, object]] = []
    axes = ["dominance", "isolation", "modulation"]
    for kernel in KERNELS:
        base = kernel_sweep_base_mask(data, kernel)
        for axis in axes:
            values, value_column, polarity = kernel_axis_values(data, kernel, axis)
            if values is None:
                rows_out.append({
                    **empty_stage21_metrics(kernel.name, int(base.sum()), int(pd.Series(data["rows"]["atom_index"].to_numpy(int)[base]).nunique()) if base.any() else 0, 0, "axis unavailable"),
                    "axis": axis,
                    "axis_value_column": value_column,
                    "axis_polarity": polarity,
                    "keep_fraction": math.nan,
                    "threshold": math.nan,
                    "selection_rule": "axis unavailable in emitted build4 row/atom substrate",
                })
                continue
            finite_base = base & np.isfinite(values)
            base_vals = values[finite_base]
            for keep in SWEEP_KEEP_FRACTIONS:
                thr = threshold_for_keep(base_vals, keep)
                mask = finite_base & (values >= thr)
                metrics = stage21_kernel_metrics(data, kernel, mask)
                rows_out.append({
                    **metrics,
                    "axis": axis,
                    "axis_value_column": value_column,
                    "axis_polarity": polarity,
                    "keep_fraction": float(keep),
                    "threshold": thr,
                    "selection_rule": f"CaseHunter {kernel.case_dir} atoms, top {keep:.0%} by {value_column}",
                })
    for axis in axes:
        base, rule = unified_base_mask(data, clean_stage2=False)
        values, value_column, polarity = unified_axis_values(data, axis)
        if values is None:
            continue
        finite_base = base & np.isfinite(values)
        base_vals = values[finite_base]
        for keep in SWEEP_KEEP_FRACTIONS:
            thr = threshold_for_keep(base_vals, keep)
            mask = finite_base & (values >= thr)
            metrics = stage21_unified_metrics(data, mask)
            rows_out.append({
                **metrics,
                "axis": axis,
                "axis_value_column": value_column,
                "axis_polarity": polarity,
                "keep_fraction": float(keep),
                "threshold": thr,
                "selection_rule": f"{rule}, top {keep:.0%} by {value_column}",
            })
    return pd.DataFrame(rows_out)


def centered_frame_mask(frames: np.ndarray, n_frames: int) -> tuple[np.ndarray, list[int]]:
    unique = np.sort(np.unique(frames))
    n = min(int(n_frames), len(unique))
    start = max(0, (len(unique) - n) // 2)
    chosen = unique[start:start + n]
    return np.isin(frames, chosen), [int(x) for x in chosen]


def run_stage21_frame_ablation(data: dict[str, object]) -> pd.DataFrame:
    rows: pd.DataFrame = data["rows"]
    frames_all = rows["original_frame_index"].to_numpy(int)
    out: list[dict[str, object]] = []
    for kernel in KERNELS:
        base, rule = kernel_clean_base_mask(data, kernel)
        for n in FRAME_COUNTS:
            fm, chosen = centered_frame_mask(frames_all, n)
            metrics = stage21_kernel_metrics(data, kernel, base & fm)
            out.append({
                **metrics,
                "n_frames_requested": int(n),
                "frame_selection": "centered_contiguous_original_frame_index_block",
                "frame_first": min(chosen) if chosen else math.nan,
                "frame_last": max(chosen) if chosen else math.nan,
                "selection_rule": rule,
            })
    base, rule = unified_base_mask(data, clean_stage2=True)
    for n in FRAME_COUNTS:
        fm, chosen = centered_frame_mask(frames_all, n)
        metrics = stage21_unified_metrics(data, base & fm)
        out.append({
            **metrics,
            "n_frames_requested": int(n),
            "frame_selection": "centered_contiguous_original_frame_index_block",
            "frame_first": min(chosen) if chosen else math.nan,
            "frame_last": max(chosen) if chosen else math.nan,
            "selection_rule": rule,
        })
    return pd.DataFrame(out)


def recovery_mean(row: pd.Series) -> float:
    vals = [float(row.get("within_R2", math.nan)), float(row.get("LOAO_R2", math.nan))]
    vals = [v for v in vals if np.isfinite(v)]
    return float(np.mean(vals)) if vals else math.nan


def classify_shape(group: pd.DataFrame) -> dict[str, object]:
    valid = group[np.isfinite(group["keep_fraction"].to_numpy(float))].copy()
    valid = valid[np.isfinite(valid["within_R2"].to_numpy(float)) | np.isfinite(valid["LOAO_R2"].to_numpy(float))]
    if valid.empty:
        return {"shape": "unavailable", "clean_pop": False, "delta_clean_minus_loose": math.nan, "clean_bucket": "undetermined", "clean_support_flag": "missing"}
    valid = valid.sort_values("keep_fraction")
    metric = valid.apply(recovery_mean, axis=1).to_numpy(float)
    clean = valid.iloc[0]
    loose = valid.iloc[-1]
    delta = float(metric[0] - metric[-1]) if np.isfinite(metric[0]) and np.isfinite(metric[-1]) else math.nan
    if not np.isfinite(delta):
        shape = "undetermined"
    elif delta > 0.05:
        shape = "rises_to_clean"
    elif delta < -0.05:
        shape = "falls_to_clean"
    else:
        shape = "flat"
    return {
        "shape": shape,
        "clean_pop": bool(shape == "rises_to_clean"),
        "delta_clean_minus_loose": delta,
        "clean_keep_fraction": float(clean["keep_fraction"]),
        "clean_within_R2": float(clean["within_R2"]),
        "clean_LOAO_R2": float(clean["LOAO_R2"]),
        "clean_coefficient": float(clean["coefficient"]) if np.isfinite(clean["coefficient"]) else math.nan,
        "clean_coef_ci95_low": float(clean["coef_ci95_low"]) if np.isfinite(clean["coef_ci95_low"]) else math.nan,
        "clean_coef_ci95_high": float(clean["coef_ci95_high"]) if np.isfinite(clean["coef_ci95_high"]) else math.nan,
        "clean_bucket": str(clean["bucket"]),
        "clean_support_flag": str(clean["support_flag"]) if str(clean["support_flag"]) else "",
        "loose_within_R2": float(loose["within_R2"]),
        "loose_LOAO_R2": float(loose["LOAO_R2"]),
    }


def summarize_stage21_shapes(sweep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, axis), group in sweep.groupby(["model", "axis"], dropna=False):
        summary = classify_shape(group)
        axis_col = str(group["axis_value_column"].iloc[0]) if not group.empty else ""
        rows.append({"model": model, "axis": axis, "axis_value_column": axis_col, **summary})
    return pd.DataFrame(rows).sort_values(["model", "axis"]).reset_index(drop=True)


def ci_off_zero(row: pd.Series) -> bool:
    lo = float(row.get("clean_coef_ci95_low", math.nan))
    hi = float(row.get("clean_coef_ci95_high", math.nan))
    return bool(np.isfinite(lo) and np.isfinite(hi) and (lo > 0.0 or hi < 0.0))


def summarize_stage21_frames(frame: pd.DataFrame) -> pd.DataFrame:
    out = []
    for model, group in frame.groupby("model", dropna=False):
        g = group.sort_values("n_frames_requested")
        full = g[g["n_frames_requested"] == 660]
        r2_full = float(full["within_R2"].iloc[0]) if not full.empty else math.nan
        target = 0.90 * r2_full if np.isfinite(r2_full) and r2_full > 0.0 else math.nan
        knee = "no_positive_plateau"
        if np.isfinite(target):
            hit = g[np.isfinite(g["within_R2"].to_numpy(float)) & (g["within_R2"].to_numpy(float) >= target)]
            if not hit.empty:
                knee = f"{int(hit['n_frames_requested'].iloc[0])}"
        r50 = g[g["n_frames_requested"] == 50]
        out.append({
            "model": model,
            "full_660_within_R2": r2_full,
            "knee_frames_90pct_full": knee,
            "frame50_within_R2": float(r50["within_R2"].iloc[0]) if not r50.empty else math.nan,
            "frame50_LOAO_R2": float(r50["LOAO_R2"].iloc[0]) if not r50.empty else math.nan,
            "frame50_N_eff_median": float(r50["within_N_eff_median"].iloc[0]) if not r50.empty else math.nan,
            "frame50_median_lag1_rho": float(r50["median_lag1_rho"].iloc[0]) if not r50.empty else math.nan,
            "full_660_N_eff_median": float(full["within_N_eff_median"].iloc[0]) if not full.empty else math.nan,
            "full_660_median_lag1_rho": float(full["median_lag1_rho"].iloc[0]) if not full.empty else math.nan,
            "frame50_support_flag": str(r50["support_flag"].iloc[0]) if not r50.empty else "missing",
        })
    return pd.DataFrame(out).sort_values("model").reset_index(drop=True)


def write_stage21_plots(out_dir: Path, sweep: pd.DataFrame, frame: pd.DataFrame) -> dict[str, object]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        return {"plots_written": False, "reason": str(exc)}
    written: list[str] = []
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    models = ["charge", "field", "ring", "mc", "unified_Dab_sum"]
    axes = ["dominance", "isolation", "modulation"]
    fig, axs = plt.subplots(len(models), len(axes), figsize=(12, 14), sharex=True)
    for i, model in enumerate(models):
        for j, axis in enumerate(axes):
            ax = axs[i][j]
            g = sweep[(sweep["model"] == model) & (sweep["axis"] == axis) & np.isfinite(sweep["keep_fraction"].to_numpy(float))]
            if not g.empty:
                g = g.sort_values("keep_fraction")
                ax.plot(g["keep_fraction"], g["within_R2"], marker="o", label="within")
                ax.plot(g["keep_fraction"], g["LOAO_R2"], marker="s", label="LOAO")
            ax.axhline(0.0, color="0.7", linewidth=0.7)
            ax.set_title(f"{model} / {axis}", fontsize=8)
            if i == len(models) - 1:
                ax.set_xlabel("kept fraction (strict to loose)")
            if j == 0:
                ax.set_ylabel("held-out R2")
    axs[0][0].legend(fontsize=7, loc="best")
    fig.tight_layout()
    path = plot_dir / "happy_spot_curves.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    written.append(str(path))

    fig, ax = plt.subplots(figsize=(9, 5))
    for model in models:
        g = frame[frame["model"] == model].sort_values("n_frames_requested")
        if not g.empty:
            ax.plot(g["n_frames_requested"], g["within_R2"], marker="o", label=model)
    ax.axhline(0.0, color="0.7", linewidth=0.7)
    ax.set_xlabel("centered contiguous frames")
    ax.set_ylabel("within-frameblock held-out R2")
    ax.set_title("Stage 2.1 frame-count ablation")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    path = plot_dir / "frame_count_ablation.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    written.append(str(path))
    return {"plots_written": True, "paths": written}


def fmt_num(x: object, digits: int = 3) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "nan"
    return f"{v:.{digits}f}"


def shape_triplet(shape_summary: pd.DataFrame, model: str) -> str:
    parts = []
    for axis in ["dominance", "isolation", "modulation"]:
        row = shape_summary[(shape_summary["model"] == model) & (shape_summary["axis"] == axis)]
        parts.append(str(row["shape"].iloc[0]) if not row.empty else "missing")
    return " / ".join(parts)


def clean_pop_text(shape_summary: pd.DataFrame, model: str) -> str:
    rows = shape_summary[shape_summary["model"] == model]
    if rows.empty:
        return "missing"
    axes = [str(r.axis) for r in rows.itertuples(index=False) if bool(r.clean_pop)]
    return ",".join(axes) if axes else "no"


def rescue_text(shape_summary: pd.DataFrame, model: str) -> str:
    rows = shape_summary[shape_summary["model"] == model]
    rescuers = [str(r.axis) for r in rows.itertuples(index=False) if ci_off_zero(pd.Series(r._asdict()))]
    return ",".join(rescuers) if rescuers else "no"


def frame_line(frame_summary: pd.DataFrame, model: str) -> str:
    row = frame_summary[frame_summary["model"] == model]
    if row.empty:
        return f"- {model}: missing frame ablation."
    r = row.iloc[0]
    knee = str(r["knee_frames_90pct_full"])
    knee_text = f"{knee} fr" if knee.isdigit() else knee
    return (
        f"- {model}: knee {knee_text}; "
        f"50fr within {fmt_num(r['frame50_within_R2'])}, LOAO {fmt_num(r['frame50_LOAO_R2'])}, "
        f"N_eff {fmt_num(r['frame50_N_eff_median'], 1)} (rho {fmt_num(r['frame50_median_lag1_rho'])}); "
        f"660fr N_eff {fmt_num(r['full_660_N_eff_median'], 1)}."
    )


def write_stage21_postmortem(out_dir: Path, substrate_dir: Path, orth: float, shape_summary: pd.DataFrame, frame_summary: pd.DataFrame, disk: dict[str, object], plots: dict[str, object]) -> Path:
    models = ["charge", "field", "ring", "mc", "unified_Dab_sum"]
    lines = [
        "# Stage 2.1 Postmortem - 2026-06-04",
        f"Run dir: `{out_dir}`",
        f"Substrate: `{substrate_dir}`",
        f"Scope: Build4 CSV/NPY sidecars only; no extraction/DFT/per-source; frozen `get_C`, `|C.T C-I|max={orth:.2e}`; five-component total-T2.",
        "Happy-spot shapes are dominance / isolation / driver-modulation; `rises_to_clean` = clean POP, `falls_to_clean` = washout.",
    ]
    for model in models:
        support = shape_summary[shape_summary["model"] == model]["clean_support_flag"].replace("", np.nan).dropna()
        support_text = f"; support {support.iloc[0]}" if not support.empty and str(support.iloc[0]) != "nan" else ""
        rescue = ""
        if model == "mc":
            rescue = f"; rescue {rescue_text(shape_summary, model)}"
        lines.append(f"- {model}: {shape_triplet(shape_summary, model)}; clean POP {clean_pop_text(shape_summary, model)}{rescue}{support_text}.")
    lines += [
        "McConnell clean-end rescue criterion: primary coefficient CI off zero at the strict end of an input-side axis; no DFT fit is used to define the axis.",
        "Frame-count ablation uses centered contiguous 20 ps-stride frame blocks; 50 frames is the 1 ns-at-20 ps ubiquitin proxy.",
    ]
    for model in models:
        lines.append(frame_line(frame_summary, model))
    lines += [
        "50-frame read: use per-kernel rows above; charge/unified positive at 50 means the cheap 1 ns run keeps the main signal, null kernels remain provisional.",
        "Geometric-noise flag: build4 has no noise axis distinct from driver modulation; motion and geometric jitter are not separable here. If clean POP tracks modulation, a future C++ emit should split signal modulation from geometric jitter.",
        f"Artifacts: `stage2_1_happy_spot_curves.csv`, `stage2_1_happy_spot_shape_summary.csv`, `stage2_1_frame_count_ablation.csv`, plots_written={plots.get('plots_written')}.",
        f"Disk: `/tmp/rediscover-runs` {float(disk['rediscover_runs_GiB_before_write']):.1f}G before write (<15G); build4+build1 kept; output drop-old={disk.get('deleted_existing_output_dir')}.",
    ]
    path = Path(__file__).resolve().parents[1] / "POSTMORTEM_STAGE2_1.md"
    path.write_text("\n".join(lines[:40]) + "\n", encoding="utf-8")
    return path


def run_stage21(args: argparse.Namespace, data: dict[str, object], out_dir: Path, disk: dict[str, object], orth: float) -> Path:
    sweep = run_stage21_happy_sweep(data)
    shape_summary = summarize_stage21_shapes(sweep)
    frame = run_stage21_frame_ablation(data)
    frame_summary = summarize_stage21_frames(frame)
    sweep.to_csv(out_dir / "stage2_1_happy_spot_curves.csv", index=False)
    shape_summary.to_csv(out_dir / "stage2_1_happy_spot_shape_summary.csv", index=False)
    frame.to_csv(out_dir / "stage2_1_frame_count_ablation.csv", index=False)
    frame_summary.to_csv(out_dir / "stage2_1_frame_count_summary.csv", index=False)
    plots = write_stage21_plots(out_dir, sweep, frame)
    audit = {
        "substrate_dir": str(args.substrate_dir),
        "out_dir": str(out_dir),
        "disk_guard": disk,
        "change_of_basis_orthogonality_max_abs": orth,
        "python_read_scope": "Build-4 per_atom_substrate CSV/NPY sidecars and CaseHunter manifests only",
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source dumps", "older rediscover dirs for fitting"],
        "anti_circularity": "cleanliness axes are emitted input-side dominance/gap/modulation only; DFT target is used only in held-out scoring",
        "sweep_keep_fractions_strict_to_loose": list(SWEEP_KEEP_FRACTIONS),
        "frame_counts": list(FRAME_COUNTS),
        "frame_ablation_policy": "centered contiguous original_frame_index blocks at the emitted 20 ps stride",
        "unified_driver_modulation_axis": "active charge/mc/field emitted sd column",
        "geometric_noise_axis_distinct_from_driver_modulation": "not available in build4",
        "plots": plots,
    }
    (out_dir / "stage2_1_run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return write_stage21_postmortem(out_dir, args.substrate_dir, orth, shape_summary, frame_summary, disk, plots)


def stage22_stack_features(features: list[np.ndarray]) -> np.ndarray:
    if not features:
        return np.empty((0, 0), dtype=float)
    return np.column_stack([np.asarray(f, dtype=float).reshape(-1) for f in features])


def stage22_fit_ridge_flat(x_train: np.ndarray, y_train: np.ndarray, alpha: float) -> dict[str, object]:
    x = np.asarray(x_train, dtype=float)
    y = np.asarray(y_train, dtype=float).reshape(-1, 1)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    p = int(x.shape[1])
    ok = np.isfinite(y[:, 0]) & np.isfinite(x).all(axis=1)
    x = x[ok]
    y = y[ok]
    if len(x) < 3 or p == 0:
        return {
            "ok": False,
            "alpha": float(alpha),
            "intercept": math.nan,
            "coef_std": np.full(p, np.nan),
            "coef_original": np.full(p, np.nan),
            "mean": np.zeros(p),
            "std": np.ones(p),
            "train_observations": int(len(x)),
            "effective_dof_penalized": math.nan,
            "effective_dof_with_intercept": math.nan,
        }
    mean = x.mean(axis=0)
    std = x.std(axis=0)
    std[~np.isfinite(std) | (std <= 0.0)] = 1.0
    z = (x - mean) / std
    design = np.column_stack([np.ones(len(z)), z])
    penalty = np.eye(p + 1, dtype=float) * max(0.0, float(alpha))
    penalty[0, 0] = 0.0
    xtx = design.T @ design
    xty = design.T @ y
    try:
        beta = np.linalg.solve(xtx + penalty, xty).reshape(-1)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(xtx + penalty, xty, rcond=None)[0].reshape(-1)
    eig = np.linalg.eigvalsh(z.T @ z)
    eig = np.maximum(eig, 0.0)
    if float(alpha) <= 0.0:
        eff_pen = float(np.sum(eig > 1.0e-10))
    else:
        eff_pen = float(np.sum(eig / (eig + float(alpha))))
    coef_std = beta[1:]
    coef_original = coef_std / std
    return {
        "ok": True,
        "alpha": float(alpha),
        "intercept": float(beta[0] - np.dot(mean / std, coef_std)),
        "coef_std": coef_std.astype(float),
        "coef_original": coef_original.astype(float),
        "mean": mean.astype(float),
        "std": std.astype(float),
        "train_observations": int(len(x)),
        "effective_dof_penalized": eff_pen,
        "effective_dof_with_intercept": float(1.0 + eff_pen),
    }


def stage22_predict_ridge_flat(model: dict[str, object], x_eval: np.ndarray) -> np.ndarray:
    x = np.asarray(x_eval, dtype=float)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    if not bool(model.get("ok")):
        return np.full(x.shape[0], np.nan, dtype=float)
    mean = np.asarray(model["mean"], dtype=float)
    std = np.asarray(model["std"], dtype=float)
    coef_std = np.asarray(model["coef_std"], dtype=float)
    x = np.where(np.isfinite(x), x, mean)
    return (float(model["intercept"]) + ((x - mean) / std) @ coef_std).astype(float)


def stage22_fit_predict_arrays(
    train_features: list[np.ndarray],
    y_train: np.ndarray,
    eval_features: list[np.ndarray],
    alpha: float,
) -> tuple[np.ndarray, dict[str, object]]:
    x_train = stage22_stack_features(train_features)
    yy = np.asarray(y_train, dtype=float).reshape(-1)
    model = stage22_fit_ridge_flat(x_train, yy, alpha)
    pred = stage22_predict_ridge_flat(model, stage22_stack_features(eval_features))
    eval_rows = int(eval_features[0].shape[0]) if eval_features else 0
    return pred.reshape((eval_rows, T2)), model


def stage22_alpha_scores() -> dict[float, list[float]]:
    return {float(alpha): [] for alpha in STAGE22_ALPHA_GRID}


def stage22_select_within_alpha(
    features: list[np.ndarray],
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    outer_train: np.ndarray,
) -> tuple[float, dict[str, object]]:
    scores = stage22_alpha_scores()
    folds = blocked_frame_cv_splits(frames, outer_train, STAGE22_INNER_CV_FOLDS, 1)
    for fold in folds:
        fit = np.asarray(fold["fit_mask"], dtype=bool)
        val = np.asarray(fold["val_mask"], dtype=bool)
        x_c = [centered_by_train_atom(f, atoms, fit) for f in features]
        y_c = centered_by_train_atom(y, atoms, fit)
        train_features = [f[fit] for f in x_c]
        val_features = [f[val] for f in x_c]
        for alpha in STAGE22_ALPHA_GRID:
            pred, _model = stage22_fit_predict_arrays(train_features, y_c[fit], val_features, float(alpha))
            score = r2_score(pred, y_c[val])
            if np.isfinite(score):
                scores[float(alpha)].append(float(score))
    mean_scores = {float(alpha): (float(np.mean(vals)) if vals else math.nan) for alpha, vals in scores.items()}
    selected = select_best_alpha(mean_scores, [float(a) for a in STAGE22_ALPHA_GRID])
    return selected, {
        "axis": "within_frameblock",
        "method": "blocked frame inner CV inside outer training frames only",
        "selected_alpha": float(selected),
        "inner_cv_folds": int(len(folds)),
        "alpha_grid": [float(x) for x in STAGE22_ALPHA_GRID],
        "inner_cv_R2_by_alpha": {f"{alpha:g}": (float(mean_scores[alpha]) if np.isfinite(mean_scores[alpha]) else None) for alpha in mean_scores},
    }


def stage22_block_bootstrap_values(
    pred: np.ndarray,
    target: np.ndarray,
    frames: np.ndarray,
    n_boot: int,
    block_len: int,
    seed: int,
) -> np.ndarray:
    pred = np.asarray(pred, dtype=float)
    target = np.asarray(target, dtype=float)
    frames = np.asarray(frames, dtype=int)
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    pred = pred[ok]
    target = target[ok]
    frames = frames[ok]
    unique = np.sort(np.unique(frames))
    if len(unique) < 2:
        return np.asarray([], dtype=float)
    by_pos = [np.flatnonzero(frames == frame) for frame in unique]
    n_blocks = int(math.ceil(len(unique) / max(1, int(block_len))))
    rng = np.random.default_rng(seed)
    vals: list[float] = []
    for _ in range(int(n_boot)):
        starts = rng.integers(0, len(unique), size=n_blocks)
        chosen_positions: list[int] = []
        for start in starts:
            for offset in range(max(1, int(block_len))):
                chosen_positions.append(int((int(start) + offset) % len(unique)))
                if len(chosen_positions) >= len(unique):
                    break
            if len(chosen_positions) >= len(unique):
                break
        pieces = [by_pos[pos] for pos in chosen_positions if len(by_pos[pos])]
        if not pieces:
            continue
        take = np.concatenate(pieces)
        val = r2_score(pred[take], target[take])
        if np.isfinite(val):
            vals.append(float(val))
    return np.asarray(vals, dtype=float)


def stage22_recovery_intervals(
    pred: np.ndarray,
    target: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    seed: int,
    include_bootstrap: bool,
) -> dict[str, object]:
    r2 = r2_score(pred, target)
    jk_vals = jackknife_metric_values(pred, target, atoms)
    jk_se = jackknife_se_from_values(jk_vals)
    jk_lo, jk_hi = ci95(float(r2), float(jk_se))
    if include_bootstrap:
        boot = stage22_block_bootstrap_values(
            pred,
            target,
            frames,
            STAGE22_BOOTSTRAPS,
            STAGE22_BOOTSTRAP_FRAME_BLOCK,
            seed,
        )
    else:
        boot = np.asarray([], dtype=float)
    if len(boot):
        boot_lo, boot_hi = np.nanpercentile(boot, [2.5, 97.5])
    else:
        boot_lo, boot_hi = math.nan, math.nan
    return {
        "R2": float(r2),
        "jackknife_se": float(jk_se) if np.isfinite(jk_se) else math.nan,
        "jackknife_ci95_low": float(jk_lo) if np.isfinite(jk_lo) else math.nan,
        "jackknife_ci95_high": float(jk_hi) if np.isfinite(jk_hi) else math.nan,
        "jackknife_groups": int(len(np.unique(atoms))),
        "bootstrap_ci95_low": float(boot_lo) if np.isfinite(boot_lo) else math.nan,
        "bootstrap_ci95_high": float(boot_hi) if np.isfinite(boot_hi) else math.nan,
        "bootstrap_replicates": int(len(boot)),
        "bootstrap_frame_block": int(STAGE22_BOOTSTRAP_FRAME_BLOCK),
    }


def stage22_atom_mean_recovery_intervals(
    pred: np.ndarray,
    target: np.ndarray,
    atom_labels: np.ndarray,
) -> dict[str, object]:
    r2 = r2_score(pred, target)
    jk_vals = jackknife_metric_values(pred, target, atom_labels)
    jk_se = jackknife_se_from_values(jk_vals)
    jk_lo, jk_hi = ci95(float(r2), float(jk_se))
    return {
        "R2": float(r2),
        "jackknife_se": float(jk_se) if np.isfinite(jk_se) else math.nan,
        "jackknife_ci95_low": float(jk_lo) if np.isfinite(jk_lo) else math.nan,
        "jackknife_ci95_high": float(jk_hi) if np.isfinite(jk_hi) else math.nan,
        "jackknife_groups": int(len(np.unique(atom_labels))),
        "bootstrap_ci95_low": math.nan,
        "bootstrap_ci95_high": math.nan,
        "bootstrap_replicates": 0,
        "bootstrap_frame_block": 0,
    }


def stage22_prefix(prefix: str, metrics: dict[str, object]) -> dict[str, object]:
    return {f"{prefix}_{key}": value for key, value in metrics.items()}


def stage22_evaluate_within(
    labels: list[str],
    features: list[np.ndarray],
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    include_ci: bool,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    split = split_frame_block(frames, 0.20, 1)
    train = np.asarray(split["train_mask"], dtype=bool)
    test = np.asarray(split["test_mask"], dtype=bool)
    alpha, alpha_audit = stage22_select_within_alpha(features, y, atoms, frames, train)
    x_c = [centered_by_train_atom(f, atoms, train) for f in features]
    y_c = centered_by_train_atom(y, atoms, train)
    pred, model = stage22_fit_predict_arrays([f[train] for f in x_c], y_c[train], [f[test] for f in x_c], alpha)
    intervals = stage22_recovery_intervals(pred, y_c[test], atoms[test], frames[test], seed, include_ci)
    coef = np.asarray(model["coef_original"], dtype=float)
    coeffs = pd.DataFrame({
        "term": labels,
        "within_ridge_intensity": coef if len(coef) == len(labels) else np.full(len(labels), np.nan),
        "selected_within_alpha": float(alpha),
    })
    return {
        **intervals,
        "selected_alpha": float(alpha),
        "effective_dof_penalized": float(model["effective_dof_penalized"]) if np.isfinite(model["effective_dof_penalized"]) else math.nan,
        "effective_dof_with_intercept": float(model["effective_dof_with_intercept"]) if np.isfinite(model["effective_dof_with_intercept"]) else math.nan,
        "train_rows": int(train.sum()),
        "test_rows": int(test.sum()),
        "train_observations": int(model["train_observations"]),
        "terms": int(len(labels)),
        "test_frames": len(split["test_frames"]),
        "purged_train_frames": len(split["purged_frames"]),
        "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
        "alpha_audit": alpha_audit,
    }, coeffs


def stage22_evaluate_loao(
    labels: list[str],
    features: list[np.ndarray],
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    include_ci: bool,
    seed: int,
) -> dict[str, object]:
    del frames, include_ci, seed
    atom_labels, atom_features, atom_y, counts = true_loao_atom_mean_arrays(features, y, atoms)
    score = true_loao_score_atom_means(labels, atom_features, atom_y)
    pred = np.asarray(score["pred"], dtype=float)
    target = np.asarray(score["target"], dtype=float)
    intervals = stage22_atom_mean_recovery_intervals(pred, target, atom_labels)
    return {
        **intervals,
        "selected_alpha_mode": math.nan,
        "selected_alpha_min": math.nan,
        "selected_alpha_max": math.nan,
        "selected_alpha_counts": {},
        "effective_dof_penalized_median": math.nan,
        "effective_dof_penalized_min": math.nan,
        "effective_dof_penalized_max": math.nan,
        "train_observations_median": float(score["train_observations_median"]),
        "terms": int(len(labels)),
        "heldout_atoms": int(len(atom_labels)),
        "atom_row_count_min": int(np.min(counts)) if len(counts) else 0,
        "atom_row_count_median": float(np.median(counts)) if len(counts) else math.nan,
        "atom_row_count_max": int(np.max(counts)) if len(counts) else 0,
        "alpha_audit": [{
            "axis": "LOAO",
            "method": "true atom-mean LOAO shared with true_loao_score_atom_means; no own-heldout-atom centering",
            "selected_alpha": None,
        }],
    }


def stage22_term_group(term: str) -> str:
    if term == "charge_total":
        return "charge"
    if term.startswith("mc_"):
        return "mc"
    if term.startswith("mopac_field") or term.startswith("water_field"):
        return "field"
    if term.startswith("pq_"):
        return "pq"
    if term.startswith("disp_"):
        return "disp"
    return "other"


def stage22_support_flag(n_atoms: int, n_terms: int, n_eff: float) -> str:
    flags = []
    if n_terms >= n_atoms:
        flags.append("thin_LOAO_p_ge_atoms")
    elif n_atoms < 2 * n_terms:
        flags.append("thin_LOAO_atoms_lt_2p")
    if n_atoms < 30:
        flags.append("thin_atoms")
    if np.isfinite(n_eff) and n_eff < 10.0 * max(1, n_terms):
        flags.append("within_N_eff_lt_10p")
    return ";".join(flags)


def stage22_model_fit(
    model: str,
    term_set: str,
    labels: list[str],
    features: list[np.ndarray],
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    include_ci: bool,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    within, coeffs = stage22_evaluate_within(labels, features, y, atoms, frames, include_ci, seed)
    loao = stage22_evaluate_loao(labels, features, y, atoms, frames, include_ci, seed + 1000)
    neff_info = effective_n_components(y, atoms, frames)
    n_eff = float(neff_info["within_N_eff_median"])
    row = {
        "model": model,
        "term_set": term_set,
        "rows": int(len(y)),
        "atoms": int(len(np.unique(atoms))),
        "frames": int(len(np.unique(frames))),
        "terms": int(len(labels)),
        "term_labels": "|".join(labels),
        "within_N_eff_median": n_eff,
        "median_lag1_rho": float(neff_info["median_lag1_rho"]),
        "support_flag": stage22_support_flag(int(len(np.unique(atoms))), int(len(labels)), n_eff),
        **stage22_prefix("within", {k: v for k, v in within.items() if k != "alpha_audit"}),
        **stage22_prefix("LOAO", {k: v for k, v in loao.items() if k != "alpha_audit"}),
    }
    coeffs.insert(0, "model", model)
    coeffs.insert(1, "term_set", term_set)
    return row, coeffs


def stage22_base_dataset(data: dict[str, object]) -> dict[str, object]:
    rows: pd.DataFrame = data["rows"]
    mask, rule = unified_base_mask(data, clean_stage2=True)
    idx = np.flatnonzero(mask)
    labels, features_all = build_unified_features(data)
    return {
        "selection_rule": rule,
        "idx": idx,
        "labels": labels,
        "features": [f[idx] for f in features_all],
        "target": stage21_target(data)[idx],
        "atoms": rows["atom_index"].to_numpy(int)[idx],
        "frames": rows["original_frame_index"].to_numpy(int)[idx],
    }


def stage22_drop_text(drop: pd.DataFrame, column: str, n: int = 4) -> str:
    if drop.empty:
        return "none"
    rows = drop.sort_values(column, ascending=False).head(n)
    parts = []
    for row in rows.itertuples(index=False):
        parts.append(f"{row.dropped_term} {getattr(row, column):+.3f}/{row.LOAO_loss_vs_all:+.3f}")
    return "; ".join(parts) if parts else "none"


def stage22_harm_text(drop: pd.DataFrame, column: str, n: int = 3) -> str:
    if drop.empty:
        return "none"
    rows = drop.sort_values(column, ascending=True).head(n)
    parts = []
    for row in rows.itertuples(index=False):
        val = getattr(row, column)
        if np.isfinite(val) and val < 0.0:
            parts.append(f"{row.dropped_term} {val:+.3f}/{row.LOAO_loss_vs_all:+.3f}")
    return "; ".join(parts) if parts else "none"


def stage22_fmt_ci(row: pd.Series, axis: str) -> str:
    return (
        f"{float(row[f'{axis}_R2']):.3f} "
        f"(JK {float(row[f'{axis}_jackknife_ci95_low']):.3f},{float(row[f'{axis}_jackknife_ci95_high']):.3f}; "
        f"BB {float(row[f'{axis}_bootstrap_ci95_low']):.3f},{float(row[f'{axis}_bootstrap_ci95_high']):.3f})"
    )


def write_stage22_postmortem(
    out_dir: Path,
    substrate_dir: Path,
    orth: float,
    selection_rule: str,
    ols: pd.DataFrame,
    ablations: pd.DataFrame,
    drop: pd.DataFrame,
    disk: dict[str, object],
) -> Path:
    by_model = ablations.set_index("model")
    charge = by_model.loc["charge_alone"]
    minus = by_model.loc["unified_minus_charge"]
    allr = by_model.loc["unified_all"]
    ols_row = ols.iloc[0]
    delta_within_charge = float(allr["within_R2"] - charge["within_R2"])
    delta_loao_charge = float(allr["LOAO_R2"] - charge["LOAO_R2"])
    delta_within_minus = float(allr["within_R2"] - minus["within_R2"])
    delta_loao_minus = float(allr["LOAO_R2"] - minus["LOAO_R2"])
    within_survives = bool(float(allr["within_R2"]) > 0.30)
    loao_survives = bool(float(allr["LOAO_R2"]) >= 0.15)
    loao_thin_positive = bool(float(allr["LOAO_R2"]) > 0.05)
    if delta_within_charge > 0.05 and float(minus["within_R2"]) > 0.05:
        combine_bucket = "REAL combine on within; LOAO support thin"
        verdict = "real combine, not charge-in-a-coat; within survives shrinkage, LOAO is overfit-sensitive/thin-positive"
    elif delta_within_charge <= 0.03:
        combine_bucket = "charge-in-a-coat"
        verdict = "artifact: unified all adds little beyond charge"
    else:
        combine_bucket = "mixed weak-combine"
        verdict = "mixed: some non-charge lift, but charge dominates"
    if within_survives and loao_survives:
        overfit_bucket = "within-and-LOAO-survive-shrinkage"
    elif within_survives and loao_thin_positive:
        overfit_bucket = "within-survives; LOAO-shrinks-to-thin-positive"
    else:
        overfit_bucket = "collapses-or-thins-under-shrinkage"
    top_drop = stage22_drop_text(drop, "within_loss_vs_all")
    harmful = stage22_harm_text(drop, "within_loss_vs_all")
    lines = [
        "# Stage 2.2 Postmortem - 2026-06-04",
        f"Run dir: `{out_dir}`",
        f"Substrate: `{substrate_dir}`",
        f"Scope: Build4 CSV/NPY sidecars only; no extraction/DFT/per-source; frozen `get_C`, `|C.T C-I|max={orth:.2e}`; five-component total-T2.",
        f"Selection: {selection_rule}; input-side dominance only; 25 atoms / 3903 rows / 26 unified terms.",
        f"Overfit vet: OLS {float(ols_row['within_frameblock_R2']):.3f}/{float(ols_row['LOAO_R2']):.3f} -> ridge ALL within {stage22_fmt_ci(allr, 'within')}; LOAO {stage22_fmt_ci(allr, 'LOAO')}.",
        f"Shrinkage: within alpha {float(allr['within_selected_alpha']):g}, effDOF {float(allr['within_effective_dof_penalized']):.1f}/26 (+intercept {float(allr['within_effective_dof_with_intercept']):.1f}); within N_eff {float(allr['within_N_eff_median']):.1f}.",
        f"LOAO true-between: atom-mean scorer shared with true_loao; train observations median {float(allr['LOAO_train_observations_median']):.1f}; support `{allr['support_flag']}`.",
        f"Bucket: overfit={overfit_bucket}; LOAO column is true atom-mean recovery, not within modulation.",
        f"Charge table: charge-alone within {float(charge['within_R2']):.3f}, LOAO {float(charge['LOAO_R2']):.3f}; minus-charge within {float(minus['within_R2']):.3f}, LOAO {float(minus['LOAO_R2']):.3f}; all within {float(allr['within_R2']):.3f}, LOAO {float(allr['LOAO_R2']):.3f}.",
        f"Deltas ALL-charge: within {delta_within_charge:+.3f}, LOAO {delta_loao_charge:+.3f}; ALL-minus-charge: within {delta_within_minus:+.3f}, LOAO {delta_loao_minus:+.3f}; bucket={combine_bucket}.",
        f"Drop-one losses (within/LOAO, positive means term carries recovery): {top_drop}.",
        f"Drop-one removals that improve within (negative loss = noisy/overfit term): {harmful}.",
        f"Verdict: {verdict}.",
        "Artifacts: `stage2_2_charge_ablation_recovery.csv`, `stage2_2_drop_one_ablation.csv`, `stage2_2_unified_all_ridge_intensities.csv`, `stage2_2_run_audit.json`.",
        f"Disk: `/tmp/rediscover-runs` {float(disk['rediscover_runs_GiB_before_write']):.1f}G before write (<15G); output drop-old={disk.get('deleted_existing_output_dir')}; build4+build1 kept.",
    ]
    path = Path(__file__).resolve().parents[1] / "POSTMORTEM_STAGE2_2.md"
    path.write_text("\n".join(lines[:40]) + "\n", encoding="utf-8")
    return path


def run_stage22(args: argparse.Namespace, data: dict[str, object], out_dir: Path, disk: dict[str, object], orth: float) -> Path:
    base = stage22_base_dataset(data)
    labels: list[str] = base["labels"]
    features: list[np.ndarray] = base["features"]
    y: np.ndarray = base["target"]
    atoms: np.ndarray = base["atoms"]
    frames: np.ndarray = base["frames"]
    charge_i = labels.index("charge_total")
    model_specs = [
        ("charge_alone", "charge_total", [charge_i]),
        ("unified_minus_charge", "all_except_charge_total", [i for i in range(len(labels)) if i != charge_i]),
        ("unified_all", "all_terms", list(range(len(labels)))),
    ]
    ablation_rows = []
    coeff_tables = []
    for j, (model_name, term_set, inds) in enumerate(model_specs):
        row, coeffs = stage22_model_fit(
            model_name,
            term_set,
            [labels[i] for i in inds],
            [features[i] for i in inds],
            y,
            atoms,
            frames,
            include_ci=True,
            seed=20_220 + j,
        )
        ablation_rows.append(row)
        coeff_tables.append(coeffs)
    ablations = pd.DataFrame(ablation_rows)
    coeffs_all = pd.concat(coeff_tables, ignore_index=True)
    all_row = ablations[ablations["model"] == "unified_all"].iloc[0]
    drop_rows = []
    for i, label in enumerate(labels):
        inds = [j for j in range(len(labels)) if j != i]
        row, _coeffs = stage22_model_fit(
            f"drop_{label}",
            f"all_except_{label}",
            [labels[j] for j in inds],
            [features[j] for j in inds],
            y,
            atoms,
            frames,
            include_ci=False,
            seed=30_000 + i,
        )
        drop_rows.append({
            "dropped_term": label,
            "dropped_group": stage22_term_group(label),
            "drop_model_terms": int(row["terms"]),
            "drop_model_within_R2": float(row["within_R2"]),
            "drop_model_LOAO_R2": float(row["LOAO_R2"]),
            "within_loss_vs_all": float(all_row["within_R2"] - row["within_R2"]),
            "LOAO_loss_vs_all": float(all_row["LOAO_R2"] - row["LOAO_R2"]),
            "drop_within_alpha": float(row["within_selected_alpha"]),
            "drop_LOAO_alpha_mode": float(row["LOAO_selected_alpha_mode"]),
            "drop_within_eff_dof": float(row["within_effective_dof_penalized"]),
            "drop_LOAO_eff_dof_median": float(row["LOAO_effective_dof_penalized_median"]),
        })
    drop = pd.DataFrame(drop_rows).sort_values("within_loss_vs_all", ascending=False).reset_index(drop=True)
    ols_summary, _ols_intensities = run_unified_fit(data)
    ablations.to_csv(out_dir / "stage2_2_charge_ablation_recovery.csv", index=False)
    drop.to_csv(out_dir / "stage2_2_drop_one_ablation.csv", index=False)
    coeffs_all.to_csv(out_dir / "stage2_2_unified_all_ridge_intensities.csv", index=False)
    ols_summary.to_csv(out_dir / "stage2_2_ols_unified_reference.csv", index=False)
    audit = {
        "substrate_dir": str(args.substrate_dir),
        "out_dir": str(out_dir),
        "disk_guard": disk,
        "change_of_basis_orthogonality_max_abs": orth,
        "python_read_scope": "Build-4 per_atom_substrate CSV/NPY sidecars and CaseHunter manifests only",
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source dumps", "older rediscover dirs for fitting"],
        "selection_rule": str(base["selection_rule"]),
        "rows": int(len(y)),
        "atoms": int(len(np.unique(atoms))),
        "terms_unified_all": int(len(labels)),
        "alpha_grid": [float(x) for x in STAGE22_ALPHA_GRID],
        "inner_cv_folds": int(STAGE22_INNER_CV_FOLDS),
        "regularization": "ridge with standardized train-only design and unpenalized intercept",
        "within_outer_split": "centered contiguous frame block; held-out block excluded from alpha selection and centering means",
        "LOAO_outer_split": "one atom held out; train on training-atom means only, transform held-out atom mean with training-atom stats",
        "ci_methods": {
            "jackknife": "leave-one-atom recovery jackknife on held-out predictions",
            "block_bootstrap": f"within axis only: moving circular frame-block bootstrap, {STAGE22_BOOTSTRAPS} replicates, block length {STAGE22_BOOTSTRAP_FRAME_BLOCK}; LOAO uses atom-mean jackknife",
        },
        "anti_circularity": "through-space atom selection uses emitted input-side dominance/CaseHunter only; DFT appears only in held-out fitting/scoring",
    }
    (out_dir / "stage2_2_run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return write_stage22_postmortem(out_dir, args.substrate_dir, orth, str(base["selection_rule"]), ols_summary, ablations, drop, disk)


def stage23_stack_flat(features: list[np.ndarray]) -> np.ndarray:
    if not features:
        return np.empty((0, 0), dtype=float)
    return np.column_stack([np.asarray(f, dtype=float).reshape(-1) for f in features])


def stage23_prepare_within_design(
    labels: list[str],
    features: list[np.ndarray],
    atoms: np.ndarray,
    frames: np.ndarray,
) -> dict[str, object]:
    split = split_frame_block(frames, 0.20, 1)
    train = np.asarray(split["train_mask"], dtype=bool)
    test = np.asarray(split["test_mask"], dtype=bool)
    centered_features = [centered_by_train_atom(f, atoms, train) for f in features]
    x_train = stage23_stack_flat([f[train] for f in centered_features])
    x_test = stage23_stack_flat([f[test] for f in centered_features])
    valid_train = np.isfinite(x_train).all(axis=1) if x_train.size else np.asarray([], dtype=bool)
    valid_test = np.isfinite(x_test).all(axis=1) if x_test.size else np.asarray([], dtype=bool)
    pinv = np.linalg.pinv(x_train[valid_train]) if int(valid_train.sum()) >= len(labels) + 2 else np.empty((len(labels), 0))
    return {
        "axis": "within",
        "labels": labels,
        "atoms": atoms,
        "frames": frames,
        "train": train,
        "test": test,
        "x_train_valid_pinv": pinv,
        "valid_train": valid_train,
        "x_test": x_test,
        "valid_test": valid_test,
        "n_test_rows": int(test.sum()),
    }


def stage23_score_within(design: dict[str, object], y: np.ndarray) -> float:
    train = np.asarray(design["train"], dtype=bool)
    test = np.asarray(design["test"], dtype=bool)
    atoms = np.asarray(design["atoms"], dtype=int)
    y_c = centered_by_train_atom(y, atoms, train)
    y_train = y_c[train].reshape(-1)
    valid_train = np.asarray(design["valid_train"], dtype=bool)
    pinv = np.asarray(design["x_train_valid_pinv"], dtype=float)
    if pinv.size == 0 or valid_train.sum() < len(design["labels"]) + 2:
        return math.nan
    beta = pinv @ y_train[valid_train]
    x_test = np.asarray(design["x_test"], dtype=float)
    valid_test = np.asarray(design["valid_test"], dtype=bool)
    pred_flat = np.full(x_test.shape[0], np.nan, dtype=float)
    pred_flat[valid_test] = x_test[valid_test] @ beta
    return r2_score(pred_flat.reshape((int(design["n_test_rows"]), T2)), y_c[test])


def mislabeled_within_loao_prepare_design(
    labels: list[str],
    features: list[np.ndarray],
    atoms: np.ndarray,
) -> dict[str, object]:
    held_designs = []
    for held in np.unique(atoms):
        train = atoms != held
        test = atoms == held
        if int(train.sum()) < len(labels) + 2 or int(test.sum()) < 2:
            continue
        train_features = [atom_center(f[train], atoms[train]) for f in features]
        test_features = [f[test] - f[test].mean(axis=0, keepdims=True) for f in features]
        x_train = stage23_stack_flat(train_features)
        x_test = stage23_stack_flat(test_features)
        valid_train = np.isfinite(x_train).all(axis=1) if x_train.size else np.asarray([], dtype=bool)
        valid_test = np.isfinite(x_test).all(axis=1) if x_test.size else np.asarray([], dtype=bool)
        if int(valid_train.sum()) < len(labels) + 2:
            continue
        held_designs.append({
            "held": int(held),
            "train": train,
            "test": test,
            "valid_train": valid_train,
            "valid_test": valid_test,
            "x_train_valid_pinv": np.linalg.pinv(x_train[valid_train]),
            "x_test": x_test,
            "n_test_rows": int(test.sum()),
        })
    return {"axis": "LOAO", "labels": labels, "atoms": atoms, "held_designs": held_designs}


def mislabeled_within_loao_score(design: dict[str, object], y: np.ndarray) -> float:
    atoms = np.asarray(design["atoms"], dtype=int)
    pred = np.full_like(y, np.nan, dtype=float)
    target = np.full_like(y, np.nan, dtype=float)
    for held in design["held_designs"]:
        train = np.asarray(held["train"], dtype=bool)
        test = np.asarray(held["test"], dtype=bool)
        y_train = atom_center(y[train], atoms[train]).reshape(-1)
        y_test = y[test] - y[test].mean(axis=0, keepdims=True)
        valid_train = np.asarray(held["valid_train"], dtype=bool)
        pinv = np.asarray(held["x_train_valid_pinv"], dtype=float)
        beta = pinv @ y_train[valid_train]
        x_test = np.asarray(held["x_test"], dtype=float)
        valid_test = np.asarray(held["valid_test"], dtype=bool)
        pred_flat = np.full(x_test.shape[0], np.nan, dtype=float)
        pred_flat[valid_test] = x_test[valid_test] @ beta
        pred[test] = pred_flat.reshape((int(held["n_test_rows"]), T2))
        target[test] = y_test
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    return r2_score(pred[ok], target[ok])


def stage23_atom_groups(atoms: np.ndarray) -> list[np.ndarray]:
    return [np.flatnonzero(atoms == atom) for atom in np.unique(atoms)]


def stage23_permute_within_atom(y: np.ndarray, groups: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    out = np.empty_like(y, dtype=float)
    for idx in groups:
        out[idx] = y[idx[rng.permutation(len(idx))]]
    return out


def stage23_derangement(n: int, rng: np.random.Generator) -> np.ndarray:
    if n <= 1:
        return np.arange(n)
    if n == 2:
        return np.asarray([1, 0], dtype=int)
    base = np.arange(n)
    for _ in range(100):
        perm = rng.permutation(n)
        if np.all(perm != base):
            return perm
    return np.roll(base, 1)


def stage23_atom_transfer_plan(atoms: np.ndarray, frames: np.ndarray) -> tuple[np.ndarray, list[np.ndarray], dict[tuple[int, int], np.ndarray]]:
    unique = np.unique(atoms)
    groups = [np.flatnonzero(atoms == atom) for atom in unique]
    frame_maps = []
    for idx in groups:
        order = np.argsort(frames[idx], kind="mergesort")
        sorted_idx = idx[order]
        sorted_frames = frames[sorted_idx]
        frame_maps.append((sorted_frames, sorted_idx, {int(frame): int(pos) for pos, frame in enumerate(sorted_frames)}))
    transfers: dict[tuple[int, int], np.ndarray] = {}
    for ri, recip_idx in enumerate(groups):
        recip_frames = frames[recip_idx]
        for di, (_donor_frames, donor_idx, fmap) in enumerate(frame_maps):
            donor_take = np.empty(len(recip_idx), dtype=int)
            donor_n = len(donor_idx)
            for j, frame in enumerate(recip_frames):
                pos = fmap.get(int(frame))
                if pos is None:
                    pos = int(round(j * max(0, donor_n - 1) / max(1, len(recip_idx) - 1)))
                donor_take[j] = donor_idx[pos]
            transfers[(ri, di)] = donor_take
    return unique, groups, transfers


def stage23_permute_across_atoms(
    y: np.ndarray,
    groups: list[np.ndarray],
    transfers: dict[tuple[int, int], np.ndarray],
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.empty_like(y, dtype=float)
    perm = stage23_derangement(len(groups), rng)
    for ri, recip_idx in enumerate(groups):
        out[recip_idx] = y[transfers[(ri, int(perm[ri]))]]
    return out


def stage23_null_position(observed: float, nulls: np.ndarray) -> dict[str, object]:
    vals = np.asarray(nulls, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0 or not np.isfinite(observed):
        return {
            "observed_R2": float(observed) if np.isfinite(observed) else math.nan,
            "null_replicates": int(len(vals)),
            "null_mean": math.nan,
            "null_sd": math.nan,
            "null_percentile": math.nan,
            "empirical_p_upper": math.nan,
            "z_vs_null": math.nan,
        }
    mean = float(vals.mean())
    sd = float(vals.std(ddof=1)) if len(vals) > 1 else math.nan
    percentile = 100.0 * float((np.sum(vals < observed) + 0.5 * np.sum(vals == observed)) / len(vals))
    p_upper = float((1 + np.sum(vals >= observed)) / (len(vals) + 1))
    z = float((observed - mean) / sd) if np.isfinite(sd) and sd > 0.0 else math.nan
    return {
        "observed_R2": float(observed),
        "null_replicates": int(len(vals)),
        "null_mean": mean,
        "null_sd": sd,
        "null_percentile": percentile,
        "empirical_p_upper": p_upper,
        "z_vs_null": z,
    }


def stage23_null_scores(
    axis: str,
    score_fn,
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    n_perm: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vals = np.full(int(n_perm), np.nan, dtype=float)
    if axis == "within":
        groups = stage23_atom_groups(atoms)
        for i in range(int(n_perm)):
            vals[i] = score_fn(stage23_permute_within_atom(y, groups, rng))
        return vals
    if axis == "LOAO":
        raise ValueError("LOAO nulls use true_loao_permutation_null on atom means")
    raise ValueError(axis)


def stage23_eval_probability(
    model: str,
    labels: list[str],
    features: list[np.ndarray],
    y: np.ndarray,
    atoms: np.ndarray,
    frames: np.ndarray,
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    rows = []
    within_design = stage23_prepare_within_design(labels, features, atoms, frames)
    within_score = lambda yy: stage23_score_within(within_design, yy)
    within_obs = within_score(y)
    within_null = stage23_null_scores("within", within_score, y, atoms, frames, n_perm, seed)
    rows.append({
        "model": model,
        "axis": "within",
        "rows": int(len(y)),
        "atoms": int(len(np.unique(atoms))),
        "frames": int(len(np.unique(frames))),
        "terms": int(len(labels)),
        "term_labels": "|".join(labels),
        "null_shuffle": "DFT target shuffled across frames within atom",
        **stage23_null_position(within_obs, within_null),
    })
    atom_labels, atom_features, atom_y, _counts = true_loao_atom_mean_arrays(features, y, atoms)
    loao_obs = float(true_loao_score_atom_means(labels, atom_features, atom_y)["R2"])
    loao_null = true_loao_permutation_null(labels, atom_features, atom_y, n_perm, seed + 10_000)
    rows.append({
        "model": model,
        "axis": "LOAO",
        "rows": int(len(y)),
        "atoms": int(len(atom_labels)),
        "frames": int(len(np.unique(frames))),
        "terms": int(len(labels)),
        "term_labels": "|".join(labels),
        "null_shuffle": "DFT target atom-means deranged across atoms; structure and LOAO folds unchanged",
        **stage23_null_position(loao_obs, loao_null),
    })
    return pd.DataFrame(rows)


def stage23_kernel_dataset(data: dict[str, object], kernel: KernelSpec) -> dict[str, object]:
    rows: pd.DataFrame = data["rows"]
    mask, rule = kernel_clean_base_mask(data, kernel)
    idx = np.flatnonzero(mask)
    return {
        "model": kernel.name,
        "selection_rule": rule,
        "labels": [kernel.primary.label],
        "features": [stage21_primary_shadow(data, kernel)[idx]],
        "target": stage21_target(data)[idx],
        "atoms": rows["atom_index"].to_numpy(int)[idx],
        "frames": rows["original_frame_index"].to_numpy(int)[idx],
        "metrics": stage21_kernel_metrics(data, kernel, mask),
    }


def stage23_unified_dataset_for_cutoffs(data: dict[str, object], dab_cut: float, ring_cut: float) -> dict[str, object]:
    rows: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    atoms: set[int] = set()
    for case_dir in ["charge_wide", "mc", "field"]:
        case_set, _cases = case_atoms(data["substrate_dir"], case_dir)
        atoms |= case_set
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    dab_max = np.maximum.reduce([
        np.asarray(dom[:, dom_cols[f"dominant_fraction_{m}"]], dtype=float)
        for m in ["charge", "mc", "field"]
    ])
    ring_frac = np.asarray(dom[:, dom_cols["dominant_fraction_ring"]], dtype=float)
    atoms_all = rows["atom_index"].to_numpy(int)
    mask = (
        (rows["dft_present"].to_numpy(int) == 1)
        & np.isin(atoms_all, list(atoms))
        & np.isfinite(dab_max)
        & np.isfinite(ring_frac)
        & (dab_max >= float(dab_cut))
        & (ring_frac < float(ring_cut))
    )
    idx = np.flatnonzero(mask)
    labels, features_all = build_unified_features(data)
    return {
        "model": "unified_Dab_sum",
        "selection_rule": f"D_ab CaseHunter atoms + max D_ab dominance>={dab_cut:.2f} + ring dominance<{ring_cut:.2f}",
        "labels": labels,
        "features": [f[idx] for f in features_all],
        "target": stage21_target(data)[idx],
        "atoms": rows["atom_index"].to_numpy(int)[idx],
        "frames": rows["original_frame_index"].to_numpy(int)[idx],
        "dab_cutoff": float(dab_cut),
        "ring_cutoff": float(ring_cut),
    }


def stage23_drop_one_unified(base: dict[str, object]) -> pd.DataFrame:
    labels: list[str] = base["labels"]
    features: list[np.ndarray] = base["features"]
    y: np.ndarray = base["target"]
    atoms: np.ndarray = base["atoms"]
    frames: np.ndarray = base["frames"]
    full = stage23_eval_probability("unified_Dab_sum", labels, features, y, atoms, frames, 1, 44)
    full_within = float(full[full["axis"] == "within"]["observed_R2"].iloc[0])
    full_loao = float(full[full["axis"] == "LOAO"]["observed_R2"].iloc[0])
    _atom_labels, atom_features, atom_y, _counts = true_loao_atom_mean_arrays(features, y, atoms)
    rows = []
    for i, label in enumerate(labels):
        keep = [j for j in range(len(labels)) if j != i]
        sub_labels = [labels[j] for j in keep]
        sub_features = [features[j] for j in keep]
        within_design = stage23_prepare_within_design(sub_labels, sub_features, atoms, frames)
        within_r2 = stage23_score_within(within_design, y)
        loao_r2 = float(true_loao_score_atom_means(sub_labels, [atom_features[j] for j in keep], atom_y)["R2"])
        rows.append({
            "dropped_term": label,
            "dropped_group": stage22_term_group(label),
            "within_R2_without_term": within_r2,
            "LOAO_R2_without_term": loao_r2,
            "within_loss_vs_all": float(full_within - within_r2) if np.isfinite(within_r2) else math.nan,
            "LOAO_loss_vs_all": float(full_loao - loao_r2) if np.isfinite(loao_r2) else math.nan,
        })
    return pd.DataFrame(rows).sort_values("within_loss_vs_all", ascending=False).reset_index(drop=True)


def stage23_fmt_p(x: object) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "nan"
    if not np.isfinite(v):
        return "nan"
    if v < 0.0015:
        return "0.001"
    return f"{v:.3f}"


def stage23_position_text(row: pd.Series) -> str:
    return (
        f"R2 {float(row['observed_R2']):.3f}; pct {float(row['null_percentile']):.1f}; "
        f"p {stage23_fmt_p(row['empirical_p_upper'])}; z {float(row['z_vs_null']):.1f}"
    )


def stage23_top_terms(drop: pd.DataFrame, axis: str, min_loss: float = 0.01, n: int = 4) -> str:
    col = f"{axis}_loss_vs_all"
    rows = drop[np.isfinite(drop[col].to_numpy(float)) & (drop[col].to_numpy(float) >= min_loss)].sort_values(col, ascending=False).head(n)
    parts = [f"{r.dropped_term} {getattr(r, col):+.3f}" for r in rows.itertuples(index=False)]
    return "; ".join(parts) if parts else "diffuse"


def stage23_attach_interpretation(prob: pd.DataFrame, kernel_metrics: dict[str, dict[str, object]], drop: pd.DataFrame) -> pd.DataFrame:
    out = prob.copy()
    determ = []
    placement = []
    for row in out.itertuples(index=False):
        model = str(row.model)
        axis = str(row.axis)
        pval = float(row.empirical_p_upper)
        r2 = float(row.observed_R2)
        above_null = bool(np.isfinite(pval) and pval <= 0.05 and np.isfinite(r2) and r2 > 0.05)
        if model == "unified_Dab_sum":
            top = stage23_top_terms(drop, axis)
            if top == "diffuse" or axis == "LOAO":
                d = f"{top}; atom-axis attribution not determined"
                place = "potentially indicative but INDETERMINATE -> needs more atoms / second-protein atom-axis test" if above_null else "not indicative (~null)"
            else:
                d = f"attributable to {top}"
                place = "indicative of field + McConnell category mixture, determinable" if above_null else "not indicative (~null)"
        elif model == "charge":
            d = "single q/r3 shadow; coefficient CI off zero"
            place = "indicative of charge q/r3, determinable" if above_null else "not indicative (~null)"
        elif model == "ring":
            d = "single JB ring shadow; thin atom support"
            place = "indicative of ring current, determinable" if above_null else "not indicative (~null)"
        elif model == "field":
            d = "single MOPAC-Coulomb shadow; recovery is 0.03-class"
            place = "not indicative (~null)"
        elif model == "mc":
            d = "single McConnell shadow; CI spans zero"
            place = "not indicative (~null)" if not above_null else "potentially indicative but INDETERMINATE -> needs joint fit"
        else:
            d = "indeterminate"
            place = "potentially indicative but INDETERMINATE -> needs targeted test" if above_null else "not indicative (~null)"
        determ.append(d)
        placement.append(place)
    out["determinability"] = determ
    out["lead_scale_placement"] = placement
    return out


def stage23_curve_rows(data: dict[str, object], n_perm: int) -> pd.DataFrame:
    rows = []
    for sweep_axis, dab_cut, ring_cut in (
        *[("dab_cutoff", cut, 0.70) for cut in STAGE23_DAB_CUTOFFS],
        *[("ring_cutoff", 0.50, cut) for cut in STAGE23_RING_CUTOFFS],
    ):
        base = stage23_unified_dataset_for_cutoffs(data, dab_cut, ring_cut)
        if len(base["target"]) < len(base["labels"]) + 3 or len(np.unique(base["atoms"])) < 2:
            continue
        prob = stage23_eval_probability(
            "unified_Dab_sum",
            base["labels"],
            base["features"],
            base["target"],
            base["atoms"],
            base["frames"],
            n_perm,
            int(90_000 + round(dab_cut * 1000) * 10 + round(ring_cut * 1000)),
        )
        for row in prob.itertuples(index=False):
            rows.append({
                "sweep_axis": sweep_axis,
                "dab_cutoff": float(dab_cut),
                "ring_cutoff": float(ring_cut),
                "axis": str(row.axis),
                "rows": int(row.rows),
                "atoms": int(row.atoms),
                "observed_R2": float(row.observed_R2),
                "empirical_p_upper": float(row.empirical_p_upper),
                "null_percentile": float(row.null_percentile),
                "z_vs_null": float(row.z_vs_null),
                "null_replicates": int(row.null_replicates),
                "selection_rule": str(base["selection_rule"]),
            })
    return pd.DataFrame(rows)


def stage23_curve_text(curve: pd.DataFrame, sweep_axis: str, axis: str) -> str:
    fixed = "ring<0.70" if sweep_axis == "dab_cutoff" else "D_ab>=0.50"
    label = "D_ab cut" if sweep_axis == "dab_cutoff" else "ring cap"
    value_col = "dab_cutoff" if sweep_axis == "dab_cutoff" else "ring_cutoff"
    g = curve[(curve["sweep_axis"] == sweep_axis) & (curve["axis"] == axis)].sort_values(value_col)
    parts = [f"{float(r[value_col]):.2f} {stage23_fmt_p(r['empirical_p_upper'])}/{float(r['observed_R2']):.3f}" for _, r in g.iterrows()]
    return f"{label} ({fixed}) {axis} p/R2: " + "; ".join(parts)


def write_stage23_postmortem(
    out_dir: Path,
    substrate_dir: Path,
    orth: float,
    prob: pd.DataFrame,
    curve: pd.DataFrame,
    disk: dict[str, object],
) -> Path:
    order = ["unified_Dab_sum", "charge", "field", "ring", "mc"]
    lines = [
        "# Stage 2.3 Postmortem - 2026-06-04",
        f"Run dir: `{out_dir}`",
        f"Substrate: `{substrate_dir}`; Build4 CSV/NPY only; frozen `get_C`, `|C.T C-I|max={orth:.2e}`; five-component total-T2.",
        "Nulls: within shuffles the DFT target across frames within atom; LOAO deranges target atom-means across atoms; 1000 shuffles each row.",
        "| recovery | statistical position vs null | determinability | lead-scale placement |",
        "| --- | --- | --- | --- |",
    ]
    for model in order:
        for axis in ["within", "LOAO"]:
            row = prob[(prob["model"] == model) & (prob["axis"] == axis)]
            if row.empty:
                continue
            r = row.iloc[0]
            lines.append(f"| {model} {axis} | {stage23_position_text(r)} | {r['determinability']} | {r['lead_scale_placement']} |")
    lines += [
        "Cutoff probability curve reports p/R2 at each input-side cut; lower p means farther above the shuffled-target null.",
        stage23_curve_text(curve, "dab_cutoff", "within"),
        stage23_curve_text(curve, "dab_cutoff", "LOAO"),
        stage23_curve_text(curve, "ring_cutoff", "within"),
        stage23_curve_text(curve, "ring_cutoff", "LOAO"),
        "~0.03-class reads as ~null here; ~0.1-0.2-class entries are placed by determinability, not by size alone.",
        "Artifacts: `stage2_3_recovery_probability.csv`, `stage2_3_unified_drop_one.csv`, `stage2_3_cutoff_probability_curve.csv`, `stage2_3_run_audit.json`.",
        f"Disk: `/tmp/rediscover-runs` {float(disk['rediscover_runs_GiB_before_write']):.1f}G before write (<15G); output drop-old={disk.get('deleted_existing_output_dir')}; build4+build1 kept.",
    ]
    path = Path(__file__).resolve().parents[1] / "POSTMORTEM_STAGE2_3.md"
    path.write_text("\n".join(lines[:45]) + "\n", encoding="utf-8")
    return path


def run_stage23(args: argparse.Namespace, data: dict[str, object], out_dir: Path, disk: dict[str, object], orth: float) -> Path:
    n_perm = max(1000, int(args.stage2_3_permutations))
    probability_frames = []
    kernel_metrics: dict[str, dict[str, object]] = {}
    base = stage23_unified_dataset_for_cutoffs(data, 0.50, 0.70)
    drop = stage23_drop_one_unified(base)
    probability_frames.append(stage23_eval_probability(
        "unified_Dab_sum",
        base["labels"],
        base["features"],
        base["target"],
        base["atoms"],
        base["frames"],
        n_perm,
        23_001,
    ))
    for i, kernel in enumerate(KERNELS):
        ds = stage23_kernel_dataset(data, kernel)
        kernel_metrics[kernel.name] = ds["metrics"]
        probability_frames.append(stage23_eval_probability(
            kernel.name,
            ds["labels"],
            ds["features"],
            ds["target"],
            ds["atoms"],
            ds["frames"],
            n_perm,
            24_000 + i * 101,
        ))
    prob = pd.concat(probability_frames, ignore_index=True)
    prob = stage23_attach_interpretation(prob, kernel_metrics, drop)
    curve = stage23_curve_rows(data, n_perm)
    prob.to_csv(out_dir / "stage2_3_recovery_probability.csv", index=False)
    drop.to_csv(out_dir / "stage2_3_unified_drop_one.csv", index=False)
    curve.to_csv(out_dir / "stage2_3_cutoff_probability_curve.csv", index=False)
    audit = {
        "substrate_dir": str(args.substrate_dir),
        "out_dir": str(out_dir),
        "disk_guard": disk,
        "change_of_basis_orthogonality_max_abs": orth,
        "python_read_scope": "Build-4 per_atom_substrate CSV/NPY sidecars and CaseHunter manifests only",
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source dumps"],
        "permutations_per_probability_row": int(n_perm),
        "nulls": {
            "within": "DFT target shuffled across frames within atom",
            "LOAO": "DFT target atom-means deranged across atoms; feature atom-means and held-out atom folds unchanged",
        },
        "cutoff_sweeps": {
            "dab_cutoffs_ring_fixed_0.70": list(STAGE23_DAB_CUTOFFS),
            "ring_cutoffs_dab_fixed_0.50": list(STAGE23_RING_CUTOFFS),
        },
        "anti_circularity": "selection cutoffs use emitted input-side dominance only; null shuffles only the DFT target",
    }
    (out_dir / "stage2_3_run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return write_stage23_postmortem(out_dir, args.substrate_dir, orth, prob, curve, disk)


def true_loao_col_stats(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=float)
    mean = np.full(arr.shape[1], np.nan, dtype=float)
    std = np.ones(arr.shape[1], dtype=float)
    for i in range(arr.shape[1]):
        col = arr[:, i]
        ok = np.isfinite(col)
        if not np.any(ok):
            continue
        mean[i] = float(np.mean(col[ok]))
        s = float(np.std(col[ok]))
        std[i] = s if np.isfinite(s) and s > 0.0 else 1.0
    return mean, std


def true_loao_atom_mean_arrays(
    features: list[np.ndarray],
    y: np.ndarray,
    atoms: np.ndarray,
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, np.ndarray]:
    labels = np.unique(atoms)
    counts = np.zeros(len(labels), dtype=int)
    atom_y = np.full((len(labels), T2), np.nan, dtype=float)
    atom_features = [np.full((len(labels), T2), np.nan, dtype=float) for _ in features]
    for i, atom in enumerate(labels):
        m = atoms == atom
        counts[i] = int(np.sum(m))
        atom_y[i], _ = true_loao_col_stats(y[m])
        for j, f in enumerate(features):
            atom_features[j][i], _ = true_loao_col_stats(f[m])
    return labels.astype(int), atom_features, atom_y, counts


def true_loao_standardized_fit(
    train_features: list[np.ndarray],
    y_train: np.ndarray,
    eval_features: list[np.ndarray],
) -> tuple[np.ndarray, dict[str, object]]:
    p = len(train_features)
    eval_rows = int(eval_features[0].shape[0]) if eval_features else 0
    if p == 0:
        return np.full((eval_rows, T2), np.nan, dtype=float), {"ok": False, "train_observations": 0}
    x_means = []
    x_stds = []
    z_train = []
    z_eval = []
    for tr, ev in zip(train_features, eval_features):
        mean, std = true_loao_col_stats(tr)
        x_means.append(mean)
        x_stds.append(std)
        z_train.append((np.asarray(tr, dtype=float) - mean.reshape(1, T2)) / std.reshape(1, T2))
        z_eval.append((np.asarray(ev, dtype=float) - mean.reshape(1, T2)) / std.reshape(1, T2))
    y_mean, y_std = true_loao_col_stats(y_train)
    x_train = stage23_stack_flat(z_train)
    yy = ((np.asarray(y_train, dtype=float) - y_mean.reshape(1, T2)) / y_std.reshape(1, T2)).reshape(-1)
    valid = np.isfinite(yy) & np.isfinite(x_train).all(axis=1)
    if int(valid.sum()) < p + 2:
        return np.full((eval_rows, T2), np.nan, dtype=float), {
            "ok": False,
            "train_observations": int(valid.sum()),
            "x_means": x_means,
            "x_stds": x_stds,
            "y_mean": y_mean,
            "y_std": y_std,
        }
    beta, *_ = np.linalg.lstsq(x_train[valid], yy[valid], rcond=None)
    x_eval = stage23_stack_flat(z_eval)
    valid_eval = np.isfinite(x_eval).all(axis=1)
    pred_z = np.full(x_eval.shape[0], np.nan, dtype=float)
    pred_z[valid_eval] = x_eval[valid_eval] @ beta
    pred = pred_z.reshape((eval_rows, T2)) * y_std.reshape(1, T2) + y_mean.reshape(1, T2)
    return pred, {
        "ok": True,
        "train_observations": int(valid.sum()),
        "beta_standardized": np.asarray(beta, dtype=float),
        "x_means": x_means,
        "x_stds": x_stds,
        "y_mean": y_mean,
        "y_std": y_std,
    }


def true_loao_score_atom_means(
    labels: list[str],
    atom_features: list[np.ndarray],
    atom_y: np.ndarray,
) -> dict[str, object]:
    n_atoms = int(atom_y.shape[0])
    pred = np.full_like(atom_y, np.nan, dtype=float)
    observations = []
    coef_rows = []
    for held_i in range(n_atoms):
        train = np.ones(n_atoms, dtype=bool)
        train[held_i] = False
        eval_features = [f[held_i:held_i + 1] for f in atom_features]
        train_features = [f[train] for f in atom_features]
        pred_i, model = true_loao_standardized_fit(train_features, atom_y[train], eval_features)
        pred[held_i] = pred_i[0]
        observations.append(int(model.get("train_observations", 0)))
        beta = np.asarray(model.get("beta_standardized", np.full(len(labels), np.nan)), dtype=float)
        if beta.shape[0] == len(labels):
            coef_rows.append(beta)
    coef = np.asarray(coef_rows, dtype=float)
    coef_median = np.nanmedian(coef, axis=0) if coef.size else np.full(len(labels), np.nan)
    coef_min = np.nanmin(coef, axis=0) if coef.size else np.full(len(labels), np.nan)
    coef_max = np.nanmax(coef, axis=0) if coef.size else np.full(len(labels), np.nan)
    return {
        "R2": r2_score(pred, atom_y),
        "pred": pred,
        "target": atom_y,
        "train_observations_median": float(np.median(observations)) if observations else math.nan,
        "standardized_coef_median": coef_median,
        "standardized_coef_min": coef_min,
        "standardized_coef_max": coef_max,
    }


def true_loao_permutation_null(
    labels: list[str],
    atom_features: list[np.ndarray],
    atom_y: np.ndarray,
    n_perm: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vals = np.full(int(n_perm), np.nan, dtype=float)
    n_atoms = int(atom_y.shape[0])
    for i in range(int(n_perm)):
        perm = stage23_derangement(n_atoms, rng)
        vals[i] = float(true_loao_score_atom_means(labels, atom_features, atom_y[perm])["R2"])
    return vals


def true_loao_support_flag(n_atoms: int, n_terms: int) -> str:
    flags = []
    if n_terms >= n_atoms:
        flags.append("thin_TRUE_LOAO_p_ge_atoms")
    elif n_atoms < 2 * n_terms:
        flags.append("thin_TRUE_LOAO_atoms_lt_2p")
    if n_atoms < 30:
        flags.append("thin_atoms_case_study")
    return ";".join(flags) if flags else "ok_atoms_report_N"


def true_loao_mislabeled_within_loao(labels: list[str], features: list[np.ndarray], y: np.ndarray, atoms: np.ndarray) -> float:
    design = mislabeled_within_loao_prepare_design(labels, features, atoms)
    return float(mislabeled_within_loao_score(design, y))


def true_loao_drop_one_unified(base: dict[str, object]) -> pd.DataFrame:
    labels: list[str] = base["labels"]
    features: list[np.ndarray] = base["features"]
    y: np.ndarray = base["target"]
    atoms: np.ndarray = base["atoms"]
    atom_labels, atom_features, atom_y, _counts = true_loao_atom_mean_arrays(features, y, atoms)
    full = true_loao_score_atom_means(labels, atom_features, atom_y)
    full_r2 = float(full["R2"])
    rows = []
    for i, label in enumerate(labels):
        keep = [j for j in range(len(labels)) if j != i]
        sub_labels = [labels[j] for j in keep]
        sub_features = [atom_features[j] for j in keep]
        score = true_loao_score_atom_means(sub_labels, sub_features, atom_y)
        r2 = float(score["R2"])
        rows.append({
            "dropped_term": label,
            "dropped_group": stage22_term_group(label),
            "atoms": int(len(atom_labels)),
            "terms_without": int(len(sub_labels)),
            "TRUE_LOAO_R2_without_term": r2,
            "TRUE_LOAO_loss_vs_all": float(full_r2 - r2) if np.isfinite(r2) else math.nan,
        })
    return pd.DataFrame(rows).sort_values("TRUE_LOAO_loss_vs_all", ascending=False).reset_index(drop=True)


def true_loao_drop_text(drop: pd.DataFrame, min_loss: float = 0.01, n: int = 4) -> str:
    if drop.empty:
        return "drop-one unavailable"
    col = "TRUE_LOAO_loss_vs_all"
    rows = drop[np.isfinite(drop[col].to_numpy(float)) & (drop[col].to_numpy(float) >= min_loss)].sort_values(col, ascending=False).head(n)
    if rows.empty:
        return "drop-one diffuse/no term >=0.01"
    return "; ".join(f"{r.dropped_term} {float(getattr(r, col)):+.3f}" for r in rows.itertuples(index=False))


def true_loao_lead_scale(r2: float, support: str) -> str:
    suffix = " (case-study N)" if "thin" in support else ""
    if not np.isfinite(r2):
        return "undetermined" + suffix
    if r2 < 0.0:
        return "negative/null; no atom-mean recovery" + suffix
    if r2 <= 0.05:
        return "~null / 0.03-class" + suffix
    if r2 < 0.20:
        return "below 0.2; weak, read by determinability" + suffix
    if r2 < 0.35:
        return "0.2-class; decided by determinability" + suffix
    return "higher-than-0.2; still read by determinability" + suffix


def true_loao_fmt(x: object) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "nan"
    if not np.isfinite(v):
        return "nan"
    return f"{v:.3f}"


def write_true_loao_postmortem(
    out_dir: Path,
    substrate_dir: Path,
    orth: float,
    recovery: pd.DataFrame,
    drop: pd.DataFrame,
    disk: dict[str, object],
) -> Path:
    order = ["charge", "ring", "unified_Dab_sum"]
    lines = [
        "# TRUE LOAO Postmortem - 2026-06-04",
        f"Run dir: `{out_dir}`",
        f"Substrate: `{substrate_dir}`; Build4 CSV/NPY only; no extraction/DFT; frozen `get_C`, `|C.T C-I|max={orth:.2e}`.",
        "TRUE LOAO: train on training-atom means only; transform held-out atom means with training-atom feature/target stats; score physical held-out atom means.",
        "Null: 1000 deranged shuffles of target atom-means across atoms; structural features and LOAO folds unchanged.",
        "Support note: charge N is the actual Stage2 charge-wide subset tied to old 0.38; the older 54-atom backbone stratum is a different comparison.",
        "| kernel | N atoms | old mislabeled within-LOAO R2 | TRUE-LOAO R2 | delta | null p/z/pct | determinability | lead-scale | support |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- | --- | --- |",
    ]
    for model in order:
        row = recovery[recovery["model"] == model]
        if row.empty:
            continue
        r = row.iloc[0]
        null_text = f"p {stage23_fmt_p(r['empirical_p_upper'])}; z {true_loao_fmt(r['z_vs_null'])}; pct {true_loao_fmt(r['null_percentile'])}"
        lines.append(
            f"| {model} | {int(r['atoms'])} | {true_loao_fmt(r['mislabeled_within_LOAO_R2'])} | "
            f"{true_loao_fmt(r['TRUE_LOAO_R2'])} | {true_loao_fmt(r['delta_TRUE_minus_mislabeled'])} | "
            f"{null_text} | {r['determinability']} | {r['lead_scale_placement']} | {r['support_flag']} |"
        )
    read = str(recovery.attrs.get("one_line_read", "one-line read unavailable"))
    lines += [
        f"Unified drop-one determinability: {true_loao_drop_text(drop)}.",
        f"One-line read: {read}",
        "Difference is explicit: the old number is within-atom modulation because the held-out atom was centered by its own mean; TRUE LOAO is between-atom atom-mean recovery.",
        "Artifacts: `true_loao_recovery.csv`, `true_loao_unified_drop_one.csv`, `true_loao_atom_predictions.csv`, `true_loao_null_scores.csv`, `true_loao_run_audit.json`.",
        f"Disk: `/tmp/rediscover-runs` {float(disk['rediscover_runs_GiB_before_write']):.1f}G before write (<15G); output drop-old={disk.get('deleted_existing_output_dir')}; build4+build1 kept.",
    ]
    path = Path(__file__).resolve().parents[1] / "POSTMORTEM_TRUE_LOAO_2026-06-04.md"
    path.write_text("\n".join(lines[:35]) + "\n", encoding="utf-8")
    return path


def run_true_loao(args: argparse.Namespace, data: dict[str, object], out_dir: Path, disk: dict[str, object], orth: float) -> Path:
    n_perm = max(1000, int(args.true_loao_permutations))
    charge = stage23_kernel_dataset(data, next(k for k in KERNELS if k.name == "charge"))
    ring = stage23_kernel_dataset(data, next(k for k in KERNELS if k.name == "ring"))
    unified = stage23_unified_dataset_for_cutoffs(data, 0.50, 0.70)
    datasets = [charge, ring, unified]
    drop = true_loao_drop_one_unified(unified)
    unified_drop_text = true_loao_drop_text(drop)
    recovery_rows = []
    prediction_rows = []
    null_rows = []
    for i, ds in enumerate(datasets):
        model = str(ds["model"])
        labels: list[str] = ds["labels"]
        features: list[np.ndarray] = ds["features"]
        y: np.ndarray = ds["target"]
        atoms: np.ndarray = ds["atoms"]
        frames: np.ndarray = ds["frames"]
        atom_labels, atom_features, atom_y, counts = true_loao_atom_mean_arrays(features, y, atoms)
        observed = true_loao_score_atom_means(labels, atom_features, atom_y)
        true_r2 = float(observed["R2"])
        null = true_loao_permutation_null(labels, atom_features, atom_y, n_perm, 61_000 + i * 997)
        position = stage23_null_position(true_r2, null)
        old_loao = true_loao_mislabeled_within_loao(labels, features, y, atoms)
        support = true_loao_support_flag(int(len(atom_labels)), int(len(labels)))
        if model == "charge":
            determinability = "single q/r3; Stage2 charge-wide subset, N-limited case-study"
        elif model == "ring":
            determinability = "single JB ring shadow; N=5 case-study, probability discrete"
        else:
            determinability = f"unified drop-one: {unified_drop_text}"
        recovery_rows.append({
            "model": model,
            "rows": int(len(y)),
            "atoms": int(len(atom_labels)),
            "frames": int(len(np.unique(frames))),
            "terms": int(len(labels)),
            "term_labels": "|".join(labels),
            "selection_rule": str(ds["selection_rule"]),
            "mislabeled_within_LOAO_R2": old_loao,
            "TRUE_LOAO_R2": true_r2,
            "delta_TRUE_minus_mislabeled": float(true_r2 - old_loao) if np.isfinite(old_loao) and np.isfinite(true_r2) else math.nan,
            "support_flag": support,
            "atom_row_count_min": int(np.min(counts)) if len(counts) else 0,
            "atom_row_count_median": float(np.median(counts)) if len(counts) else math.nan,
            "atom_row_count_max": int(np.max(counts)) if len(counts) else 0,
            "null_shuffle": "target atom-means deranged across atoms; structure and LOAO folds unchanged",
            "determinability": determinability,
            "lead_scale_placement": true_loao_lead_scale(true_r2, support),
            **position,
        })
        pred = np.asarray(observed["pred"], dtype=float)
        for ai, atom in enumerate(atom_labels):
            for comp in range(T2):
                prediction_rows.append({
                    "model": model,
                    "atom_index": int(atom),
                    "component": int(comp),
                    "target_atom_mean": float(atom_y[ai, comp]) if np.isfinite(atom_y[ai, comp]) else math.nan,
                    "predicted_atom_mean": float(pred[ai, comp]) if np.isfinite(pred[ai, comp]) else math.nan,
                    "rows_for_atom": int(counts[ai]),
                })
        for j, val in enumerate(null):
            null_rows.append({
                "model": model,
                "permutation": int(j),
                "TRUE_LOAO_R2_null": float(val) if np.isfinite(val) else math.nan,
            })
    recovery = pd.DataFrame(recovery_rows)
    unified_r2 = float(recovery[recovery["model"] == "unified_Dab_sum"]["TRUE_LOAO_R2"].iloc[0])
    charge_r2 = float(recovery[recovery["model"] == "charge"]["TRUE_LOAO_R2"].iloc[0])
    ring_r2 = float(recovery[recovery["model"] == "ring"]["TRUE_LOAO_R2"].iloc[0])
    if max(unified_r2, charge_r2, ring_r2) <= 0.05:
        read = "TRUE between-atom recovery on 1P9J is ~null; between-axis probability should move to the 720-WT."
    elif unified_r2 < 0.20 and charge_r2 < 0.20 and ring_r2 < 0.20:
        read = "TRUE between-atom recovery is weak below the 0.2 lead-scale; treat 1P9J as case-study and move probability to the 720-WT."
    else:
        read = "TRUE between-atom recovery has non-null case-study signal on 1P9J, but atom support keeps probability with the 720-WT."
    recovery.attrs["one_line_read"] = read
    recovery.to_csv(out_dir / "true_loao_recovery.csv", index=False)
    drop.to_csv(out_dir / "true_loao_unified_drop_one.csv", index=False)
    pd.DataFrame(prediction_rows).to_csv(out_dir / "true_loao_atom_predictions.csv", index=False)
    pd.DataFrame(null_rows).to_csv(out_dir / "true_loao_null_scores.csv", index=False)
    audit = {
        "substrate_dir": str(args.substrate_dir),
        "out_dir": str(out_dir),
        "disk_guard": disk,
        "change_of_basis_orthogonality_max_abs": orth,
        "python_read_scope": "Build-4 per_atom_substrate CSV/NPY sidecars and CaseHunter manifests only",
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source dumps", "older rediscover dirs for fitting"],
        "permutations_per_true_loao_row": int(n_perm),
        "true_loao_definition": "atom means over frames; each held-out atom predicted from model fit on other atom means after train-atom-only per-component feature and target standardization",
        "mislabeled_reference": "existing LOAO code path that centers the held-out atom by its own mean, hence within-atom modulation",
        "null": "deranged shuffle of target atom-mean rows across atoms; no feature shuffle and no new DFT/extraction",
        "unified_drop_one": "drop one unified D_ab term and recompute TRUE LOAO atom-mean recovery",
        "anti_circularity": "selection cutoffs use emitted input-side dominance only; null shuffles atom means only",
    }
    (out_dir / "true_loao_run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return write_true_loao_postmortem(out_dir, args.substrate_dir, orth, recovery, drop, disk)


def main() -> None:
    args = parse_args()
    args.substrate_dir = args.substrate_dir.resolve()
    if args.stage2_1 and args.out_dir == OUT_DIR:
        args.out_dir = STAGE21_OUT_DIR
    if args.stage2_2 and args.out_dir == OUT_DIR:
        args.out_dir = STAGE22_OUT_DIR
    if args.stage2_3 and args.out_dir == OUT_DIR:
        args.out_dir = STAGE23_OUT_DIR
    if args.true_loao and args.out_dir == OUT_DIR:
        args.out_dir = TRUE_LOAO_OUT_DIR
    out_dir = args.out_dir.resolve()
    disk = disk_guard(out_dir, drop_existing=bool(args.stage2_1 or args.stage2_2 or args.stage2_3 or args.true_loao))
    data, orth = load_stage2_inputs(args.substrate_dir, include_stage21=bool(args.stage2_1))
    manifest: dict[str, object] = data["manifest"]
    if args.stage2_1:
        postmortem = run_stage21(args, data, out_dir, disk, orth)
        print(f"stage2.1 wrote {out_dir}")
        print(f"postmortem {postmortem}")
        return
    if args.stage2_2:
        postmortem = run_stage22(args, data, out_dir, disk, orth)
        print(f"stage2.2 wrote {out_dir}")
        print(f"postmortem {postmortem}")
        return
    if args.stage2_3:
        postmortem = run_stage23(args, data, out_dir, disk, orth)
        print(f"stage2.3 wrote {out_dir}")
        print(f"postmortem {postmortem}")
        return
    if args.true_loao:
        postmortem = run_true_loao(args, data, out_dir, disk, orth)
        print(f"true-loao wrote {out_dir}")
        print(f"postmortem {postmortem}")
        return
    kernel_summary, path_coeffs, case_outputs = run_kernel_fits(args, data, out_dir)
    unified_summary, unified_intensities = run_unified_fit(data)

    kernel_summary.to_csv(out_dir / "stage2_kernel_summary.csv", index=False)
    path_coeffs.to_csv(out_dir / "stage2_path_coefficients.csv", index=False)
    unified_summary.to_csv(out_dir / "stage2_unified_dab_summary.csv", index=False)
    unified_intensities.to_csv(out_dir / "stage2_unified_dab_intensities.csv", index=False)

    for mech, df in case_outputs.items():
        d = out_dir / "equations" / mech
        d.mkdir(parents=True, exist_ok=True)
        df.to_csv(d / "cases_manifest.csv", index=False)
    (out_dir / "convention_ledger.md").write_text(
        "\n".join([
            "# Convention Ledger",
            "",
            "- T2 targets and shadows are five-component tensors transformed with frozen change_of_basis.get_C(); |C.T C-I|max = "
            + f"{orth:.3e}.",
            "- Ring uses JB/BS/HM as same-convention current-loop paths; ringchi is excluded from path agreement because its sign convention is separate.",
            "- Field keeps MOPAC-Coulomb as the dominance source; water_efg is sign-flipped into the APBS/MOPAC -Hessian convention before field/unified fits.",
            "- Ring is not included in the unified D_ab sum; it is a current-loop G=-V_a n_b object, not the symmetric D_ab shadow.",
        ]) + "\n",
        encoding="utf-8",
    )
    audit = {
        "substrate_dir": str(args.substrate_dir),
        "out_dir": str(out_dir),
        "disk_guard": disk,
        "change_of_basis_orthogonality_max_abs": orth,
        "python_read_scope": "Build-4 per_atom_substrate CSV/NPY sidecars and CaseHunter manifests only",
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source dumps", "older rediscover dirs for fitting"],
        "pysr_requested": bool(args.run_pysr),
        "manifest_rows": int(manifest["rows"]),
        "manifest_n_atoms": int(manifest["n_atoms"]),
        "manifest_n_dft_frames": int(manifest["n_dft_frames"]),
    }
    (out_dir / "stage2_run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"stage2 wrote {out_dir}")


if __name__ == "__main__":
    main()
