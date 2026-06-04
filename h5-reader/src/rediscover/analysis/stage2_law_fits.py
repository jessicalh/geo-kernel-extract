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
    centered_by_train_atom,
    ci95,
    effective_n_components,
    finite_fmt,
    json_sanitize,
    r2_score,
    split_frame_block,
)


SUBSTRATE_DIR = Path("/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4")
OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-stage2-fits")
MIN_DISK_FREE_BYTES = 20 * 1024**3
MAX_REDISCOVER_BYTES = 15 * 1024**3
T2 = 5
PYSR_SAMPLE = 3000


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
            ShadowSpec("mc_category_sum", "per_atom_substrate_features_method_paths", "mc_category_T2", groups=5),
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
    KernelSpec(
        "hbond",
        "hbond",
        "hbond",
        ShadowSpec("hbond_geometric", "per_atom_substrate_features_classical", "hbond_T2"),
        (
            ShadowSpec("hbond_geometric", "per_atom_substrate_features_classical", "hbond_T2"),
            ShadowSpec("larsen_ppm_sum", "per_atom_substrate_features_hbond_conditioning", "larsen_hbond", groups=4),
        ),
        "CaseHunter hbond atoms plus dominant_fraction_hbond >= 0.5",
    ),
]


UNIFIED_SPECS = [
    ShadowSpec("charge_total", "per_atom_substrate_features_classical", "charge_q_over_r3_T2"),
    *[ShadowSpec(f"mc_category_{i}", "per_atom_substrate_features_method_paths", f"mc_category_T2_{5*i}", groups=0) for i in range(5)],
    ShadowSpec("mopac_field_backbone", "per_atom_substrate_features_method_paths", "mopac_coulomb_efg_backbone"),
    ShadowSpec("mopac_field_aromatic", "per_atom_substrate_features_method_paths", "mopac_coulomb_efg_aromatic"),
    ShadowSpec("water_field_efg_reconciled", "per_atom_substrate_features_method_paths", "water_efg", sign=-1.0),
    ShadowSpec("hbond_geometric", "per_atom_substrate_features_classical", "hbond_T2"),
    *[ShadowSpec(f"pq_type_{i}", "per_atom_substrate_features_ring_paths", f"pq_per_type_T2_{5*i}", groups=0) for i in range(8)],
    *[ShadowSpec(f"disp_type_{i}", "per_atom_substrate_features_ring_paths", f"disp_per_type_T2_{5*i}", groups=0) for i in range(8)],
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate-dir", type=Path, default=SUBSTRATE_DIR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
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


def disk_guard(out_dir: Path) -> dict[str, object]:
    out_dir = out_dir.resolve()
    runs_root = Path("/tmp/rediscover-runs").resolve()
    if not path_under(out_dir, runs_root):
        raise SystemExit(f"FATAL: stage2 output must be under {runs_root}, got {out_dir}")
    if path_under(out_dir, Path("/shared").resolve()):
        raise SystemExit(f"FATAL: refusing to write result data under /shared: {out_dir}")
    usage = shutil.disk_usage(out_dir.parent if out_dir.parent.exists() else runs_root)
    rediscover_bytes = directory_size_bytes(runs_root)
    if usage.free < MIN_DISK_FREE_BYTES:
        raise SystemExit(f"FATAL: disk guard: {usage.free / 1024**3:.1f} GiB free, need >=20 GiB")
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


def loao_modulation_fit(x: np.ndarray, y: np.ndarray, atoms: np.ndarray) -> dict[str, object]:
    labels = np.unique(atoms)
    pred = np.full_like(y, np.nan, dtype=float)
    target = np.full_like(y, np.nan, dtype=float)
    coefs = []
    for held in labels:
        train_atoms = atoms != held
        test_atoms = atoms == held
        if train_atoms.sum() < 3 or test_atoms.sum() < 2:
            continue
        x_train = atom_center(x[train_atoms], atoms[train_atoms])
        y_train = atom_center(y[train_atoms], atoms[train_atoms])
        k = fit_scalar_no_intercept(x_train, y_train)
        coefs.append(k)
        x_test = x[test_atoms] - x[test_atoms].mean(axis=0, keepdims=True)
        y_test = y[test_atoms] - y[test_atoms].mean(axis=0, keepdims=True)
        pred[test_atoms] = scalar_predict(x_test, k)
        target[test_atoms] = y_test
    ok = np.isfinite(pred).all(axis=1) & np.isfinite(target).all(axis=1)
    return {
        "loao_R2": r2_score(pred[ok], target[ok]),
        "loao_atoms": int(len(labels)),
        "loao_coeff_median": float(np.nanmedian(coefs)) if coefs else math.nan,
        "loao_coeff_min": float(np.nanmin(coefs)) if coefs else math.nan,
        "loao_coeff_max": float(np.nanmax(coefs)) if coefs else math.nan,
    }


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
            loao = loao_modulation_fit(x, y, atoms)
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
    for case_dir in ["charge_wide", "mc", "field", "hbond"]:
        case_set, _ = case_atoms(data["substrate_dir"], case_dir)
        atoms |= case_set
    dab_max = np.maximum.reduce([
        np.asarray(dom[:, dom_cols[f"dominant_fraction_{m}"]], dtype=float)
        for m in ["charge", "mc", "field", "hbond"]
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
    loao_pred = np.full_like(y, np.nan, dtype=float)
    loao_tgt = np.full_like(y, np.nan, dtype=float)
    for held in np.unique(atoms_s):
        tr = atoms_s != held
        te = atoms_s == held
        cf = [atom_center(f[tr], atoms_s[tr]) for f in features]
        cy = atom_center(y[tr], atoms_s[tr])
        b, _ = fit_linear_multi(cf, cy, np.ones(tr.sum(), dtype=bool))
        hf = [f[te] - f[te].mean(axis=0, keepdims=True) for f in features]
        loao_pred[te] = (np.column_stack([f.reshape(-1) for f in hf]) @ b).reshape((-1, T2))
        loao_tgt[te] = y[te] - y[te].mean(axis=0, keepdims=True)
    ok = np.isfinite(loao_pred).all(axis=1)
    loao_r2 = r2_score(loao_pred[ok], loao_tgt[ok])
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
        "selection_rule": "charge/mc/field/hbond CaseHunter atoms, max D_ab dominance >=0.5, ring dominance <0.7",
        "within_frameblock_R2": within_r2,
        "LOAO_R2": loao_r2,
        "terms": int(len(labels)),
        "test_frames": len(split["test_frames"]),
        "purged_train_frames": len(split["purged_frames"]),
        "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
    }])
    return summary, pd.DataFrame(inten_rows)


def main() -> None:
    args = parse_args()
    args.substrate_dir = args.substrate_dir.resolve()
    out_dir = args.out_dir.resolve()
    disk = disk_guard(out_dir)
    specs = load_json(args.substrate_dir / "per_atom_substrate_column_specs.json")
    manifest = load_json(args.substrate_dir / "per_atom_substrate_manifest.json")
    rows = pd.read_csv(args.substrate_dir / "per_atom_substrate_rows.csv")
    arrays = {
        name: load_array(args.substrate_dir, name)
        for name in [
            "per_atom_substrate_features_classical",
            "per_atom_substrate_features_method_paths",
            "per_atom_substrate_features_hbond_conditioning",
            "per_atom_substrate_features_ring_paths",
            "per_atom_substrate_features_dominance",
            "per_atom_substrate_dominance_bins",
        ]
    }
    target_lib = load_array(args.substrate_dir, "per_atom_substrate_target_T2")
    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(T2)).max())
    if orth >= 1.0e-12:
        raise SystemExit(f"FATAL: change_of_basis.get_C orthogonality failed: {orth:.3e}")
    data = {
        "manifest": manifest,
        "specs": specs,
        "rows": rows,
        "arrays": arrays,
        "target_lib": target_lib,
        "C": C,
        "substrate_dir": args.substrate_dir,
    }
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
            "- H-bond unified fit uses the geometric hbond_T2 shadow only; Larsen ppm tensors are reported as a separate per-kernel path and not double-counted.",
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
