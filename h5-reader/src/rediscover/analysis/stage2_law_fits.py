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
STAGE21_OUT_DIR = Path("/tmp/rediscover-runs/2026-06-04-stage2_1-sweep")
MIN_DISK_FREE_BYTES = 20 * 1024**3
MAX_REDISCOVER_BYTES = 15 * 1024**3
T2 = 5
PYSR_SAMPLE = 3000
SWEEP_KEEP_FRACTIONS = (0.10, 0.20, 0.35, 0.50, 0.70, 1.00)
FRAME_COUNTS = (20, 50, 100, 200, 400, 660)


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
    ap.add_argument("--stage2-1", action="store_true", help="run the Stage 2.1 happy-spot sweep and frame-count ablation")
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
    for case_dir in ["charge_wide", "mc", "field", "hbond"]:
        case_set, _cases = case_atoms(data["substrate_dir"], case_dir)
        atoms |= case_set
    dom = arrays["per_atom_substrate_features_dominance"]
    dom_cols = specs_for_array(specs, "per_atom_substrate_features_dominance")
    dab_max = np.maximum.reduce([
        np.asarray(dom[:, dom_cols[f"dominant_fraction_{m}"]], dtype=float)
        for m in ["charge", "mc", "field", "hbond"]
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
        for m in ["charge", "mc", "field", "hbond"]
    ])
    with np.errstate(all="ignore"):
        active = np.nanargmax(np.where(np.isfinite(dom_stack), dom_stack, -np.inf), axis=1)
    gap_stack = np.column_stack([
        np.asarray(dom[:, dom_cols["gap_to_2nd_charge_r"]], dtype=float),
        np.asarray(dom[:, dom_cols["gap_to_2nd_bond_r"]], dtype=float),
        np.asarray(dom[:, dom_cols["gap_to_2nd_field_r"]], dtype=float),
        np.asarray(dom[:, dom_cols["gap_to_2nd_hbond_r"]], dtype=float),
    ])
    mod_cols = specs_for_array(specs, "per_atom_substrate_driver_modulation_by_atom")
    mod = arrays["per_atom_substrate_driver_modulation_by_atom"]
    atom_index = rows["atom_index"].to_numpy(int)
    mod_stack = np.column_stack([
        np.asarray(mod[atom_index, mod_cols["sd_charge_T2_by_atom"]], dtype=float),
        np.asarray(mod[atom_index, mod_cols["sd_mc_lit_T2_by_atom"]], dtype=float),
        np.asarray(mod[atom_index, mod_cols["sd_mopac_coulomb_T2_by_atom"]], dtype=float),
        np.full(len(rows), np.nan, dtype=float),
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
            "hbond": "gap_to_2nd_hbond_r",
        }[kernel.name]
        return np.asarray(dom[:, dom_cols[gap_col]], dtype=float), gap_col, "higher_is_cleaner"
    if axis == "modulation":
        mod_name = {
            "ring": "sd_ring_jb_T2_by_atom",
            "charge": "sd_charge_T2_by_atom",
            "mc": "sd_mc_lit_T2_by_atom",
            "field": "sd_mopac_coulomb_T2_by_atom",
            "hbond": "",
        }[kernel.name]
        if not mod_name:
            return None, "not_emitted_for_hbond", "unavailable_in_build4"
        mod_cols = specs_for_array(specs, "per_atom_substrate_driver_modulation_by_atom")
        mod = arrays["per_atom_substrate_driver_modulation_by_atom"]
        atom_index = rows["atom_index"].to_numpy(int)
        return np.asarray(mod[atom_index, mod_cols[mod_name]], dtype=float), mod_name, "higher_driver_exercise"
    raise ValueError(axis)


def unified_axis_values(data: dict[str, object], axis: str) -> tuple[np.ndarray | None, str, str]:
    dab_max, active_gap, active_mod = active_dab_arrays(data)
    if axis == "dominance":
        return dab_max, "max_dominant_fraction_charge_mc_field_hbond", "higher_is_cleaner"
    if axis == "isolation":
        return active_gap, "gap_to_2nd_of_active_Dab_mechanism", "higher_is_cleaner"
    if axis == "modulation":
        return active_mod, "sd_of_active_Dab_mechanism_no_hbond_sd", "partial_no_hbond_modulation_column"
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
    loao = loao_modulation_fit(x, y, atoms)
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
    loao_pred = np.full_like(y, np.nan, dtype=float)
    loao_tgt = np.full_like(y, np.nan, dtype=float)
    for held in np.unique(atoms):
        tr = atoms != held
        te = atoms == held
        if tr.sum() < len(labels) + 3:
            continue
        cf = [atom_center(f[tr], atoms[tr]) for f in features]
        cy = atom_center(y[tr], atoms[tr])
        b, _ = fit_linear_multi(cf, cy, np.ones(tr.sum(), dtype=bool))
        if not np.isfinite(b).all():
            continue
        hf = [f[te] - f[te].mean(axis=0, keepdims=True) for f in features]
        loao_pred[te] = (np.column_stack([f.reshape(-1) for f in hf]) @ b).reshape((-1, T2))
        loao_tgt[te] = y[te] - y[te].mean(axis=0, keepdims=True)
    ok = np.isfinite(loao_pred).all(axis=1) & np.isfinite(loao_tgt).all(axis=1)
    loao_r2 = r2_score(loao_pred[ok], loao_tgt[ok])
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
    models = ["charge", "field", "ring", "mc", "hbond", "unified_Dab_sum"]
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
    models = ["charge", "field", "ring", "mc", "hbond", "unified_Dab_sum"]
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
        if model in {"mc", "hbond"}:
            rescue = f"; rescue {rescue_text(shape_summary, model)}"
        lines.append(f"- {model}: {shape_triplet(shape_summary, model)}; clean POP {clean_pop_text(shape_summary, model)}{rescue}{support_text}.")
    lines += [
        "McConnell/H-bond clean-end rescue criterion: primary coefficient CI off zero at the strict end of an input-side axis; no DFT fit is used to define the axis.",
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
        "hbond_driver_modulation_axis": "not emitted in per_atom_substrate_driver_modulation_by_atom; not fabricated",
        "unified_driver_modulation_axis": "active charge/mc/field emitted sd column; rows whose active mechanism is hbond have missing modulation",
        "geometric_noise_axis_distinct_from_driver_modulation": "not available in build4",
        "plots": plots,
    }
    (out_dir / "stage2_1_run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return write_stage21_postmortem(out_dir, args.substrate_dir, orth, shape_summary, frame_summary, disk, plots)


def main() -> None:
    args = parse_args()
    args.substrate_dir = args.substrate_dir.resolve()
    if args.stage2_1 and args.out_dir == OUT_DIR:
        args.out_dir = STAGE21_OUT_DIR
    out_dir = args.out_dir.resolve()
    disk = disk_guard(out_dir, drop_existing=bool(args.stage2_1))
    data, orth = load_stage2_inputs(args.substrate_dir, include_stage21=bool(args.stage2_1))
    manifest: dict[str, object] = data["manifest"]
    if args.stage2_1:
        postmortem = run_stage21(args, data, out_dir, disk, orth)
        print(f"stage2.1 wrote {out_dir}")
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
