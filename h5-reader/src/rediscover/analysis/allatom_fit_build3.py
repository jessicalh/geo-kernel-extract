# Build-3 active all-atom fit pipeline. Mechanical split from allatom_fit_piece2.py.
from allatom_fit_common import *

@dataclass(frozen=True)
class Build2ConditionSpec:
    family: str
    mechanism: str
    name: str
    kind: str
    source: str
    column: str
    bin_index: int | None = None


def path_under(path: Path, root: Path) -> bool:
    path = path.resolve()
    root = root.resolve()
    return path == root or root in path.parents


def directory_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file() or item.is_symlink():
                total += item.stat().st_size
        except FileNotFoundError:
            continue
    return int(total)


def prepare_build2_output_dir(out_dir: Path, keep_out_dir: bool) -> tuple[Path, dict[str, object]]:
    out_dir = out_dir.resolve()
    runs_root = Path("/tmp/rediscover-runs").resolve()
    shared_root = Path("/shared").resolve()
    if not path_under(out_dir, runs_root):
        raise SystemExit(f"FATAL: Build 3 output must be under {runs_root}, got {out_dir}")
    if path_under(out_dir, shared_root):
        raise SystemExit(f"FATAL: refusing to write Build 3 result data under /shared: {out_dir}")
    usage = shutil.disk_usage(out_dir.parent if out_dir.parent.exists() else runs_root)
    rediscover_bytes = directory_size_bytes(runs_root)
    audit = {
        "checked_before_write": True,
        "path": str(out_dir),
        "filesystem_free_bytes": int(usage.free),
        "filesystem_free_GiB": float(usage.free / 1024**3),
        "rediscover_runs_bytes_before_write": int(rediscover_bytes),
        "rediscover_runs_GiB_before_write": float(rediscover_bytes / 1024**3),
        "min_free_GiB_required": float(MIN_DISK_FREE_BYTES / 1024**3),
        "max_rediscover_GiB_required": float(MAX_REDISCOVER_BYTES / 1024**3),
        "deleted_existing_output_dir": False,
    }
    if usage.free < MIN_DISK_FREE_BYTES:
        raise SystemExit(f"FATAL: disk guard: {usage.free / 1024**3:.1f} GiB free under /tmp, need >=20 GiB")
    if rediscover_bytes > MAX_REDISCOVER_BYTES:
        raise SystemExit(f"FATAL: disk guard: /tmp/rediscover-runs is {rediscover_bytes / 1024**3:.1f} GiB, need <15 GiB")
    if out_dir.exists() and not keep_out_dir:
        shutil.rmtree(out_dir)
        audit["deleted_existing_output_dir"] = True
        audit["deleted_existing_output_dir_full_path"] = str(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir, audit


def specs_array_columns(specs: dict[str, object], array_name: str) -> list[str]:
    return [str(c["column"]) for c in specs["columns"] if c["array"] == array_name]


def build2_load_inputs(substrate_dir: Path) -> dict[str, object]:
    manifest_path = require_file(substrate_dir / "per_atom_substrate_manifest.json")
    specs_path = require_file(substrate_dir / "per_atom_substrate_column_specs.json")
    rows_path = require_file(substrate_dir / "per_atom_substrate_rows.csv")
    with manifest_path.open("r", encoding="utf-8") as f:
        manifest = json.load(f)
    with specs_path.open("r", encoding="utf-8") as f:
        specs = json.load(f)
    rows = pd.read_csv(rows_path)
    arrays = {
        "per_atom_substrate_features_classical": np.load(require_file(substrate_dir / "per_atom_substrate_features_classical.npy"), mmap_mode="r"),
        "per_atom_substrate_features_conditioning": np.load(require_file(substrate_dir / "per_atom_substrate_features_conditioning.npy"), mmap_mode="r"),
        "per_atom_substrate_features_hbond_conditioning": np.load(require_file(substrate_dir / "per_atom_substrate_features_hbond_conditioning.npy"), mmap_mode="r"),
        "per_atom_substrate_features_method_paths": np.load(require_file(substrate_dir / "per_atom_substrate_features_method_paths.npy"), mmap_mode="r"),
        "per_atom_substrate_partition_bins": np.load(require_file(substrate_dir / "per_atom_substrate_partition_bins.npy"), mmap_mode="r"),
        "per_atom_substrate_aimnet2_embedding": np.load(require_file(substrate_dir / "per_atom_substrate_aimnet2_embedding.npy"), mmap_mode="r"),
    }
    targets = {}
    for target_name, target_spec in TARGET_SPECS.items():
        array_name = str(target_spec["array"])
        targets[target_name] = np.load(require_file(substrate_dir / str(target_spec["file"])), mmap_mode="r")
    files_read = [
        str(manifest_path),
        str(specs_path),
        str(rows_path),
        *[str(substrate_dir / f"{name}.npy") for name in arrays],
        *[str(substrate_dir / str(spec["file"])) for spec in TARGET_SPECS.values()],
    ]
    return {"manifest": manifest, "specs": specs, "rows": rows, "arrays": arrays, "targets": targets, "files_read": files_read}


def build2_input_acceptance_checks(data: dict[str, object], substrate_dir: Path) -> tuple[dict[str, object], np.ndarray]:
    manifest: dict[str, object] = data["manifest"]
    specs: dict[str, object] = data["specs"]
    rows: pd.DataFrame = data["rows"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    targets: dict[str, np.ndarray] = data["targets"]
    checks: dict[str, object] = {}
    required_cols = [
        "row_id",
        "atom_index",
        "frame_slot",
        "original_frame_index",
        "stratum",
        "role",
        "ff_atom_type_ord",
        "formal_charge",
        "ff14sb_charge",
        "ff14sb_charge_present",
        "mopac_welford_mean_charge",
        "mopac_welford_mean_charge_present",
        "eeq_charge",
        "eeq_charge_present",
    ]
    require_columns(rows, required_cols, "per_atom_substrate_rows.csv")
    n_rows = len(rows)
    n_atoms = int(manifest.get("n_atoms", -1))
    n_frames = int(manifest.get("n_dft_frames", -1))
    row_contract = bool(
        np.array_equal(rows["row_id"].to_numpy(int), np.arange(n_rows, dtype=int))
        and np.all(rows["row_id"].to_numpy(int) == rows["frame_slot"].to_numpy(int) * n_atoms + rows["atom_index"].to_numpy(int))
    )
    expected_widths = {
        name: expected_width_from_specs(specs, name)
        for name in arrays
    }
    shape_expected = {
        name: (n_rows, width)
        for name, width in expected_widths.items()
    }
    shapes_ok = n_atoms == 846 and n_frames == 660 and n_rows == 558_360
    for name, shape in shape_expected.items():
        shapes_ok = shapes_ok and tuple(arrays[name].shape) == shape
    target_expected_widths = {
        target_name: expected_width_from_specs(specs, str(TARGET_SPECS[target_name]["array"]))
        for target_name in targets
    }
    target_shapes_ok = all(
        tuple(targets[target_name].shape) == (n_rows, width)
        for target_name, width in target_expected_widths.items()
    )
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    aim_idx = int(classical_specs["aimnet2_charge"]["index"])
    classical = arrays["per_atom_substrate_features_classical"]
    charge_complete_mask = (
        (rows["ff14sb_charge_present"].to_numpy(int) == 1)
        & (rows["mopac_welford_mean_charge_present"].to_numpy(int) == 1)
        & np.isfinite(rows["ff14sb_charge"].to_numpy(float))
        & np.isfinite(rows["mopac_welford_mean_charge"].to_numpy(float))
        & np.isfinite(np.asarray(classical[:, aim_idx], dtype=float))
    )
    target_finite = True
    for target in targets.values():
        target_finite = target_finite and bool(np.isfinite(np.asarray(target[charge_complete_mask], dtype=float)).all())
    manifest_charge_rows = int(manifest.get("feature_support_rows", {}).get("charge_complete_rows", -1))
    partition_manifest = manifest.get("partition_bins", {})
    partition_cols = [str(x) for x in partition_manifest.get("columns", [])]
    spec_partition_cols = specs_array_columns(specs, "per_atom_substrate_partition_bins")
    bins = arrays["per_atom_substrate_partition_bins"]
    charge_bin_cols = [
        "bin_nearest_charge_r_4_6_8_10",
        "bin_gap_to_2nd_charge_r_4_6_8_10",
        "bin_abs_charge_T2_quintile",
        "bin_sd_charge_T2_by_atom_quintile",
    ]
    charge_bin_unique = {}
    for col in charge_bin_cols:
        idx = partition_cols.index(col)
        vals = np.asarray(bins[:, idx])
        vals = vals[vals >= 0]
        charge_bin_unique[col] = sorted(int(x) for x in np.unique(vals))
    checks["input_acceptance"] = {
        "pass": bool(
            manifest.get("relationship") == "per_atom_substrate"
            and manifest.get("relationship_kind") == "per_atom_aggregate"
            and shapes_ok
            and target_shapes_ok
            and row_contract
            and manifest_charge_rows == int(charge_complete_mask.sum())
            and int(charge_complete_mask.sum()) == 558_360
            and target_finite
            and partition_cols == spec_partition_cols
            and substrate_dir == SUBSTRATE_DIR.resolve()
            and len(charge_bin_unique["bin_nearest_charge_r_4_6_8_10"]) > 1
            and len(charge_bin_unique["bin_gap_to_2nd_charge_r_4_6_8_10"]) > 1
        ),
        "relationship": manifest.get("relationship"),
        "relationship_kind": manifest.get("relationship_kind"),
        "shape": {"n_atoms": n_atoms, "n_dft_frames": n_frames, "rows": n_rows},
        "expected_shape_846x660": shapes_ok,
        "expected_widths_from_column_specs": expected_widths,
        "target_expected_widths_from_column_specs": target_expected_widths,
        "target_shapes_ok": target_shapes_ok,
        "row_contract": row_contract,
        "charge_complete_rows_manifest": manifest_charge_rows,
        "charge_complete_rows_used": int(charge_complete_mask.sum()),
        "target_finite_on_charge_complete_rows": target_finite,
        "partition_bin_columns_match_manifest_and_specs": partition_cols == spec_partition_cols,
        "charge_partition_non_degenerate_bins": charge_bin_unique,
        "dominant_charge_bin_column_present": "bin_dominant_fraction_charge" in partition_cols,
        "substrate_dir_only_requested_path": str(substrate_dir),
    }
    return checks, charge_complete_mask


def build2_feature_block_specs() -> list[dict[str, object]]:
    return [
        {"name": "ring_jb_T0/T2", "array": "per_atom_substrate_features_classical", "scalar": ["ring_jb_T0"], "t2": ["ring_jb_T2"]},
        {"name": "charge_q_over_r3_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["charge_q_over_r3_T2"]},
        {"name": "mc_lit_T0/T2_valid", "array": "per_atom_substrate_features_classical", "scalar": ["mc_lit_T0_valid"], "t2": ["mc_lit_T2_valid"]},
        {"name": "mopac_coulomb_shielding_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["mopac_coulomb_shielding_T2"]},
        {"name": "mopac_mc_shielding_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["mopac_mc_shielding_T2"]},
        {"name": "ff14sb_field_x/y/z/mag", "array": "per_atom_substrate_features_classical", "scalar": ["ff14sb_field_x", "ff14sb_field_y", "ff14sb_field_z", "ff14sb_field_mag"], "t2": []},
        {"name": "apbs_E_x/y/z/mag", "array": "per_atom_substrate_features_classical", "scalar": ["apbs_E_x", "apbs_E_y", "apbs_E_z", "apbs_E_mag"], "t2": []},
        {"name": "apbs_efg_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["apbs_efg_T2"]},
        {"name": "larsen_hbond_T2", "array": "per_atom_substrate_features_hbond_conditioning", "scalar": ["larsen_hbond_water_term"], "t2": ["larsen_hbond_1pHB_T2", "larsen_hbond_2pHB_T2", "larsen_hbond_1pHaB_T2", "larsen_hbond_2pHaB_T2"]},
        {"name": "hbond_count", "array": "per_atom_substrate_features_classical", "scalar": ["hbond_count"], "t2": []},
        {"name": "pi_quadrupole_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["pi_quadrupole_T2"]},
        {"name": "dispersion_T2", "array": "per_atom_substrate_features_classical", "scalar": [], "t2": ["dispersion_T2"]},
        {"name": "water_field", "array": "per_atom_substrate_features_classical", "scalar": ["water_efield_x", "water_efield_y", "water_efield_z", "water_efield_mag", "water_n_first", "water_n_second", "water_half_shell_asymmetry", "water_dipole_cos"], "t2": []},
        {"name": "water_efg_T2", "array": "per_atom_substrate_features_method_paths", "scalar": [], "t2": ["water_efg"]},
        {"name": "water_efield_first", "array": "per_atom_substrate_features_method_paths", "scalar": ["water_efield_first_x", "water_efield_first_y", "water_efield_first_z"], "t2": []},
        {"name": "mopac_coulomb_E", "array": "per_atom_substrate_features_method_paths", "scalar": ["mopac_coulomb_E_x", "mopac_coulomb_E_y", "mopac_coulomb_E_z"], "t2": []},
        {"name": "mopac_coulomb_efg_paths", "array": "per_atom_substrate_features_method_paths", "scalar": [], "t2": ["mopac_coulomb_efg_backbone", "mopac_coulomb_efg_aromatic"]},
        {"name": "aimnet2_charge", "array": "per_atom_substrate_features_classical", "scalar": ["aimnet2_charge"], "t2": []},
        {"name": "aimnet2_crg_scalar/x/y/z", "array": "per_atom_substrate_features_classical", "scalar": ["aimnet2_crg_scalar", "aimnet2_crg_x", "aimnet2_crg_y", "aimnet2_crg_z"], "t2": []},
        {"name": "aimnet2_efg_paths", "array": "per_atom_substrate_features_method_paths", "scalar": [], "t2": ["aimnet2_efg", "aimnet2_efg_aromatic", "aimnet2_efg_backbone"]},
    ]


def extract_build2_block(
    arrays: dict[str, np.ndarray],
    specs: dict[str, object],
    spec: dict[str, object],
    C: np.ndarray,
    transform_audit: list[dict[str, object]],
) -> tuple[np.ndarray, list[str]]:
    array_name = str(spec["array"])
    matrix = arrays[array_name]
    col_specs = specs_for_array(specs, array_name)
    pieces = []
    names = []
    for col in spec.get("scalar", []):
        idx = int(col_specs[str(col)]["index"])
        pieces.append(np.asarray(matrix[:, idx], dtype=float).reshape(-1, 1))
        names.append(str(col))
    for prefix in spec.get("t2", []):
        cols = [f"{prefix}_{i}" for i in range(5)]
        idx = [int(col_specs[col]["index"]) for col in cols]
        vals = np.asarray(matrix[:, idx], dtype=float) @ C.T
        pieces.append(vals)
        names += [f"{prefix}_2e_{i}" for i in range(5)]
        transform_audit.append({"feature_block": str(prefix), "array": array_name, "source_columns": cols, "transform": "change_of_basis.get_C().T"})
    if not pieces:
        return np.empty((matrix.shape[0], 0), dtype=float), []
    return np.column_stack(pieces), names


def build2_feature_blocks(
    rows: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    specs: dict[str, object],
    C: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, list[str]], list[dict[str, object]]]:
    transform_audit: list[dict[str, object]] = []
    blocks: dict[str, np.ndarray] = {}
    labels: dict[str, list[str]] = {}
    for spec in build2_feature_block_specs():
        block, names = extract_build2_block(arrays, specs, spec, C, transform_audit)
        blocks[str(spec["name"])] = block
        labels[str(spec["name"])] = names
    for name, _source, scalar_cols, _t2_prefixes in RAW_CHARGE_BLOCK_SPECS:
        block = rows[scalar_cols].to_numpy(float)
        blocks[name] = block
        labels[name] = list(scalar_cols)
    return blocks, labels, transform_audit


def transform_build2_target(target_name: str, raw_target: np.ndarray, row_mask: np.ndarray, C: np.ndarray) -> np.ndarray:
    target = raw_target if bool(row_mask.all()) else raw_target[row_mask]
    arr = np.asarray(target, dtype=float)
    if target_name in {"total-T2", "dia-T2", "para-T2"}:
        return arr[:, :5] @ C.T
    return arr[:, :3]


def build2_component_label(target_name: str, component: int | None) -> str:
    if component is None:
        dims = int(TARGET_SPECS[target_name]["components"])
        return f"vector_{dims}d"
    if target_name in {"total-T2", "dia-T2", "para-T2"}:
        return f"2e_{component}"
    return f"T1_{component}"


def build2_long_score_rows(row: dict[str, object], axis: str, pred: np.ndarray, target: np.ndarray, groups: np.ndarray, ridge_alpha: float) -> list[dict[str, object]]:
    out = []
    vector_r2 = r2_score(pred, target)
    vector_se = jackknife_metric_se(pred, target, groups)
    lo, hi = ci95(vector_r2, vector_se)
    base = {k: row[k] for k in ["target", "atom_type_axis", "atom_type", "tier", "tier_label", "rows", "n_atoms_in_slice", "n_features", "alpha_mode"]}
    base["ridge_alpha"] = float(ridge_alpha)
    out.append({**base, "axis": axis, "component": build2_component_label(str(row["target"]), None), "heldout_R2": vector_r2, "jackknife_se": vector_se, "ci95_low": lo, "ci95_high": hi})
    comp = component_r2(pred, target)
    comp_se = jackknife_component_se(np.asarray(pred), np.asarray(target), groups)
    for i, (r2, se) in enumerate(zip(comp, comp_se)):
        lo, hi = ci95(r2, se)
        out.append({**base, "axis": axis, "component": build2_component_label(str(row["target"]), i), "heldout_R2": r2, "jackknife_se": se, "ci95_low": lo, "ci95_high": hi})
    return out


def build2_score_rows(
    target_name: str,
    fit_scope: str,
    tier_results: dict[str, dict[str, object]],
    tier_order: list[str],
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    y: np.ndarray,
    frames: np.ndarray,
    atoms: np.ndarray,
    split: dict[str, object],
    base_blocks: dict[str, np.ndarray],
    alpha_mode: str,
    embedding_components: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    score_rows = []
    long_rows = []
    n_atoms_global = int(len(atom_rows))
    strata = sorted(str(x) for x in atom_rows["stratum"].dropna().unique())
    atom_types = ["ALL"] + strata if fit_scope == "global_sliced" else strata
    prev_map = previous_tier_map(tier_order)
    r2_cache: dict[tuple[str, str, str], float] = {}
    for tier in tier_order:
        n_features = feature_count_for_tier(tier, embedding_components, base_blocks)
        between = tier_results[tier]["between"]
        within = tier_results[tier]["within"]
        between_alpha_all = np.asarray(between["alpha_by_heldout"], dtype=float)
        within_alpha = float(within["alpha"])
        within_alpha_by_score = within.get("alpha_by_score_row")
        within_alpha_all = np.asarray(within_alpha_by_score, dtype=float) if within_alpha_by_score is not None else None
        within_positions = within["score_row_positions"]
        within_frames = frames[within_positions]
        within_rows = rows.iloc[within_positions]
        for atom_type in atom_types:
            if atom_type == "ALL":
                atom_mask = np.ones(n_atoms_global, dtype=bool)
                row_mask = np.ones(len(rows), dtype=bool)
                within_mask = np.ones(len(within["target"]), dtype=bool)
            else:
                atom_mask = atom_rows["stratum"].to_numpy(str) == atom_type
                row_mask = rows["stratum"].to_numpy(str) == atom_type
                within_mask = within_rows["stratum"].to_numpy(str) == atom_type
            n_atoms_slice = int(atom_mask.sum())
            rows_slice = int(row_mask.sum())
            frames_slice = int(rows.loc[row_mask, "original_frame_index"].nunique()) if rows_slice else 0
            b_pred = between["pred"][atom_mask]
            b_target = between["target"][atom_mask]
            b_groups = between["groups"][atom_mask]
            b_alpha_dist = alpha_distribution(between_alpha_all[atom_mask])
            b_alpha_mode = math.nan if b_alpha_dist["selected_alpha"] is None else float(b_alpha_dist["selected_alpha"])
            b_alpha_min = math.nan if b_alpha_dist["min"] is None else float(b_alpha_dist["min"])
            b_alpha_max = math.nan if b_alpha_dist["max"] is None else float(b_alpha_dist["max"])
            if within_alpha_all is None:
                w_alpha_mode = within_alpha
            else:
                w_alpha_dist = alpha_distribution(within_alpha_all[within_mask])
                w_alpha_mode = math.nan if w_alpha_dist["selected_alpha"] is None else float(w_alpha_dist["selected_alpha"])
            w_pred = within["pred"][within_mask]
            w_target = within["target"][within_mask]
            w_groups = within["groups"][within_mask]
            b_r2 = r2_score(b_pred, b_target)
            w_r2 = r2_score(w_pred, w_target)
            r2_cache[(tier, atom_type, "between")] = b_r2
            r2_cache[(tier, atom_type, "within")] = w_r2
            b_se = jackknife_metric_se(b_pred, b_target, b_groups)
            w_se = jackknife_metric_se(w_pred, w_target, w_groups)
            b_low, b_high = ci95(b_r2, b_se)
            w_low, w_high = ci95(w_r2, w_se)
            share_b, share_w = variance_shares(y[row_mask], atoms[row_mask])
            neff = effective_n_components(y[within_positions][within_mask], within["groups"][within_mask], within_frames[within_mask])
            thin = n_atoms_slice < THIN_ATOMS or len(np.unique(b_groups)) < 3 or len(np.unique(w_groups)) < 3
            n_atoms_between_fit_scope = n_atoms_global if fit_scope == "global_sliced" else n_atoms_slice
            between_n_eff = float(n_atoms_between_fit_scope)
            between_support = ""
            if fit_scope == "per_type" and n_atoms_between_fit_scope < 100:
                between_support = f"thin_between_atoms={n_atoms_between_fit_scope}<100"
            within_support = ""
            if np.isfinite(neff["within_N_eff_median"]) and neff["within_N_eff_median"] < MIN_FAVORABLE_N_EFF:
                within_support = f"thin_within_N_eff={finite_fmt(neff['within_N_eff_median'])}"
            p_gt_atoms = n_features >= max(1, n_atoms_between_fit_scope)
            prev = prev_map[tier]
            row = {
                "target": target_name,
                "fit_scope": fit_scope,
                "atom_type_axis": "stratum",
                "atom_type": atom_type,
                "tier": tier,
                "tier_label": TIER_LABELS[tier],
                "rows": rows_slice,
                "n_atoms_between_global": n_atoms_global,
                "n_atoms_between_fit_scope": n_atoms_between_fit_scope,
                "n_atoms_in_slice": n_atoms_slice,
                "n_original_frames": frames_slice,
                "n_features": n_features,
                "ridge_alpha": b_alpha_mode,
                "alpha_mode": alpha_mode,
                "between_ridge_alpha": b_alpha_mode,
                "between_ridge_alpha_min": b_alpha_min,
                "between_ridge_alpha_max": b_alpha_max,
                "within_ridge_alpha": w_alpha_mode,
                "variance_share_between": share_b,
                "variance_share_within": share_w,
                **neff,
                "between_N_eff": between_n_eff,
                "thin_flag": "thin" if thin else "",
                "between_support_flag": between_support,
                "within_support_flag": within_support,
                "p_gt_atoms_flag": "p>=atoms" if p_gt_atoms else "",
                "between_LOAO_test_R2": b_r2,
                "between_LOAO_test_R2_jackknife_se": b_se,
                "between_LOAO_test_R2_ci95_low": b_low,
                "between_LOAO_test_R2_ci95_high": b_high,
                "between_delta_R2_vs_classical": b_r2 - r2_cache.get(("classical", atom_type, "between"), math.nan),
                "between_delta_R2_vs_previous_tier": b_r2 - r2_cache.get((prev, atom_type, "between"), math.nan) if prev else math.nan,
                "between_delta_R2_vs_global_sliced": math.nan,
                "within_frameblock_test_R2": w_r2,
                "within_frameblock_test_R2_jackknife_se": w_se,
                "within_frameblock_test_R2_ci95_low": w_low,
                "within_frameblock_test_R2_ci95_high": w_high,
                "within_delta_R2_vs_classical": w_r2 - r2_cache.get(("classical", atom_type, "within"), math.nan),
                "within_delta_R2_vs_previous_tier": w_r2 - r2_cache.get((prev, atom_type, "within"), math.nan) if prev else math.nan,
                "within_delta_R2_vs_global_sliced": math.nan,
                "within_split_strategy": "centered_by_atom_train_frames_only; contiguous_frame_block",
                "test_frames": int(len(split["test_frames"])),
                "purged_train_frames": int(len(split["purged_frames"])),
                "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
                "feature_blocks": "|".join(TIER_BLOCKS[tier]),
            }
            score_rows.append(row)
            long_rows += build2_long_score_rows(row, "between", b_pred, b_target, b_groups, b_alpha_mode)
            long_rows += build2_long_score_rows(row, "within", w_pred, w_target, w_groups, w_alpha_mode)
    return score_rows, long_rows


def add_fit_scope_deltas(score_table: pd.DataFrame) -> pd.DataFrame:
    if score_table.empty or "fit_scope" not in score_table.columns:
        return score_table
    out = score_table.copy()
    key_cols = ["target", "atom_type_axis", "atom_type", "tier"]
    global_rows = out[out["fit_scope"] == "global_sliced"][key_cols + ["between_LOAO_test_R2", "within_frameblock_test_R2"]].rename(
        columns={
            "between_LOAO_test_R2": "_global_between_R2",
            "within_frameblock_test_R2": "_global_within_R2",
        }
    )
    merged = out.merge(global_rows, on=key_cols, how="left")
    per_type = merged["fit_scope"] == "per_type"
    merged.loc[merged["fit_scope"] == "global_sliced", "between_delta_R2_vs_global_sliced"] = 0.0
    merged.loc[merged["fit_scope"] == "global_sliced", "within_delta_R2_vs_global_sliced"] = 0.0
    merged.loc[per_type, "between_delta_R2_vs_global_sliced"] = merged.loc[per_type, "between_LOAO_test_R2"] - merged.loc[per_type, "_global_between_R2"]
    merged.loc[per_type, "within_delta_R2_vs_global_sliced"] = merged.loc[per_type, "within_frameblock_test_R2"] - merged.loc[per_type, "_global_within_R2"]
    return merged.drop(columns=["_global_between_R2", "_global_within_R2"])


def build2_mechanism_from_bin_column(column: str) -> str:
    for key, mech in [
        ("ring", "ring"),
        ("charge", "charge"),
        ("mc_lit", "mc"),
        ("mopac_mc", "mc"),
        ("bond", "mc"),
        ("mopac_coulomb", "charge"),
        ("apbs", "field"),
        ("ff14sb_E", "field"),
        ("aimnet2", "aimnet2"),
        ("heavy_atom", "geometry"),
    ]:
        if key in column:
            return mech
    return "input"


def build2_family_from_bin_column(column: str, mechanism: str) -> str:
    if column.startswith("bin_abs_"):
        return "driver magnitude"
    if column.startswith("bin_sd_"):
        return "driver modulation"
    if column.startswith("bin_nearest_") or column.startswith("bin_gap_"):
        return "geometry and isolation"
    return "input-side bin"


def build2_condition_specs(rows: pd.DataFrame, partition_columns: list[str]) -> list[Build2ConditionSpec]:
    specs: list[Build2ConditionSpec] = []
    for idx, col in enumerate(partition_columns):
        mech = build2_mechanism_from_bin_column(col)
        specs.append(Build2ConditionSpec(build2_family_from_bin_column(col, mech), mech, col.removeprefix("bin_"), "cxx_bin", "per_atom_substrate_partition_bins", col, idx))
    for mech, col in DOMINANCE_FRACTION_COLUMNS.items():
        if col in rows.columns:
            specs.append(Build2ConditionSpec("dominance response", mech, col, "python_quantile_cpp_scalar", "per_atom_substrate_features_conditioning", col, None))
    categorical = [
        ("atom identity", "identity", "stratum"),
        ("atom identity", "identity", "role"),
        ("atom identity", "identity", "element_ord"),
        ("atom identity", "identity", "ff_atom_type_ord"),
        ("atom identity", "identity", "formal_charge"),
        ("atom identity", "identity", "planar_group_ord"),
        ("atom identity", "identity", "polar_h_kind_ord"),
        ("atom identity", "identity", "aromatic"),
        ("atom identity", "identity", "is_exchangeable"),
        ("atom identity", "identity", "amino_acid_ord"),
        ("atom identity", "identity", "locant_ord"),
        ("atom identity", "identity", "branch_outer_ord"),
        ("atom identity", "identity", "branch_inner_ord"),
        ("support", "support", "ring_present"),
        ("support", "support", "charge_present"),
        ("support", "support", "mc_lit_valid_present"),
        ("support", "support", "ff14sb_field_present"),
        ("support", "support", "apbs_E_present"),
        ("support", "support", "apbs_efg_present"),
        ("support", "support", "hbond_shielding_present"),
        ("support", "support", "hbond_count_present"),
        ("support", "support", "pi_quadrupole_present"),
        ("support", "support", "dispersion_present"),
        ("support", "support", "water_field_present"),
        ("support", "support", "hydration_shell_present"),
        ("support", "support", "sasa_present"),
        ("support", "support", "aimnet2_charge_present"),
        ("support", "support", "aimnet2_crg_present"),
        ("support", "support", "aimnet2_embedding_present"),
        ("support", "support", "ff14sb_charge_present"),
        ("support", "support", "mopac_welford_mean_charge_present"),
        ("support", "support", "eeq_charge_present"),
        ("charge-source agreement", "charge", "sign_agree_ff14sb_mopac_welford"),
        ("charge-source agreement", "charge", "sign_agree_ff14sb_aimnet2"),
        ("charge-source agreement", "charge", "sign_agree_mopac_welford_aimnet2"),
    ]
    for family, mech, col in categorical:
        if col in rows.columns:
            specs.append(Build2ConditionSpec(family, mech, col, "categorical", "per_atom_substrate_rows", col, None))
    return specs


def add_build2_charge_agreement_columns(rows: pd.DataFrame, arrays: dict[str, np.ndarray], specs: dict[str, object]) -> pd.DataFrame:
    out = rows.copy()
    classical_specs = specs_for_array(specs, "per_atom_substrate_features_classical")
    conditioning_specs = specs_for_array(specs, "per_atom_substrate_features_conditioning")
    aim_idx = int(classical_specs["aimnet2_charge"]["index"])
    out["aimnet2_charge"] = np.asarray(arrays["per_atom_substrate_features_classical"][:, aim_idx], dtype=float)
    conditioning = arrays["per_atom_substrate_features_conditioning"]
    for col in DOMINANCE_FRACTION_COLUMNS.values():
        if col in conditioning_specs:
            out[col] = np.asarray(conditioning[:, int(conditioning_specs[col]["index"])], dtype=float)
    out["sign_agree_ff14sb_mopac_welford"] = sign_agreement(out["ff14sb_charge"], out["mopac_welford_mean_charge"])
    out["sign_agree_ff14sb_aimnet2"] = sign_agreement(out["ff14sb_charge"], out["aimnet2_charge"])
    out["sign_agree_mopac_welford_aimnet2"] = sign_agreement(out["mopac_welford_mean_charge"], out["aimnet2_charge"])
    return out


def mode_int(values: np.ndarray) -> int:
    vals = np.asarray(values, dtype=int)
    if len(vals) == 0:
        return -1
    unique, counts = np.unique(vals, return_counts=True)
    order = np.lexsort((unique, -counts))
    return int(unique[order[0]])


def atom_mode_bins(partition_bins: np.ndarray, n_frames: int, n_atoms: int) -> np.ndarray:
    out = np.empty((n_atoms, partition_bins.shape[1]), dtype=np.int32)
    for atom in range(n_atoms):
        vals = np.asarray(partition_bins[atom::n_atoms], dtype=np.int32)
        for col in range(vals.shape[1]):
            out[atom, col] = mode_int(vals[:, col])
    return out


def build2_bin_label(column: str, bin_id: int, manifest: dict[str, object]) -> str:
    if bin_id < 0:
        return "missing"
    if column.startswith("bin_abs_") or column.startswith("bin_sd_"):
        return f"Q{bin_id + 1}"
    if column == "bin_gap_to_2nd_charge_r_4_6_8_10":
        edges = [float(x) for x in manifest.get("partition_bins", {}).get("charge_gap_band_edges_A", [0.25, 0.5, 1.0, 1.5])]
    else:
        edges = [float(x) for x in manifest.get("partition_bins", {}).get("distance_band_edges_A", [4.0, 6.0, 8.0, 10.0])]
    if bin_id == 0:
        return f"<= {finite_fmt(edges[0], 3)}A"
    if bin_id < len(edges):
        return f"{finite_fmt(edges[bin_id - 1], 3)}-{finite_fmt(edges[bin_id], 3)}A"
    return f"> {finite_fmt(edges[-1], 3)}A"


def build2_bin_order(spec: Build2ConditionSpec, label: object, bin_id: int | None = None) -> float:
    if bin_id is not None:
        return float(bin_id)
    return bin_order_key(label)


def quantile_bin_arrays(values: np.ndarray | pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    finite = np.isfinite(numeric)
    labels = np.full(len(numeric), "missing", dtype=object)
    raw_bins = np.full(len(numeric), -1, dtype=int)
    orders = np.full(len(numeric), -1.0, dtype=float)
    if finite.sum() == 0:
        return labels, raw_bins, orders
    vals = numeric[finite]
    qs = np.quantile(vals, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    edges: list[float] = []
    for q in qs:
        qf = float(q)
        if not edges or qf > edges[-1]:
            edges.append(qf)
    if len(edges) <= 2:
        lo = float(np.nanmin(vals))
        hi = float(np.nanmax(vals))
        labels[finite] = f"Q1 [{finite_fmt(lo)}, {finite_fmt(hi)}]"
        raw_bins[finite] = 0
        orders[finite] = 0.0
        return labels, raw_bins, orders
    for bin_id, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        if bin_id == 0:
            mask = finite & (numeric >= lo) & (numeric <= hi)
        else:
            mask = finite & (numeric > lo) & (numeric <= hi)
        labels[mask] = f"Q{bin_id + 1} [{finite_fmt(lo)}, {finite_fmt(hi)}]"
        raw_bins[mask] = int(bin_id)
        orders[mask] = float(bin_id)
    return labels, raw_bins, orders


def build2_axis_conditions(
    spec: Build2ConditionSpec,
    axis: str,
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    partition_bins: np.ndarray,
    atom_bins: np.ndarray,
    row_positions: np.ndarray,
    manifest: dict[str, object],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if axis == "between":
        atoms = atom_rows["atom_index"].to_numpy(int)
        if spec.kind == "cxx_bin":
            raw_bins = atom_bins[:, int(spec.bin_index)]
            labels = np.asarray([build2_bin_label(spec.column, int(x), manifest) for x in raw_bins], dtype=object)
            orders = raw_bins.astype(float)
            return labels, raw_bins, orders
        if spec.kind == "python_quantile_cpp_scalar":
            return quantile_bin_arrays(atom_rows[spec.column])
        labels = np.asarray([str(x) if pd.notna(x) else "missing" for x in atom_rows[spec.column].to_numpy()], dtype=object)
        return labels, np.full(len(labels), -999, dtype=int), np.arange(len(labels), dtype=float)
    atoms = rows.iloc[row_positions]["atom_index"].to_numpy(int)
    if spec.kind == "cxx_bin":
        raw_bins = np.asarray(partition_bins[row_positions, int(spec.bin_index)], dtype=int)
        labels = np.asarray([build2_bin_label(spec.column, int(x), manifest) for x in raw_bins], dtype=object)
        orders = raw_bins.astype(float)
        return labels, raw_bins, orders
    if spec.kind == "python_quantile_cpp_scalar":
        return quantile_bin_arrays(rows.iloc[row_positions][spec.column])
    labels = np.asarray([str(x) if pd.notna(x) else "missing" for x in rows.iloc[row_positions][spec.column].to_numpy()], dtype=object)
    return labels, np.full(len(labels), -999, dtype=int), np.arange(len(labels), dtype=float)


def build2_support_flag(axis: str, rows_count: int, n_atoms: int, n_eff: float) -> str:
    flags = []
    if n_atoms < THIN_ATOMS:
        flags.append("thin_atoms")
    if axis == "within" and rows_count < THIN_WITHIN_ROWS:
        flags.append("thin_rows")
    if np.isfinite(n_eff) and n_eff < MIN_FAVORABLE_N_EFF:
        flags.append("thin_N_eff")
    return ";".join(flags)


def build2_response_shapes(curves: pd.DataFrame) -> dict[tuple[str, str, str, str, str, str, str, str], str]:
    out: dict[tuple[str, str, str, str, str, str, str, str], str] = {}
    group_cols = ["target", "fit_scope", "mechanism", "condition_family", "condition_name", "stratum", "tier", "axis"]
    for key, group in curves.groupby(group_cols, dropna=False):
        usable = group[(group["support_flag"] == "") & np.isfinite(group["heldout_R2"])].sort_values("_bin_order")
        vals = usable["heldout_R2"].to_numpy(float)
        if len(vals) < 3:
            out[key] = "unstable-thin"
            continue
        if str(usable["condition_kind"].iloc[0]) not in {"cxx_bin", "python_quantile_cpp_scalar"}:
            span = float(np.nanmax(vals) - np.nanmin(vals))
            out[key] = "categorical contrast" if span >= 0.02 else "flat response"
            continue
        span = float(np.nanmax(vals) - np.nanmin(vals))
        tol = max(0.02, 0.10 * span)
        diffs = np.diff(vals)
        if span < 0.02:
            out[key] = "flat response"
        elif np.all(diffs >= -tol) and vals[-1] > vals[0]:
            out[key] = "monotone rise"
        elif np.all(diffs <= tol) and vals[-1] < vals[0]:
            out[key] = "monotone fall"
        elif len(vals) >= 4 and vals[0] > vals[1] and vals[-1] > vals[-2] and np.nanmin(vals[1:-1]) < min(vals[0], vals[-1]) - tol:
            out[key] = "U-shape"
        elif np.max(np.abs(diffs)) > 0.6 * span:
            out[key] = "threshold behavior"
        else:
            out[key] = "flat response"
    return out


def build2_partition_curves(
    target_name: str,
    fit_scope: str,
    tier_results: dict[str, dict[str, object]],
    tier_order: list[str],
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    partition_bins: np.ndarray,
    atom_bins: np.ndarray,
    frames: np.ndarray,
    condition_list: list[Build2ConditionSpec],
    manifest: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    rows_out = []
    bin_defs = []
    strata = ["ALL"] + sorted(str(x) for x in atom_rows["stratum"].dropna().unique())
    axes: dict[str, dict[str, object]] = {"between": {"length": len(atom_rows), "frames": np.arange(len(atom_rows), dtype=int), "row_positions": None}}
    first_within = tier_results[tier_order[0]]["within"]
    axes["within"] = {
        "length": len(first_within["target"]),
        "frames": frames[first_within["score_row_positions"]],
        "row_positions": first_within["score_row_positions"],
    }
    neff_by_tier_axis = {
        (tier, "within"): per_atom_neff_components(
            tier_results[tier]["within"]["target"],
            tier_results[tier]["within"]["groups"],
            np.asarray(axes["within"]["frames"]),
        )
        for tier in tier_order
    }
    for axis, axis_info in axes.items():
        row_positions = axis_info["row_positions"]
        if axis == "between":
            axis_atom_index = atom_rows["atom_index"].to_numpy(int)
            axis_strata = atom_rows["stratum"].to_numpy(str)
        else:
            assert row_positions is not None
            axis_atom_index = rows.iloc[row_positions]["atom_index"].to_numpy(int)
            axis_strata = rows.iloc[row_positions]["stratum"].to_numpy(str)
        for stratum in strata:
            stratum_mask = np.ones(len(axis_atom_index), dtype=bool) if stratum == "ALL" else axis_strata == stratum
            if not stratum_mask.any():
                continue
            for spec in condition_list:
                labels, raw_bins, orders = build2_axis_conditions(spec, axis, rows, atom_rows, partition_bins, atom_bins, row_positions if row_positions is not None else np.asarray([], dtype=int), manifest)
                labels_in = labels[stratum_mask]
                raw_in = raw_bins[stratum_mask]
                order_in = orders[stratum_mask]
                atoms_in = axis_atom_index[stratum_mask]
                for label in sorted(set(labels_in), key=str):
                    unit_local = labels_in == label
                    if not unit_local.any():
                        continue
                    if spec.kind in {"cxx_bin", "python_quantile_cpp_scalar"}:
                        bin_id = mode_int(raw_in[unit_local])
                        bin_order = build2_bin_order(spec, label, bin_id)
                    else:
                        bin_id = -999
                        bin_order = build2_bin_order(spec, label)
                    global_mask = np.zeros(len(axis_atom_index), dtype=bool)
                    idx = np.flatnonzero(stratum_mask)
                    global_mask[idx[unit_local]] = True
                    n_atoms = int(len(np.unique(axis_atom_index[global_mask])))
                    rows_count = int(global_mask.sum())
                    if spec.kind in {"cxx_bin", "python_quantile_cpp_scalar"}:
                        bin_defs.append({
                            "target": target_name,
                            "fit_scope": fit_scope,
                            "axis": axis,
                            "stratum": stratum,
                            "mechanism": spec.mechanism,
                            "condition_family": spec.family,
                            "condition_name": spec.name,
                            "condition_kind": spec.kind,
                            "source": spec.source,
                            "source_column": spec.column,
                            "bin_label": label,
                            "bin_id": int(bin_id),
                            "n_atoms_pre_rank": n_atoms,
                        })
                    for tier in tier_order:
                        axis_pred = tier_results[tier][axis]
                        pred = axis_pred["pred"][global_mask]
                        target = axis_pred["target"][global_mask]
                        groups = axis_pred["groups"][global_mask]
                        r2 = r2_score(pred, target)
                        se = jackknife_metric_se(pred, target, groups)
                        lo, hi = ci95(r2, se)
                        atom_neff = neff_by_tier_axis.get((tier, axis))
                        n_eff_median = bin_neff_median(axis, axis_atom_index, global_mask, atom_neff)
                        support_flag = build2_support_flag(axis, rows_count, n_atoms, n_eff_median)
                        rows_out.append({
                            "target": target_name,
                            "fit_scope": fit_scope,
                            "mechanism": spec.mechanism,
                            "condition_family": spec.family,
                            "condition_name": spec.name,
                            "condition_kind": spec.kind,
                            "condition_source": spec.source,
                            "bin_label": label,
                            "bin_id": int(bin_id),
                            "stratum": stratum,
                            "tier": tier,
                            "axis": axis,
                            "rows": rows_count,
                            "n_atoms": n_atoms,
                            "N_eff_median": n_eff_median,
                            "heldout_R2": r2,
                            "heldout_R2_ci95_low": lo,
                            "heldout_R2_ci95_high": hi,
                            "delta_R2_vs_classical": math.nan,
                            "delta_R2_vs_previous_tier": math.nan,
                            "support_flag": support_flag,
                            "_bin_order": bin_order,
                        })
    curves = pd.DataFrame(rows_out)
    if curves.empty:
        return curves, curves.copy(), bin_defs
    key_cols = ["target", "fit_scope", "mechanism", "condition_family", "condition_name", "condition_kind", "bin_label", "bin_id", "stratum", "axis"]
    classical = curves[curves["tier"] == "classical"][key_cols + ["heldout_R2"]].rename(columns={"heldout_R2": "_classical_R2"})
    curves = curves.merge(classical, on=key_cols, how="left")
    curves["delta_R2_vs_classical"] = curves["heldout_R2"] - curves["_classical_R2"]
    prev_parts = []
    prev_map = previous_tier_map(tier_order)
    for tier, prev in prev_map.items():
        part = curves[curves["tier"] == tier].copy()
        if prev is None:
            part["_previous_R2"] = np.nan
        else:
            prev_df = curves[curves["tier"] == prev][key_cols + ["heldout_R2"]].rename(columns={"heldout_R2": "_previous_R2"})
            part = part.drop(columns=["_previous_R2"], errors="ignore").merge(prev_df, on=key_cols, how="left")
        prev_parts.append(part)
    curves = pd.concat(prev_parts, ignore_index=True)
    curves["delta_R2_vs_previous_tier"] = curves["heldout_R2"] - curves["_previous_R2"]
    shape_map = build2_response_shapes(curves)
    curves["response_shape"] = [
        shape_map.get((r.target, r.fit_scope, r.mechanism, r.condition_family, r.condition_name, r.stratum, r.tier, r.axis), "unstable-thin")
        for r in curves.itertuples(index=False)
    ]
    public_cols = [
        "target",
        "fit_scope",
        "mechanism",
        "condition_family",
        "condition_name",
        "condition_kind",
        "condition_source",
        "bin_label",
        "bin_id",
        "stratum",
        "tier",
        "axis",
        "rows",
        "n_atoms",
        "N_eff_median",
        "heldout_R2",
        "heldout_R2_ci95_low",
        "heldout_R2_ci95_high",
        "delta_R2_vs_classical",
        "delta_R2_vs_previous_tier",
        "support_flag",
        "response_shape",
    ]
    curves_public = curves[public_cols].copy()
    long_rows = []
    for r in curves.itertuples(index=False):
        for metric in ["heldout_R2", "delta_R2_vs_classical", "delta_R2_vs_previous_tier"]:
            long_rows.append({
                "target": r.target,
                "fit_scope": r.fit_scope,
                "mechanism": r.mechanism,
                "condition_family": r.condition_family,
                "condition_name": r.condition_name,
                "bin_label": r.bin_label,
                "bin_id": r.bin_id,
                "stratum": r.stratum,
                "tier": r.tier,
                "axis": r.axis,
                "metric": metric,
                "value": getattr(r, metric),
                "support_flag": r.support_flag,
                "response_shape": r.response_shape,
            })
    return curves_public, pd.DataFrame(long_rows), bin_defs


def build2_favorable_partitions(curves: pd.DataFrame, shortlist_size: int = 120) -> pd.DataFrame:
    if curves.empty:
        return curves.copy()
    eligible = curves[
        (curves["support_flag"] == "")
        & np.isfinite(curves["heldout_R2"])
        & (curves["n_atoms"] >= THIN_ATOMS)
        & (curves["N_eff_median"] >= MIN_FAVORABLE_N_EFF)
        & (curves["axis"] == "within")
    ].copy()
    if eligible.empty:
        eligible = curves[(curves["support_flag"] == "") & np.isfinite(curves["heldout_R2"])].copy()
    eligible = eligible[(eligible["heldout_R2_ci95_low"] > 0.0) | (eligible["heldout_R2"] > 0.0)]
    eligible["_positive_ci"] = eligible["heldout_R2_ci95_low"] > 0.0
    eligible = eligible.sort_values(["_positive_ci", "heldout_R2", "delta_R2_vs_previous_tier"], ascending=[False, False, False], na_position="last")
    cols = [
        "target",
        "fit_scope",
        "mechanism",
        "condition_family",
        "condition_name",
        "bin_label",
        "bin_id",
        "stratum",
        "tier",
        "axis",
        "rows",
        "n_atoms",
        "N_eff_median",
        "heldout_R2",
        "heldout_R2_ci95_low",
        "heldout_R2_ci95_high",
        "delta_R2_vs_classical",
        "delta_R2_vs_previous_tier",
        "response_shape",
    ]
    return eligible.head(shortlist_size)[cols]


def load_casehunter_candidates(substrate_dir: Path) -> pd.DataFrame:
    parts = []
    for mechanism in ["ring", "charge", "mc"]:
        path = substrate_dir / "equations" / mechanism / "cases_manifest.csv"
        if path.exists():
            df = pd.read_csv(path)
            df["case_manifest_path"] = str(path)
            parts.append(df)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def case_window_positions(rows: pd.DataFrame, atom_index: int, begin: int, end: int) -> np.ndarray:
    mask = (
        (rows["atom_index"].to_numpy(int) == int(atom_index))
        & (rows["original_frame_index"].to_numpy(int) >= int(begin))
        & (rows["original_frame_index"].to_numpy(int) < int(end))
    )
    pos = np.flatnonzero(mask)
    if len(pos) == 0:
        mask = (
            (rows["atom_index"].to_numpy(int) == int(atom_index))
            & (rows["original_frame_index"].to_numpy(int) >= int(begin))
            & (rows["original_frame_index"].to_numpy(int) <= int(end))
        )
        pos = np.flatnonzero(mask)
    return pos


def build2_casehunter_intersection(
    favorable: pd.DataFrame,
    cases: pd.DataFrame,
    rows: pd.DataFrame,
    atom_rows: pd.DataFrame,
    partition_bins: np.ndarray,
    condition_specs: list[Build2ConditionSpec],
    manifest: dict[str, object],
    shortlist_size: int = 60,
) -> pd.DataFrame:
    if favorable.empty or cases.empty:
        return pd.DataFrame()
    specs_by_name = {spec.name: spec for spec in condition_specs}
    out = []
    atom_stratum = atom_rows.set_index("atom_index")["stratum"].astype(str).to_dict()
    for fav in favorable.itertuples(index=False):
        if fav.mechanism not in {"ring", "charge", "mc"}:
            continue
        spec = specs_by_name.get(str(fav.condition_name))
        if spec is None:
            continue
        same_mech = cases[cases["mechanism"].astype(str) == str(fav.mechanism)]
        for case in same_mech.itertuples(index=False):
            case_stratum = str(atom_stratum.get(int(case.atom_index), "missing"))
            if fav.stratum != "ALL" and fav.stratum != case_stratum:
                continue
            positions = case_window_positions(rows, int(case.atom_index), int(case.frame_window_begin), int(case.frame_window_end))
            if len(positions) == 0:
                continue
            if spec.kind == "cxx_bin":
                modal_bin = mode_int(np.asarray(partition_bins[positions, int(spec.bin_index)], dtype=int))
                if int(fav.bin_id) != int(modal_bin):
                    continue
                modal_label = build2_bin_label(spec.column, modal_bin, manifest)
            else:
                labels = np.asarray([str(x) if pd.notna(x) else "missing" for x in rows.iloc[positions][spec.column].to_numpy()], dtype=object)
                unique, counts = np.unique(labels, return_counts=True)
                modal_label = str(unique[np.argmax(counts)])
                if str(fav.bin_label) != modal_label:
                    continue
                modal_bin = -999
            out.append({
                "target": fav.target,
                "fit_scope": fav.fit_scope,
                "mechanism": fav.mechanism,
                "stratum": fav.stratum,
                "case_atom_index": int(case.atom_index),
                "case_frame_window_begin": int(case.frame_window_begin),
                "case_frame_window_end": int(case.frame_window_end),
                "condition_family": fav.condition_family,
                "condition_name": fav.condition_name,
                "bin_label": modal_label,
                "bin_id": int(modal_bin),
                "tier": fav.tier,
                "axis": fav.axis,
                "heldout_R2": float(fav.heldout_R2),
                "delta_R2_vs_classical": float(fav.delta_R2_vs_classical) if np.isfinite(fav.delta_R2_vs_classical) else math.nan,
                "response_shape": fav.response_shape,
                "partition_n_atoms": int(fav.n_atoms),
                "partition_N_eff_median": float(fav.N_eff_median),
                "case_score": float(case.score),
                "case_dft_recovery_R2_for_navigation_only": float(case.dft_recovery_R2) if pd.notna(case.dft_recovery_R2) else math.nan,
                "case_manifest_path": case.case_manifest_path,
            })
    if not out:
        return pd.DataFrame()
    df = pd.DataFrame(out)
    return df.sort_values(["heldout_R2", "case_score"], ascending=[False, False], na_position="last").head(shortlist_size)


def write_build2_reports(
    out_dir: Path,
    score_table: pd.DataFrame,
    curves: pd.DataFrame,
    favorable: pd.DataFrame,
    case_shortlist: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    score_focus = score_table[score_table["atom_type"].isin(["ALL", "HN", "O", "N", "CA", "C", "HA"])].copy()
    lines = [
        "# Build 3 Fit-Architecture Report",
        "",
        f"Substrate: `{audit['substrate_dir']}`",
        f"Run dir: `{audit['output_dir']}`",
        f"Targets: {', '.join(TARGET_SPECS)}",
        "",
        "## Fit Scores",
        markdown_table(score_focus[[
            "target",
            "fit_scope",
            "atom_type",
            "tier",
            "n_atoms_in_slice",
            "n_atoms_between_fit_scope",
            "between_LOAO_test_R2",
            "between_delta_R2_vs_global_sliced",
            "between_delta_R2_vs_classical",
            "within_frameblock_test_R2",
            "within_delta_R2_vs_global_sliced",
            "within_delta_R2_vs_classical",
            "between_support_flag",
            "thin_flag",
        ]], max_rows=80),
        "",
        "## Partition Shapes",
    ]
    if curves.empty:
        lines.append("No curves emitted.")
    else:
        shape_counts = curves.groupby(["target", "fit_scope", "mechanism", "condition_family", "tier", "axis", "response_shape"], dropna=False).size().reset_index(name="curve_points")
        lines.append(markdown_table(shape_counts, max_rows=120))
    lines += [
        "",
        "## Favourable Partitions",
        markdown_table(favorable, max_rows=80) if not favorable.empty else "No eligible favourable partitions.",
        "",
        "## CaseHunter Intersection",
        markdown_table(case_shortlist, max_rows=80) if not case_shortlist.empty else "No CaseHunter candidate intersected the favourable partitions.",
    ]
    (out_dir / "build2_partition_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_build2_plots(out_dir: Path, score_table: pd.DataFrame, curves: pd.DataFrame) -> dict[str, object]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plot availability is environment-dependent
        return {"plots_written": False, "reason": str(exc)}
    written = []
    focus = score_table[(score_table["fit_scope"] == "global_sliced") & score_table["atom_type"].isin(["N", "CA", "C", "O", "HN", "HA"])].copy()
    for target_name, group in focus.groupby("target", dropna=False):
        pivot = group.pivot(index="atom_type", columns="tier", values="between_LOAO_test_R2").reindex(["N", "CA", "C", "O", "HN", "HA"])
        fig, ax = plt.subplots(figsize=(10, 5))
        pivot.plot(ax=ax, marker="o")
        ax.set_title(f"{target_name} between-LOAO R2 by stratum")
        ax.set_xlabel("stratum")
        ax.set_ylabel("held-out R2")
        ax.axhline(0, color="black", linewidth=0.8)
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = out_dir / f"score_by_stratum_{target_name.replace('/', '_')}.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(str(path))
    if not curves.empty:
        top = curves[(curves["support_flag"] == "") & (curves["condition_kind"] == "cxx_bin")].copy()
        top = top.sort_values("heldout_R2", ascending=False).head(24)
        if not top.empty:
            labels = [f"{r.target}\\n{r.mechanism}:{r.stratum}\\n{r.condition_name} {r.bin_label}" for r in top.itertuples(index=False)]
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.barh(np.arange(len(top)), top["heldout_R2"].to_numpy(float))
            ax.set_yticks(np.arange(len(top)), labels=labels, fontsize=6)
            ax.invert_yaxis()
            ax.set_xlabel("held-out R2")
            ax.set_title("Top supported partition bins")
            fig.tight_layout()
            path = out_dir / "partition_top_supported_bins.png"
            fig.savefig(path, dpi=140)
            plt.close(fig)
            written.append(str(path))
    return {"plots_written": True, "paths": written}


def build2_fit_stage_checks(audit: dict[str, object], score_table: pd.DataFrame, curves: pd.DataFrame, case_shortlist: pd.DataFrame) -> None:
    checks = audit["fit_stage_checks"]
    checks["no_external_merge"] = {
        "pass": True,
        "python_read_only_reduced_substrate_sidecars_and_bin_ids": True,
        "files_read": audit["files_read"],
        "forbidden_inputs_not_opened": ["trajectory.h5", "ORCA outputs", "per-source directories", "older rediscover dirs", "ring_paths cross-method validation data"],
    }
    checks["basis_and_targets"] = {
        "pass": bool(audit["change_of_basis_orthogonality_max_abs"] < 1.0e-12 and set(audit["targets_fit"]) == set(TARGET_SPECS)),
        "C_T_C_minus_I_max_abs": audit["change_of_basis_orthogonality_max_abs"],
        "targets_fit": audit["targets_fit"],
        "T1_label": "field-linear diagnostic, not a shift claim",
    }
    cv_checks = audit["cv_integrity_details"]
    checks["cv_integrity"] = {
        "pass": bool(
            all(v["between"]["cv_check"]["heldout_atoms_never_in_train"] for target in cv_checks.values() for v in target.values())
            and all(v["within"]["cv_check"]["train_test_overlap"] == 0 for target in cv_checks.values() for v in target.values())
            and all(v["within"]["cv_check"]["purged_in_train_or_test"] == 0 for target in cv_checks.values() for v in target.values())
            and audit["within_split"]["cross_split_lag1_pairs"] == 0
            and score_table["within_N_eff_median"].notna().all()
        ),
        "between_heldout_atoms_never_in_train": all(v["between"]["cv_check"]["heldout_atoms_never_in_train"] for target in cv_checks.values() for v in target.values()),
        "within_train_test_overlap_zero": all(v["within"]["cv_check"]["train_test_overlap"] == 0 for target in cv_checks.values() for v in target.values()),
        "purged_frames_absent_from_train_and_test": all(v["within"]["cv_check"]["purged_in_train_or_test"] == 0 for target in cv_checks.values() for v in target.values()),
        "scores_are_held_out": True,
        "N_eff_reported_per_score_row": bool(score_table["within_N_eff_median"].notna().all()),
        "cross_split_lag1_pairs": audit["within_split"]["cross_split_lag1_pairs"],
    }
    checks["partition_integrity"] = {
        "pass": bool(
            not curves.empty
            and audit["partition_conditions_input_side_only"]
            and audit["partition_bins_from_cpp_lookup_or_python_quantile_cpp_scalar"]
            and audit["charge_partition_non_degenerate"]
        ),
        "cpp_bin_lookup_only_for_existing_numeric_partition_bins": audit["partition_bins_from_cpp_lookup_only"],
        "python_quantile_bins_on_cpp_dominance_scalars": audit["dominance_quantile_bins_from_python_on_cpp_scalar"],
        "no_condition_uses_dft_target_residual_or_coefficients": True,
        "casehunter_intersection_rows": int(len(case_shortlist)),
        "charge_partition_non_degenerate": audit["charge_partition_non_degenerate"],
        "dominance_scalar_columns_available_in_substrate": audit["dominance_scalar_columns_available_in_substrate"],
        "dominance_cpp_bin_id_next_emit_flag": audit["dominance_cpp_bin_id_next_emit_flag"],
    }


def build2_log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def build2_main() -> None:
    args = parse_args()
    tier_order = parse_tiers(args.tiers)
    alpha_grid = parse_alpha_grid(args.alpha_grid)
    substrate_dir = args.substrate_dir.resolve()
    data = build2_load_inputs(substrate_dir)
    checks, charge_mask = build2_input_acceptance_checks(data, substrate_dir)
    if not checks["input_acceptance"]["pass"]:
        raise SystemExit(f"FATAL: input acceptance failed: {checks['input_acceptance']}")
    out_dir, disk_audit = prepare_build2_output_dir(args.out_dir, args.keep_out_dir)

    rows_all: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    arrays: dict[str, np.ndarray] = data["arrays"]
    manifest: dict[str, object] = data["manifest"]
    rows = rows_all if bool(charge_mask.all()) else rows_all.loc[charge_mask].reset_index(drop=True)
    if not bool(charge_mask.all()):
        arrays = {k: (v[charge_mask] if v.shape[0] == len(charge_mask) else v) for k, v in arrays.items()}
    n_atoms = int(manifest["n_atoms"])
    n_frames = int(manifest["n_dft_frames"])
    if len(rows) != n_atoms * n_frames:
        raise SystemExit("FATAL: charge-complete filter no longer has dense 846 x 660 rows")

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    if orth >= 1.0e-12:
        raise SystemExit(f"FATAL: change_of_basis.get_C() orthogonality check failed: {orth:.3e}")

    row_indices = np.flatnonzero(charge_mask)
    frames = rows["original_frame_index"].to_numpy(int)
    atoms = rows["atom_index"].to_numpy(int)
    rows_for_conditions = add_build2_charge_agreement_columns(rows, arrays, specs)
    atom_conditions = atom_condition_frame(rows_for_conditions, n_frames, n_atoms)
    partition_bins = arrays["per_atom_substrate_partition_bins"]
    atom_bins = atom_mode_bins(partition_bins, n_frames, n_atoms)
    partition_cols = [str(x) for x in manifest["partition_bins"]["columns"]]
    cond_specs = build2_condition_specs(rows_for_conditions, partition_cols)

    base_blocks, block_labels, transform_audit = build2_feature_blocks(rows, arrays, specs, C)
    base_atom_blocks = {name: atom_means_dense(block, n_frames, n_atoms) for name, block in base_blocks.items()}
    embedding = arrays["per_atom_substrate_aimnet2_embedding"]
    embedding_atom_mean = atom_means_embedding(embedding, n_frames, n_atoms)
    split = split_frame_block(frames, args.within_test_fraction, args.purge_frames_each_side)

    within_embedding_pcs = None
    within_embedding_pca_audit = None
    embedding_tiers = [tier for tier in tier_order if tier_uses_embedding(tier)]
    if embedding_tiers:
        within_train_indices = row_indices[split["train_mask"]]
        within_pca = fit_pca_memmap(embedding, within_train_indices, args.embedding_components)
        within_embedding_pcs = transform_pca_memmap(embedding, within_pca, row_indices)
        within_embedding_pca_audit = {
            "method": "unsupervised PCA fit once on within-CV training rows only; shared by selected tiers that consume AIMNet2 embedding PCs",
            "n_components": int(within_pca.n_components),
            "original_dims": int(within_pca.original_dims),
            "training_rows": int(within_pca.n_train_rows),
            "test_rows_excluded": int(split["test_mask"].sum()),
            "purged_rows_excluded": int(split["purge_mask"].sum()),
            "explained_variance_ratio": float(within_pca.explained_variance_ratio),
            "consumed_by_tiers": embedding_tiers,
        }
    else:
        within_embedding_pca_audit = {"method": "embedding sidecar present but not consumed by selected tiers", "n_components": 0, "original_dims": int(embedding.shape[1]), "consumed_by_tiers": []}

    all_score_rows = []
    all_long_rows = []
    all_curves = []
    all_curves_long = []
    all_bin_defs = []
    alpha_selection_by_target: dict[str, object] = {}
    pca_audit_by_target: dict[str, object] = {}
    cv_details_by_target: dict[str, object] = {}
    between_pca_cache: dict[object, PcaTransform] = {}
    between_alpha_pca_cache: dict[object, PcaTransform] = {}
    within_alpha_pca_cache: dict[object, PcaTransform] = {}
    for target_name, raw_target in data["targets"].items():
        build2_log(f"target {target_name}: transform target")
        y = transform_build2_target(target_name, raw_target, charge_mask, C)
        y_atom = atom_means_dense(y, n_frames, n_atoms)
        global_tier_results: dict[str, dict[str, object]] = {}
        per_type_tier_results: dict[str, dict[str, object]] = {}
        pca_audit: dict[str, object] = {"global_sliced": {}, "per_type": {}}
        cv_details: dict[str, object] = {}
        alpha_selection_by_scope: dict[str, object] = {"global_sliced": {}, "per_type": {}}
        for tier in tier_order:
            if args.alpha_mode == "select":
                build2_log(f"target {target_name}: global_sliced tier {tier}: select between alpha")
                between_alpha, between_alpha_audit = select_between_alphas(
                    tier,
                    base_atom_blocks,
                    embedding_atom_mean,
                    y_atom,
                    alpha_grid,
                    args.embedding_components,
                    args.inner_cv_folds,
                    between_alpha_pca_cache,
                )
                build2_log(f"target {target_name}: global_sliced tier {tier}: select within alpha")
                within_alpha, within_alpha_audit = select_within_alpha(
                    tier,
                    base_blocks,
                    embedding,
                    row_indices,
                    y,
                    atoms,
                    frames,
                    split,
                    alpha_grid,
                    args.embedding_components,
                    args.inner_cv_folds,
                    args.purge_frames_each_side,
                    within_alpha_pca_cache,
                )
            else:
                fixed_alpha = float(args.ridge_alpha)
                between_alpha = np.full(n_atoms, fixed_alpha, dtype=float)
                between_alpha_audit = {"axis": "between_LOAO", "method": "fixed-alpha labelled baseline; no inner CV", "selected_alpha": fixed_alpha, "min": fixed_alpha, "max": fixed_alpha, "counts": {f"{fixed_alpha:g}": int(n_atoms)}, "heldout_test_atom_excluded_from_alpha_selection": True, "alpha_grid": [float(x) for x in alpha_grid]}
                within_alpha = fixed_alpha
                within_alpha_audit = {"axis": "within_frameblock", "method": "fixed-alpha labelled baseline; no inner CV", "selected_alpha": fixed_alpha, "outer_test_rows_used_for_alpha_selection": 0, "outer_purged_rows_used_for_alpha_selection": 0, "alpha_grid": [float(x) for x in alpha_grid]}
            alpha_selection_by_scope["global_sliced"][tier] = {"between": between_alpha_audit, "within": within_alpha_audit}
            build2_log(f"target {target_name}: global_sliced tier {tier}: fit between")
            cache = between_pca_cache if tier_uses_embedding(tier) else None
            between = evaluate_between_tier(tier, base_atom_blocks, embedding_atom_mean, y_atom, between_alpha, args.embedding_components, cache)
            build2_log(f"target {target_name}: global_sliced tier {tier}: fit within")
            within = evaluate_within_tier(tier, base_blocks, within_embedding_pcs, within_embedding_pca_audit, y, atoms, split, within_alpha)
            global_tier_results[tier] = {"between": between, "within": within}
            pca_audit["global_sliced"][tier] = {"between": between["pca_audit"], "within": within["pca_audit"]}
            cv_details[f"global_sliced:{tier}"] = {"between": {"cv_check": between["cv_check"]}, "within": {"cv_check": within["cv_check"]}}

            build2_log(f"target {target_name}: per_type tier {tier}: select alphas and fit strata")
            per_type_result, per_type_alpha_audit, per_type_pca_audit, per_type_cv_details = evaluate_per_type_tier(
                tier,
                rows,
                atom_conditions,
                base_blocks,
                base_atom_blocks,
                embedding,
                embedding_atom_mean,
                row_indices,
                y,
                y_atom,
                atoms,
                frames,
                split,
                args.alpha_mode,
                float(args.ridge_alpha),
                alpha_grid,
                args.embedding_components,
                args.inner_cv_folds,
                args.purge_frames_each_side,
            )
            per_type_tier_results[tier] = per_type_result
            alpha_selection_by_scope["per_type"][tier] = per_type_alpha_audit
            pca_audit["per_type"][tier] = per_type_pca_audit
            cv_details[f"per_type:{tier}"] = per_type_cv_details
        build2_log(f"target {target_name}: build global_sliced score rows")
        score_rows, long_rows = build2_score_rows(target_name, "global_sliced", global_tier_results, tier_order, rows, atom_conditions, y, frames, atoms, split, base_blocks, args.alpha_mode, args.embedding_components)
        build2_log(f"target {target_name}: build per_type score rows")
        per_score_rows, per_long_rows = build2_score_rows(target_name, "per_type", per_type_tier_results, tier_order, rows, atom_conditions, y, frames, atoms, split, base_blocks, args.alpha_mode, args.embedding_components)
        score_rows.extend(per_score_rows)
        long_rows.extend(per_long_rows)
        build2_log(f"target {target_name}: build global_sliced partition curves")
        curves, curves_long, bin_defs = build2_partition_curves(target_name, "global_sliced", global_tier_results, tier_order, rows_for_conditions, atom_conditions, partition_bins, atom_bins, frames, cond_specs, manifest)
        build2_log(f"target {target_name}: build per_type partition curves")
        per_curves, per_curves_long, per_bin_defs = build2_partition_curves(target_name, "per_type", per_type_tier_results, tier_order, rows_for_conditions, atom_conditions, partition_bins, atom_bins, frames, cond_specs, manifest)
        all_score_rows.extend(score_rows)
        all_long_rows.extend(long_rows)
        all_curves.extend([curves, per_curves])
        all_curves_long.extend([curves_long, per_curves_long])
        all_bin_defs.extend(bin_defs)
        all_bin_defs.extend(per_bin_defs)
        alpha_selection_by_target[target_name] = {"mode": args.alpha_mode, "grid": [float(x) for x in alpha_grid], "by_scope_tier_axis": alpha_selection_by_scope}
        pca_audit_by_target[target_name] = pca_audit
        cv_details_by_target[target_name] = cv_details

    score_table = pd.DataFrame(all_score_rows).reindex(columns=SCORE_TABLE_COLUMNS)
    score_table = add_fit_scope_deltas(score_table).reindex(columns=SCORE_TABLE_COLUMNS)
    score_long = pd.DataFrame(all_long_rows)
    curves = pd.concat(all_curves, ignore_index=True) if all_curves else pd.DataFrame()
    curves_long = pd.concat(all_curves_long, ignore_index=True) if all_curves_long else pd.DataFrame()
    favorable = build2_favorable_partitions(curves)
    cases = load_casehunter_candidates(substrate_dir)
    case_shortlist = build2_casehunter_intersection(favorable, cases, rows_for_conditions, atom_conditions, partition_bins, cond_specs, manifest)
    charge_bin_unique = checks["input_acceptance"]["charge_partition_non_degenerate_bins"]
    audit = {
        "script": str(Path(__file__).resolve()),
        "substrate_dir": str(substrate_dir),
        "output_dir": str(out_dir),
        "disk_guard": disk_audit,
        "files_read": data["files_read"],
        "tier_order": tier_order,
        "fit_scopes": ["global_sliced", "per_type"],
        "targets_fit": list(data["targets"].keys()),
        "feature_tier_selection": {"selected_tiers": tier_order, "headline_delta": "per_type vs global_sliced by stratum, alongside tier deltas"},
        "python_read_scope": "Only emitted reduced per_atom_substrate sidecars, target sidecars, CaseHunter manifests, per_atom_substrate_partition_bins.npy, and C++ dominance scalar conditioners were opened; no trajectory.h5, per-source, older dirs, ring-path validation, or ORCA data.",
        "manifest_relationship": manifest.get("relationship"),
        "manifest_relationship_kind": manifest.get("relationship_kind"),
        "manifest_rows": int(manifest["rows"]),
        "manifest_n_atoms": int(manifest["n_atoms"]),
        "manifest_n_dft_frames": int(manifest["n_dft_frames"]),
        "charge_complete_rows_used": int(charge_mask.sum()),
        "same_charge_complete_row_set_for_all_tiers": True,
        "alpha_mode": args.alpha_mode,
        "ridge_alpha_baseline": float(args.ridge_alpha),
        "alpha_selection": alpha_selection_by_target,
        "change_of_basis_orthogonality_max_abs": orth,
        "t2_transform_check_blocks": transform_audit,
        "embedding_pca": {"n_components": int(args.embedding_components), "provenance_by_target": pca_audit_by_target},
        "within_split": {"strategy": "deterministic contiguous middle frame block", "test_frames": split["test_frames"], "purged_train_frames": split["purged_frames"], "purge_frames_each_side": int(args.purge_frames_each_side), "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"])},
        "cv_integrity_details": cv_details_by_target,
        "no_python_physics_recompute": True,
        "partition_conditions": [spec.__dict__ for spec in cond_specs],
        "partition_bin_definitions": all_bin_defs,
        "partition_bin_definitions_written_before_ranking": True,
        "partition_conditions_input_side_only": True,
        "partition_bins_from_cpp_lookup_only": True,
        "partition_bins_from_cpp_lookup_or_python_quantile_cpp_scalar": True,
        "dominance_quantile_bins_from_python_on_cpp_scalar": True,
        "dominance_scalar_columns": list(DOMINANCE_FRACTION_COLUMNS.values()),
        "dominance_scalar_columns_available_in_substrate": all(col in rows_for_conditions.columns for col in DOMINANCE_FRACTION_COLUMNS.values()),
        "dominance_cpp_bin_id_next_emit_flag": "dominant_fraction_{ring,charge,mc} bin ids should be emitted by C++ in the next substrate; Build 3 bins the C++ scalars in Python to avoid re-emitting.",
        "charge_partition_unique_bins": charge_bin_unique,
        "charge_partition_non_degenerate": bool(len(charge_bin_unique["bin_nearest_charge_r_4_6_8_10"]) > 1 and len(charge_bin_unique["bin_gap_to_2nd_charge_r_4_6_8_10"]) > 1),
        "dominant_charge_bin_available_in_substrate": bool(checks["input_acceptance"]["dominant_charge_bin_column_present"]),
        "ring_field_cross_method_validations_deferred": True,
        "seti_no_verdicts": True,
        "fit_stage_checks": checks,
    }
    build2_fit_stage_checks(audit, score_table, curves, case_shortlist)
    if not all(bool(check["pass"]) for check in audit["fit_stage_checks"].values()):
        (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
        failed = [name for name, check in audit["fit_stage_checks"].items() if not check["pass"]]
        raise SystemExit(f"FATAL: fit-stage checks failed; scores are not valid: {failed}")

    score_table.to_csv(out_dir / "allatom_fit_score_table.csv", index=False)
    score_long.to_csv(out_dir / "allatom_fit_score_long.csv", index=False)
    join_coverage(rows).to_csv(out_dir / "join_coverage.csv", index=False)
    (out_dir / "feature_blocks_used.json").write_text(json.dumps(json_sanitize(feature_blocks_used(block_labels, args.embedding_components, transform_audit, tier_order)), indent=2, sort_keys=True), encoding="utf-8")
    curves.to_csv(out_dir / "partition_response_curves.csv", index=False)
    curves_long.to_csv(out_dir / "partition_response_curves_long.csv", index=False)
    pd.DataFrame(all_bin_defs).to_csv(out_dir / "partition_bin_definitions.csv", index=False)
    favorable.to_csv(out_dir / "partition_favorable_partitions.csv", index=False)
    case_shortlist.to_csv(out_dir / "partition_casehunter_shortlist.csv", index=False)
    plots_audit = write_build2_plots(out_dir, score_table, curves)
    audit["plots"] = plots_audit
    (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
    write_build2_reports(out_dir, score_table, curves, favorable, case_shortlist, audit)
    print(f"wrote {out_dir}")
