# Legacy Piece-2 entrypoint retained for import compatibility.
from allatom_fit_common import *

def piece2_main() -> None:
    args = parse_args()
    tier_order = parse_tiers(args.tiers)
    alpha_grid = parse_alpha_grid(args.alpha_grid)
    substrate_dir = args.substrate_dir.resolve()
    out_dir = args.out_dir
    data = load_inputs(substrate_dir)
    checks, charge_mask = input_acceptance_checks(data, substrate_dir)
    rows_all: pd.DataFrame = data["rows"]
    specs: dict[str, object] = data["specs"]
    classical_all: np.ndarray = data["classical"]
    conditioning_all: np.ndarray = data["conditioning"]
    modulation: np.ndarray = data["modulation"]
    embedding: np.ndarray = data["embedding"]
    target_T2: np.ndarray = data["target_T2"]
    manifest: dict[str, object] = data["manifest"]

    row_indices = np.flatnonzero(charge_mask)
    rows = rows_all.loc[charge_mask].reset_index(drop=True)
    classical = classical_all[charge_mask]
    conditioning = conditioning_all[charge_mask]
    n_atoms = int(manifest["n_atoms"])
    n_frames = int(manifest["n_dft_frames"])
    if len(rows) != n_atoms * n_frames:
        raise SystemExit("FATAL: charge-complete filter no longer has dense 846 x 660 rows")

    C = cob.get_C()
    orth = float(np.abs(C.T @ C - np.eye(5)).max())
    if orth >= 1.0e-12:
        raise SystemExit(f"FATAL: change_of_basis.get_C() orthogonality check failed: {orth:.3e}")
    y = target_T2[charge_mask, :5] @ C.T
    frames = rows["original_frame_index"].to_numpy(int)
    atoms = rows["atom_index"].to_numpy(int)

    base_blocks, block_labels, transform_audit = build_non_embedding_blocks(rows, classical, specs, C)
    row_conditions, formulas = conditioning_frame(rows, conditioning, modulation, classical, specs)
    atom_conditions = atom_condition_frame(row_conditions, n_frames, n_atoms)
    cond_specs = condition_specs(row_conditions)
    split = split_frame_block(frames, args.within_test_fraction, args.purge_frames_each_side)

    base_atom_blocks = {name: atom_means_dense(block, n_frames, n_atoms) for name, block in base_blocks.items()}
    y_atom = atom_means_dense(y, n_frames, n_atoms)
    embedding_atom_mean = atom_means_embedding(embedding, n_frames, n_atoms)
    print("loaded substrate and built feature blocks", flush=True)

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
        print("fit within-split training-only embedding PCA", flush=True)
    else:
        within_embedding_pca_audit = {
            "method": "embedding sidecar present in substrate but not consumed by the selected fit tiers",
            "n_components": 0,
            "original_dims": int(embedding.shape[1]),
            "consumed_by_tiers": [],
        }

    tier_results: dict[str, dict[str, object]] = {}
    pca_audit: dict[str, object] = {}
    cv_details: dict[str, object] = {}
    alpha_selection_by_tier_axis: dict[str, dict[str, object]] = {}
    between_pca_cache: dict[object, PcaTransform] = {}
    between_alpha_pca_cache: dict[object, PcaTransform] = {}
    within_alpha_pca_cache: dict[object, PcaTransform] = {}
    for tier in tier_order:
        if args.alpha_mode == "select":
            print(f"selecting alpha for tier {tier}: between LOAO inner CV", flush=True)
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
            print(f"selecting alpha for tier {tier}: within frame-block inner CV", flush=True)
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
            between_alpha_audit = {
                "axis": "between_LOAO",
                "method": "fixed-alpha labelled baseline; no inner CV",
                "selected_alpha": fixed_alpha,
                "min": fixed_alpha,
                "max": fixed_alpha,
                "counts": {f"{fixed_alpha:g}": int(n_atoms)},
                "heldout_test_atom_excluded_from_alpha_selection": True,
                "alpha_grid": [float(x) for x in alpha_grid],
            }
            within_alpha = fixed_alpha
            within_alpha_audit = {
                "axis": "within_frameblock",
                "method": "fixed-alpha labelled baseline; no inner CV",
                "selected_alpha": fixed_alpha,
                "outer_test_rows_used_for_alpha_selection": 0,
                "outer_purged_rows_used_for_alpha_selection": 0,
                "alpha_grid": [float(x) for x in alpha_grid],
            }
        alpha_selection_by_tier_axis[tier] = {"between": between_alpha_audit, "within": within_alpha_audit}

        print(f"fitting tier {tier}: between LOAO", flush=True)
        cache = between_pca_cache if tier_uses_embedding(tier) else None
        between = evaluate_between_tier(tier, base_atom_blocks, embedding_atom_mean, y_atom, between_alpha, args.embedding_components, cache)
        print(f"fitting tier {tier}: within frame block", flush=True)
        within = evaluate_within_tier(
            tier,
            base_blocks,
            within_embedding_pcs,
            within_embedding_pca_audit,
            y,
            atoms,
            split,
            within_alpha,
        )
        tier_results[tier] = {"between": between, "within": within}
        pca_audit[tier] = {"between": between["pca_audit"], "within": within["pca_audit"]}
        cv_details[tier] = {
            "between": {"cv_check": between["cv_check"]},
            "within": {"cv_check": within["cv_check"]},
        }

    print("building score tables", flush=True)
    score_rows, long_rows = build_score_rows(
        tier_results,
        tier_order,
        rows,
        atom_conditions,
        y,
        frames,
        atoms,
        split,
        base_blocks,
        args.alpha_mode,
        args.embedding_components,
    )
    score_table = pd.DataFrame(score_rows)[SCORE_TABLE_COLUMNS]
    score_long = pd.DataFrame(long_rows)
    print("building partition curves", flush=True)
    curves, curves_long, bin_defs = partition_curves(tier_results, tier_order, atom_conditions, row_conditions, frames, cond_specs)
    audit = {
        "script": str(Path(__file__).resolve()),
        "substrate_dir": str(substrate_dir),
        "output_dir": str(out_dir),
        "files_read": data["files_read"],
        "tier_order": tier_order,
        "feature_tier_selection": {
            "selected_tiers": tier_order,
            "scope": "analysis-side fit feature consumption only; emitted substrate keeps AIMNet2 columns present",
        },
        "python_read_scope": "Only files inside the emitted per_atom_substrate directory listed in files_read; no trajectory.h5, ORCA, older broad-backbone, MOPAC, all-atom-equivariant, source, pair, or external merge directories were opened.",
        "manifest_relationship": manifest.get("relationship"),
        "manifest_relationship_kind": manifest.get("relationship_kind"),
        "manifest_rows": int(manifest["rows"]),
        "manifest_n_atoms": int(manifest["n_atoms"]),
        "manifest_n_dft_frames": int(manifest["n_dft_frames"]),
        "charge_complete_rows_used": int(charge_mask.sum()),
        "same_charge_complete_row_set_for_all_tiers": True,
        "alpha_mode": args.alpha_mode,
        "ridge_alpha_baseline": float(args.ridge_alpha),
        "alpha_selection": {
            "mode": args.alpha_mode,
            "grid": [float(x) for x in alpha_grid],
            "fixed_alpha_baseline": float(args.ridge_alpha),
            "inner_cv_folds_requested": int(args.inner_cv_folds),
            "by_tier_axis": alpha_selection_by_tier_axis,
        },
        "target": "DFT total T2, per_atom_substrate_target_T2.npy[:, 0:5] @ change_of_basis.get_C().T",
        "target_T2_original_shape": list(target_T2.shape),
        "change_of_basis_orthogonality_max_abs": orth,
        "t2_transform_check_blocks": [
            {"block": item["feature_block"], "components": 5, "transform": item["transform"]}
            for item in transform_audit
        ],
        "embedding_pca": {
            "owner_override": "training-only PCA is primary; full 256d embedding is not used as primary",
            "n_components": int(args.embedding_components),
            "provenance_by_tier": pca_audit,
        },
        "within_split": {
            "strategy": "deterministic contiguous middle frame block",
            "test_frames": split["test_frames"],
            "purged_train_frames": split["purged_frames"],
            "purge_frames_each_side": int(args.purge_frames_each_side),
            "cross_split_lag1_pairs": int(split["cross_split_lag1_pairs"]),
        },
        "cv_integrity_details": cv_details,
        "no_python_physics_recompute": True,
        "algebraic_partition_conditioners": formulas,
        "partition_conditions": [spec.__dict__ for spec in cond_specs],
        "partition_bin_definitions": bin_defs,
        "partition_bin_definitions_written_before_ranking": True,
        "partition_conditions_input_side_only": True,
        "fit_stage_checks": checks,
    }
    fit_stage_checks(audit, score_table, curves, favorable_cases(curves))
    if not all(bool(check["pass"]) for check in audit["fit_stage_checks"].values()):
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
        failed = [name for name, check in audit["fit_stage_checks"].items() if not check["pass"]]
        raise SystemExit(f"FATAL: fit-stage checks failed; scores are not valid: {failed}")

    fav = favorable_cases(curves)
    if out_dir.exists() and not args.keep_out_dir:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    score_table.to_csv(out_dir / "allatom_fit_score_table.csv", index=False)
    score_long.to_csv(out_dir / "allatom_fit_score_long.csv", index=False)
    join_coverage(rows).to_csv(out_dir / "join_coverage.csv", index=False)
    (out_dir / "feature_blocks_used.json").write_text(
        json.dumps(json_sanitize(feature_blocks_used(block_labels, args.embedding_components, transform_audit, tier_order)), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    curves.to_csv(out_dir / "partition_response_curves.csv", index=False)
    curves_long.to_csv(out_dir / "partition_response_curves_long.csv", index=False)
    fav.to_csv(out_dir / "partition_favorable_cases.csv", index=False)
    (out_dir / "run_audit.json").write_text(json.dumps(json_sanitize(audit), indent=2, sort_keys=True), encoding="utf-8")
    write_reports(out_dir, score_table, curves, fav, audit)
    print(f"wrote {out_dir}")
    print("fit-stage checks: " + ", ".join(f"{name}=pass" for name in audit["fit_stage_checks"]))


