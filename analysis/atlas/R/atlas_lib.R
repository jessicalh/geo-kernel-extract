options(stringsAsFactors = FALSE)

atlas_required_packages <- c(
  "arrow", "data.table", "jsonlite", "RcppCNPy", "mgcv", "lme4",
  "relaimpo", "circular", "ggplot2", "patchwork", "ggrepel", "viridis",
  "ragg", "dplyr"
)

atlas_optional_packages <- c(
  "ggtext", "ggdist", "ggridges", "ggnewscale", "scico", "MetBrewer",
  "sandwich", "clubSandwich", "fst", "qs"
)

atlas_require_packages <- function(required = atlas_required_packages) {
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      "\nInstall with Rscript analysis/atlas/install_atlas_packages.R",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

atlas_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

atlas_default_config <- function() {
  root <- "/shared/2026Thesis/nmr-shielding-emit-build"
  input <- file.path(root, "output", "consolidated_fp_emit_20260612T000000Z")
  output <- Sys.getenv(
    "ATLAS_OUTPUT_DIR",
    file.path(root, "output", "contribution_atlas_20260612T000000Z")
  )
  list(
    repo_root = root,
    input_dir = Sys.getenv("ATLAS_INPUT_DIR", input),
    output_dir = output,
    table_dir = file.path(output, "tables"),
    figure_dir = file.path(output, "figures"),
    checkpoint_dir = file.path(output, "checkpoints"),
    log_dir = file.path(output, "logs"),
    max_block_cols = as.integer(Sys.getenv("ATLAS_MAX_BLOCK_COLS", "16")),
    ridge_lambda = as.numeric(Sys.getenv("ATLAS_RIDGE_LAMBDA", "1e-6")),
    min_model_rows = as.integer(Sys.getenv("ATLAS_MIN_MODEL_ROWS", "20")),
    min_model_units = as.integer(Sys.getenv("ATLAS_MIN_MODEL_UNITS", "3")),
    null_sign_permutations = as.integer(Sys.getenv("ATLAS_NULL_SIGN_PERMUTATIONS", "49")),
    relaimpo_min_rows = as.integer(Sys.getenv("ATLAS_RELAIMPO_MIN_ROWS", "20")),
    relaimpo_min_units = as.integer(Sys.getenv("ATLAS_RELAIMPO_MIN_UNITS", "3")),
    render_figures = Sys.getenv("ATLAS_RENDER_FIGURES", "1") != "0",
    render_cell_figures = Sys.getenv("ATLAS_RENDER_CELL_FIGURES", "1") != "0",
    dry_run = Sys.getenv("ATLAS_DRY_RUN", "0") == "1"
  )
}

atlas_prepare_dirs <- function(cfg) {
  dirs <- c(cfg$output_dir, cfg$table_dir, cfg$figure_dir, cfg$checkpoint_dir, cfg$log_dir)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$figure_dir, "contributor_matrix"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$figure_dir, "dihedral"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$figure_dir, "secondary_structure"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$figure_dir, "shrinkage"), recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

atlas_manifest_path <- function(cfg, name) file.path(cfg$input_dir, name)
atlas_checkpoint <- function(cfg, name) file.path(cfg$checkpoint_dir, paste0(name, ".rds"))
atlas_table <- function(cfg, name) file.path(cfg$table_dir, name)

atlas_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  data.table::fwrite(x, tmp)
  file.rename(tmp, path)
  invisible(path)
}

atlas_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(x, tmp, compress = "xz")
  file.rename(tmp, path)
  invisible(path)
}

atlas_parse_shape <- function(shape) {
  as.integer(regmatches(shape, gregexpr("[0-9]+", shape))[[1]])
}

atlas_read_json <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

atlas_rows_schema <- function(rows_path) {
  cols <- strsplit(readLines(rows_path, n = 1L), ",", fixed = TRUE)[[1]]
  string_cols <- c(
    "global_row_key", "dataset_id", "protein_id", "pose_id", "conformation_id",
    "pose_kind", "split_group_id", "element", "residue_type", "iupac_name",
    "bmrb_name", "chemical_category", "ring_type", "polar_h_kind",
    "backbone_role", "prev_residue_type", "next_residue_type",
    "prev_residue_class", "next_residue_class", "terminal_state", "ss8", "ss3",
    "region_def_id", "rama_region_label", "rotamer_id", "rotamer_label",
    "folding_rule_id", "prop_cutoff_set_id"
  )
  int_cols <- c(
    "row_id", "frame_slot", "h5_row", "original_frame_index", "atom_index",
    "residue_index", "residue_number", "ring_membership_id", "in_aromatic_ring",
    "formal_charge", "is_exchangeable", "pre_proline", "is_pro", "is_gly",
    "phi_bin", "psi_bin", "omega_bin", "chi1_exists", "chi2_exists",
    "chi3_exists", "chi4_exists", "rama_region_id", "prop_self_or_bonded_atom_count",
    "prop_self_or_bonded_bond_count", "prop_has_self_or_bonded_driver",
    "prop_nearest_ring_id", "prop_nearest_charge_atom", "prop_nearest_bond_id",
    "prop_nearest_atom_index", "dft_present", "target_total_present",
    "target_dia_present", "target_para_present", "aimnet2_embedding_present",
    "prop_n_atoms_4A", "prop_n_rings_4A", "prop_n_charges_4A", "prop_n_bonds_4A",
    "prop_n_atoms_6A", "prop_n_rings_6A", "prop_n_charges_6A", "prop_n_bonds_6A",
    "prop_n_atoms_8A", "prop_n_rings_8A", "prop_n_charges_8A", "prop_n_bonds_8A",
    "prop_n_atoms_10A", "prop_n_rings_10A", "prop_n_charges_10A", "prop_n_bonds_10A",
    "locant", "branch_outer", "branch_inner", "di_index", "prochiral",
    "equivalence_class"
  )
  types <- lapply(cols, function(nm) {
    if (nm %in% string_cols) return(arrow::utf8())
    if (nm %in% int_cols) return(arrow::int64())
    arrow::float64()
  })
  do.call(arrow::schema, stats::setNames(types, cols))
}

atlas_sidecar_manifest <- function(cfg) {
  data.table::as.data.table(atlas_read_json(atlas_manifest_path(cfg, "sidecar_manifest.json"))$sidecars)
}

atlas_family_manifest <- function(cfg) {
  fam <- atlas_read_json(atlas_manifest_path(cfg, "contributor_family_manifest.json"))
  rows <- data.table::as.data.table(fam$rows)
  if (!"family_id" %in% names(rows)) rows[, family_id := ""]
  if (!"legal_permutation_block" %in% names(rows)) rows[, legal_permutation_block := ""]
  rows
}

atlas_read_rows <- function(cfg) {
  cp <- atlas_checkpoint(cfg, "rows_index_with_targets")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }

  atlas_msg("opening rows.csv through arrow")
  rows_path <- file.path(cfg$input_dir, "rows.csv")
  ds <- arrow::open_dataset(rows_path, format = "csv", schema = atlas_rows_schema(rows_path), skip = 1)
  needed <- c(
    "row_id", "dataset_id", "protein_id", "pose_id", "split_group_id",
    "frame_slot", "original_frame_index", "time_ps", "atom_index",
    "element", "residue_index", "residue_number", "residue_type",
    "iupac_name", "chemical_category", "backbone_role", "locant",
    "branch_outer", "branch_inner", "di_index", "prochiral",
    "equivalence_class", "ss8", "ss3", "phi", "psi", "omega",
    "phi_bin", "psi_bin", "omega_bin", "chi1_sin", "chi1_cos",
    "chi1_exists", "chi2_sin", "chi2_cos", "chi2_exists", "chi3_sin",
    "chi3_cos", "chi3_exists", "chi4_sin", "chi4_cos", "chi4_exists",
    "region_def_id", "rama_region_id", "rama_region_label",
    "rama_region_prob", "rotamer_id", "rotamer_label", "rotamer_prob",
    "folding_rule_id", "dft_present", "target_total_present",
    "target_dia_present", "target_para_present"
  )
  rows <- ds |>
    dplyr::select(dplyr::all_of(needed)) |>
    dplyr::collect()
  data.table::setDT(rows)
  data.table::setorder(rows, row_id)

  atlas_msg("loading ORCA scalar targets with RcppCNPy")
  target_dir <- file.path(cfg$input_dir, "targets")
  row_idx <- rows$row_id + 1L
  rows[, target_present_sidecar := as.logical(RcppCNPy::npyLoad(file.path(target_dir, "target_present.npy"))[row_idx])]
  rows[, sigma_iso := as.numeric(RcppCNPy::npyLoad(file.path(target_dir, "orca_total_T0.npy"))[row_idx])]
  dia <- as.numeric(RcppCNPy::npyLoad(file.path(target_dir, "orca_diamagnetic_T0.npy"))[row_idx])
  para <- as.numeric(RcppCNPy::npyLoad(file.path(target_dir, "orca_paramagnetic_T0.npy"))[row_idx])
  rows[, dft_floor_resid := sigma_iso - dia - para]
  rows[, supervised := target_present_sidecar & target_total_present == 1 & is.finite(sigma_iso)]

  make_safe <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "."
    x
  }
  rows[, nucleus_id := paste(
    make_safe(element), make_safe(iupac_name), make_safe(chemical_category),
    make_safe(backbone_role), make_safe(locant), make_safe(branch_outer),
    make_safe(branch_inner), make_safe(di_index), make_safe(prochiral),
    make_safe(equivalence_class),
    sep = "|"
  )]
  rows[, nucleus_label := paste0(make_safe(element), ":", make_safe(iupac_name))]
  rows[, residue_scope := make_safe(residue_type)]
  rows[, residue_nucleus_id := paste(residue_scope, nucleus_id, sep = "||")]
  rows[, pooled_nucleus_id := paste("ALL", nucleus_id, sep = "||")]
  rows[, physical_atom_key := paste(dataset_id, protein_id, pose_id, atom_index, sep = "|")]
  rows[, physical_residue_key := paste(dataset_id, protein_id, pose_id, residue_index, sep = "|")]
  rows[, time_block := data.table::fifelse(
    dataset_id == "1p9j",
    paste0("frame_block_", floor(as.numeric(original_frame_index) / 20)),
    protein_id
  )]
  rows[, fold_unit := data.table::fifelse(dataset_id == "720_static", protein_id, time_block)]
  rows[, fold_id := paste(dataset_id, fold_unit, sep = "::")]

  atlas_save_rds(rows, cp)
  atlas_write_csv(rows[, .(
    row_id, dataset_id, protein_id, pose_id, atom_index, element, residue_type,
    iupac_name, nucleus_id, residue_nucleus_id, pooled_nucleus_id, supervised
  )], atlas_table(cfg, "row_index_audit.csv"))
  rows
}

atlas_view_subset <- function(dt, view) {
  if (view == "within_1p9j") return(dt[dataset_id == "1p9j"])
  if (view == "between_720") return(dt[dataset_id == "720_static"])
  if (view == "pooled") return(dt)
  stop("unknown view: ", view, call. = FALSE)
}

atlas_ar1 <- function(y, id, ord) {
  o <- order(id, ord)
  y <- y[o]
  id <- id[o]
  same <- id[-1L] == id[-length(id)]
  if (sum(same & is.finite(y[-1L]) & is.finite(y[-length(y)])) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(y[-1L][same], y[-length(y)][same], use = "complete.obs"))
}

atlas_support_gate <- function(view, n_rows, n_proteins, n_physical_atoms, n_eff) {
  n_eff_safe <- ifelse(is.na(n_eff), n_rows, n_eff)
  if (n_rows <= 1L || n_physical_atoms <= 1L) return("singleton_shrinkage_only")
  if (view %in% c("between_720", "pooled") && n_proteins < 5L) return("shrinkage_only")
  if (view %in% c("within_1p9j", "pooled") && !is.na(n_eff) && n_eff < 20) return("descriptive_only")
  if (n_rows < 20L) return("descriptive_only")
  if ((view %in% c("between_720", "pooled") && n_proteins < 20L) || n_eff_safe < 100L) return("partial_pool_only")
  "fixed_effect_eligible"
}

atlas_make_support <- function(rows, cfg) {
  cp <- atlas_checkpoint(cfg, "analysis_cell_support")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }
  atlas_msg("building support table")
  sup_rows <- rows[supervised == TRUE]
  out <- list()
  k <- 1L
  for (view in c("within_1p9j", "between_720", "pooled")) {
    dtv <- atlas_view_subset(sup_rows, view)
    for (granularity in c("per_residue_type", "pooled_over_residue")) {
      if (granularity == "per_residue_type") {
        dtv[, `:=`(cell_id_tmp = residue_nucleus_id, residue_scope_tmp = residue_scope)]
      } else {
        dtv[, `:=`(cell_id_tmp = pooled_nucleus_id, residue_scope_tmp = "ALL")]
      }
      base <- dtv[, .(
        n_rows = .N,
        n_physical_atoms = data.table::uniqueN(physical_atom_key),
        n_residues = data.table::uniqueN(physical_residue_key),
        n_proteins = data.table::uniqueN(protein_id),
        n_frames = data.table::uniqueN(original_frame_index[dataset_id == "1p9j"]),
        n_datasets = data.table::uniqueN(dataset_id),
        sigma_mean = mean(sigma_iso, na.rm = TRUE),
        sigma_sd = stats::sd(sigma_iso, na.rm = TRUE),
        dft_floor_mse = mean(dft_floor_resid^2, na.rm = TRUE)
      ), by = .(cell_id = cell_id_tmp, residue_scope = residue_scope_tmp, nucleus_id, nucleus_label)]

      ar <- dtv[dataset_id == "1p9j", .(
        lag1_rho = atlas_ar1(sigma_iso, physical_atom_key, original_frame_index)
      ), by = .(cell_id = cell_id_tmp)]
      base <- merge(base, ar, by = "cell_id", all.x = TRUE)
      base[, n_eff := data.table::fifelse(
        is.na(lag1_rho), as.numeric(n_rows),
        pmax(1, as.numeric(n_rows) * pmax(0.001, 1 - lag1_rho) / pmax(0.001, 1 + lag1_rho))
      )]
      base[, `:=`(dataset_scope = view, granularity = granularity)]
      base[, support_gate := mapply(
        atlas_support_gate, dataset_scope, n_rows, n_proteins, n_physical_atoms, n_eff,
        USE.NAMES = FALSE
      )]
      out[[k]] <- base
      k <- k + 1L
    }
  }
  support <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  data.table::setcolorder(support, c(
    "dataset_scope", "granularity", "cell_id", "residue_scope", "nucleus_id",
    "nucleus_label", "support_gate"
  ))
  atlas_save_rds(support, cp)
  atlas_write_csv(support, atlas_table(cfg, "analysis_cell_support.csv"))

  thresholds <- data.table::data.table(
    gate = c("singleton_shrinkage_only", "shrinkage_only", "descriptive_only", "partial_pool_only", "fixed_effect_eligible"),
    rule = c(
      "n_rows <= 1 or n_physical_atoms <= 1",
      "720/pooled support has fewer than 5 proteins",
      "n_rows < 20 or 1P9J effective support < 20",
      "720/pooled support has fewer than 20 proteins or effective support < 100",
      "passes all support gates"
    )
  )
  atlas_write_csv(thresholds, atlas_table(cfg, "support_gate_thresholds.csv"))
  support
}

atlas_contributors <- function(cfg) {
  sidecars <- atlas_sidecar_manifest(cfg)
  fam <- atlas_family_manifest(cfg)
  contrib <- sidecars[route == "measure"]
  contrib[, row_aligned := vapply(shape, function(s) {
    sh <- atlas_parse_shape(s)
    length(sh) >= 1L && sh[1L] == 1744962L
  }, logical(1))]
  contrib <- contrib[row_aligned == TRUE]
  contrib <- merge(
    contrib, fam,
    by.x = "stem", by.y = "contributor_stem",
    all.x = TRUE, sort = FALSE
  )
  for (nm in c("family_id", "legal_permutation_block")) {
    if (!nm %in% names(contrib)) contrib[, (nm) := ""]
    contrib[is.na(get(nm)), (nm) := ""]
  }
  contrib[, contributor_id := stem]
  contrib[]
}

atlas_as_matrix <- function(x) {
  d <- dim(x)
  if (is.null(d)) return(matrix(as.numeric(x), ncol = 1L))
  if (length(d) == 1L) return(matrix(as.numeric(x), ncol = 1L))
  matrix(as.numeric(x), nrow = d[1L], ncol = prod(d[-1L]))
}

atlas_load_contributor_block <- function(cfg, meta, row_ids, max_cols) {
  path <- file.path(cfg$input_dir, meta$sidecar)
  atlas_msg("loading contributor ", meta$stem, " from ", meta$sidecar)
  arr <- RcppCNPy::npyLoad(path)
  x <- atlas_as_matrix(arr)
  x <- x[row_ids + 1L, , drop = FALSE]
  presence_path <- sub("\\.npy$", "_presence.npy", path)
  if (file.exists(presence_path)) {
    present <- as.logical(RcppCNPy::npyLoad(presence_path)[row_ids + 1L])
    x[!present, ] <- NA_real_
  }
  atlas_reduce_block(x, meta, max_cols)
}

atlas_reduce_block <- function(x, meta, max_cols) {
  note <- "as_emitted"
  keep <- vapply(seq_len(ncol(x)), function(j) {
    z <- x[, j]
    any(is.finite(z)) && stats::sd(z[is.finite(z)]) > 0
  }, logical(1))
  x <- x[, keep, drop = FALSE]
  if (!ncol(x)) {
    return(list(x = matrix(NA_real_, nrow = nrow(x), ncol = 1L), note = "no_finite_nonconstant_columns"))
  }
  if (ncol(x) > max_cols) {
    note <- paste0("all_rows_predictor_pca_to_", max_cols)
    mu <- colMeans(x, na.rm = TRUE)
    sdv <- apply(x, 2L, stats::sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    z <- sweep(sweep(x, 2L, mu, "-"), 2L, sdv, "/")
    z[!is.finite(z)] <- 0
    eig <- eigen(crossprod(z) / max(1, nrow(z) - 1), symmetric = TRUE)
    k <- min(max_cols, ncol(z), ncol(eig$vectors))
    x <- z %*% eig$vectors[, seq_len(k), drop = FALSE]
    colnames(x) <- paste0("PC", seq_len(k))
  }
  list(x = x, note = note)
}

atlas_scale_train <- function(x_train, x_test) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sdv <- apply(x_train, 2L, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  xt <- sweep(sweep(x_train, 2L, mu, "-"), 2L, sdv, "/")
  xv <- sweep(sweep(x_test, 2L, mu, "-"), 2L, sdv, "/")
  xt[!is.finite(xt)] <- 0
  xv[!is.finite(xv)] <- 0
  list(train = xt, test = xv)
}

atlas_ridge_beta <- function(x, y, lambda) {
  p <- ncol(x)
  xtx <- crossprod(x)
  diag(xtx) <- diag(xtx) + lambda
  rhs <- crossprod(x, y)
  out <- tryCatch(
    solve(xtx, rhs),
    error = function(e) qr.solve(xtx, rhs)
  )
  matrix(out, nrow = p)
}

atlas_crossfit_block <- function(y, x, fold, floor_resid, baseline_factor = NULL,
                                 lambda = 1e-6, n_perm = 49L, min_rows = 20L,
                                 min_units = 3L) {
  ok <- is.finite(y) & is.finite(floor_resid) & !is.na(fold) & apply(x, 1L, function(z) all(is.finite(z)))
  if (!is.null(baseline_factor)) ok <- ok & !is.na(baseline_factor)
  y <- y[ok]
  x <- x[ok, , drop = FALSE]
  fold <- as.character(fold[ok])
  floor_resid <- floor_resid[ok]
  baseline_factor <- if (is.null(baseline_factor)) rep("baseline", length(y)) else as.character(baseline_factor[ok])
  folds <- unique(fold)
  if (length(y) < min_rows || length(folds) < min_units) {
    return(list(
      n_model_rows = length(y), n_fold_units = length(folds), R_oof = NA_real_,
      R2_oof = NA_real_, R2_oof_floor = NA_real_, dft_floor_mse = mean(floor_resid^2),
      permutation_null_p = NA_real_, fold_scheme = "insufficient_support"
    ))
  }
  pred <- rep(NA_real_, length(y))
  base_oof <- rep(NA_real_, length(y))
  for (f in folds) {
    test <- fold == f
    train <- !test
    train_base <- data.table::data.table(y = y[train], b = baseline_factor[train])[, .(base = mean(y)), by = b]
    base_map <- stats::setNames(train_base$base, train_base$b)
    global_base <- mean(y[train])
    base_oof[test] <- base_map[baseline_factor[test]]
    base_oof[test][is.na(base_oof[test])] <- global_base
    y_res <- y[train] - (base_map[baseline_factor[train]] %||% global_base)
    missing_base <- is.na(y_res)
    if (any(missing_base)) y_res[missing_base] <- y[train][missing_base] - global_base
    sc <- atlas_scale_train(x[train, , drop = FALSE], x[test, , drop = FALSE])
    beta <- atlas_ridge_beta(sc$train, y_res, lambda)
    pred[test] <- base_oof[test] + as.numeric(sc$test %*% beta)
  }
  valid <- is.finite(pred) & is.finite(base_oof)
  if (sum(valid) < min_rows) {
    return(list(
      n_model_rows = sum(valid), n_fold_units = length(folds), R_oof = NA_real_,
      R2_oof = NA_real_, R2_oof_floor = NA_real_, dft_floor_mse = mean(floor_resid^2),
      permutation_null_p = NA_real_, fold_scheme = "failed_fit"
    ))
  }
  y <- y[valid]
  pred <- pred[valid]
  base_oof <- base_oof[valid]
  fold <- fold[valid]
  floor_resid <- floor_resid[valid]
  contrib_pred <- pred - base_oof
  y_res <- y - base_oof
  r <- suppressWarnings(stats::cor(y_res, contrib_pred, use = "complete.obs"))
  sse <- sum((y - pred)^2)
  sst <- sum((y - base_oof)^2)
  floor_mse <- mean(floor_resid^2, na.rm = TRUE)
  r2 <- 1 - sse / max(sst, .Machine$double.eps)
  r2_floor <- 1 - max(sse - length(y) * floor_mse, 0) / max(sst - length(y) * floor_mse, .Machine$double.eps)
  p_null <- NA_real_
  if (n_perm > 0L && is.finite(r)) {
    set.seed(7349L)
    fold_levels <- unique(fold)
    null_r <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      signs <- sample(c(-1, 1), length(fold_levels), replace = TRUE)
      names(signs) <- fold_levels
      null_r[i] <- suppressWarnings(stats::cor(y_res, contrib_pred * signs[fold], use = "complete.obs"))
    }
    p_null <- (1 + sum(abs(null_r) >= abs(r), na.rm = TRUE)) / (n_perm + 1)
  }
  list(
    n_model_rows = length(y), n_fold_units = length(folds), R_oof = r,
    R2_oof = r2, R2_oof_floor = r2_floor, dft_floor_mse = floor_mse,
    permutation_null_p = p_null, fold_scheme = "leave_fold_unit_out"
  )
}

`%||%` <- function(a, b) {
  ifelse(is.na(a), b, a)
}

atlas_solve <- function(a, b) {
  tryCatch(
    solve(a, b),
    error = function(e1) {
      tryCatch(
        qr.solve(a, b),
        error = function(e2) {
          scale <- max(abs(diag(a)), 1, na.rm = TRUE)
          solve(a + diag(scale * 1e-8, nrow(a)), b)
        }
      )
    }
  )
}

atlas_oof_refit_predictions <- function(design, y, fold, penalty_diag) {
  n <- nrow(design)
  y_was_vector <- is.null(dim(y))
  y_mat <- if (y_was_vector) matrix(y, ncol = 1L) else as.matrix(y)
  pred <- matrix(NA_real_, nrow = n, ncol = ncol(y_mat))
  xtx_full <- crossprod(design)
  rhs_full <- crossprod(design, y_mat)
  penalized_full <- xtx_full
  diag(penalized_full) <- diag(penalized_full) + penalty_diag
  split_idx <- split(seq_len(n), fold, drop = TRUE)
  for (idx in split_idx) {
    xg <- design[idx, , drop = FALSE]
    yg <- y_mat[idx, , drop = FALSE]
    a <- penalized_full - crossprod(xg)
    rhs <- rhs_full - crossprod(xg, yg)
    beta <- atlas_solve(a, rhs)
    pred[idx, ] <- xg %*% beta
  }
  if (y_was_vector) as.numeric(pred[, 1L]) else pred
}

atlas_baseline_design <- function(baseline_factor) {
  baseline_factor <- as.factor(baseline_factor)
  if (nlevels(baseline_factor) < 2L) {
    out <- matrix(1, nrow = length(baseline_factor), ncol = 1L)
    colnames(out) <- "baseline_intercept"
    return(out)
  }
  stats::model.matrix(~ 0 + baseline_factor)
}

atlas_crossfit_block <- function(y, x, fold, floor_resid, baseline_factor = NULL,
                                 lambda = 1e-6, n_perm = 49L, min_rows = 20L,
                                 min_units = 3L) {
  row_ok <- is.finite(y) & is.finite(floor_resid) & !is.na(fold)
  x_ok <- apply(x, 1L, function(z) all(is.finite(z)))
  if (!is.null(baseline_factor)) row_ok <- row_ok & !is.na(baseline_factor)
  ok <- row_ok & x_ok
  y <- y[ok]
  x <- x[ok, , drop = FALSE]
  fold <- as.character(fold[ok])
  floor_resid <- floor_resid[ok]
  baseline_factor <- if (is.null(baseline_factor)) rep("baseline", length(y)) else as.character(baseline_factor[ok])
  folds <- unique(fold)
  if (length(y) < min_rows || length(folds) < min_units) {
    return(list(
      n_model_rows = length(y), n_fold_units = length(folds), R_oof = NA_real_,
      R2_oof = NA_real_, R2_oof_floor = NA_real_, dft_floor_mse = mean(floor_resid^2),
      permutation_null_p = NA_real_, fold_scheme = "insufficient_support"
    ))
  }

  mu <- colMeans(x, na.rm = TRUE)
  sdv <- apply(x, 2L, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  xs <- sweep(sweep(x, 2L, mu, "-"), 2L, sdv, "/")
  xs[!is.finite(xs)] <- 0
  base_design <- atlas_baseline_design(baseline_factor)
  full_design <- cbind(base_design, xs)
  pred_base <- atlas_oof_refit_predictions(
    base_design, y, fold,
    penalty_diag = rep(0, ncol(base_design))
  )
  pred_full <- atlas_oof_refit_predictions(
    full_design, y, fold,
    penalty_diag = c(rep(0, ncol(base_design)), rep(lambda, ncol(xs)))
  )
  valid <- is.finite(pred_base) & is.finite(pred_full)
  if (sum(valid) < min_rows) {
    return(list(
      n_model_rows = sum(valid), n_fold_units = length(folds), R_oof = NA_real_,
      R2_oof = NA_real_, R2_oof_floor = NA_real_, dft_floor_mse = mean(floor_resid^2),
      permutation_null_p = NA_real_, fold_scheme = "failed_press_refit"
    ))
  }
  y <- y[valid]
  pred_base <- pred_base[valid]
  pred_full <- pred_full[valid]
  floor_resid <- floor_resid[valid]
  fold <- fold[valid]
  y_res <- y - pred_base
  contrib_pred <- pred_full - pred_base
  r <- suppressWarnings(stats::cor(y_res, contrib_pred, use = "complete.obs"))
  sse <- sum((y - pred_full)^2)
  sst <- sum((y - pred_base)^2)
  floor_mse <- mean(floor_resid^2, na.rm = TRUE)
  r2 <- 1 - sse / max(sst, .Machine$double.eps)
  r2_floor <- 1 - max(sse - length(y) * floor_mse, 0) / max(sst - length(y) * floor_mse, .Machine$double.eps)
  p_null <- NA_real_
  if (n_perm > 0L && is.finite(r)) {
    set.seed(7349L)
    fold_levels <- unique(fold)
    null_r <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      signs <- sample(c(-1, 1), length(fold_levels), replace = TRUE)
      names(signs) <- fold_levels
      null_r[i] <- suppressWarnings(stats::cor(y_res, contrib_pred * signs[fold], use = "complete.obs"))
    }
    p_null <- (1 + sum(abs(null_r) >= abs(r), na.rm = TRUE)) / (n_perm + 1)
  }
  list(
    n_model_rows = length(y), n_fold_units = length(folds), R_oof = r,
    R2_oof = r2, R2_oof_floor = r2_floor, dft_floor_mse = floor_mse,
    permutation_null_p = p_null, fold_scheme = "leave_fold_unit_out_crossproduct_refit"
  )
}

atlas_crossfit_t2 <- function(y, x, fold, lambda = 1e-6, min_rows = 20L, min_units = 3L) {
  ok <- stats::complete.cases(y) & !is.na(fold) & apply(x, 1L, function(z) all(is.finite(z)))
  y <- y[ok, , drop = FALSE]
  x <- x[ok, , drop = FALSE]
  fold <- as.character(fold[ok])
  folds <- unique(fold)
  if (nrow(y) < min_rows || length(folds) < min_units) {
    return(list(n_model_rows = nrow(y), n_fold_units = length(folds), vector_R2_oof = NA_real_,
                angular_cosine = NA_real_, norm_recovery_slope = NA_real_, norm_rmse = NA_real_))
  }
  mu <- colMeans(x, na.rm = TRUE)
  sdv <- apply(x, 2L, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  xs <- sweep(sweep(x, 2L, mu, "-"), 2L, sdv, "/")
  xs[!is.finite(xs)] <- 0
  base_design <- matrix(1, nrow = nrow(y), ncol = 1L)
  colnames(base_design) <- "baseline_intercept"
  full_design <- cbind(base_design, xs)
  base <- atlas_oof_refit_predictions(base_design, y, fold, penalty_diag = 0)
  pred <- atlas_oof_refit_predictions(full_design, y, fold, penalty_diag = c(0, rep(lambda, ncol(xs))))
  okp <- stats::complete.cases(pred)
  y <- y[okp, , drop = FALSE]
  pred <- pred[okp, , drop = FALSE]
  base <- base[okp, , drop = FALSE]
  sse <- sum((y - pred)^2)
  sst <- sum((y - base)^2)
  yn <- sqrt(rowSums(y^2))
  pn <- sqrt(rowSums(pred^2))
  cosine <- rowSums(y * pred) / pmax(yn * pn, .Machine$double.eps)
  slope <- as.numeric(stats::coef(stats::lm(pn ~ 0 + yn))[1])
  list(
    n_model_rows = nrow(y),
    n_fold_units = length(folds),
    vector_R2_oof = 1 - sse / max(sst, .Machine$double.eps),
    angular_cosine = mean(cosine[is.finite(cosine)], na.rm = TRUE),
    norm_recovery_slope = slope,
    norm_rmse = sqrt(mean((pn - yn)^2, na.rm = TRUE))
  )
}

atlas_score_contributors <- function(rows, support, cfg) {
  cp <- atlas_checkpoint(cfg, "contributor_score_by_nucleus")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    scores <- readRDS(cp)
    if (!"contributor_id" %in% names(scores) || all(is.na(scores$contributor_id))) {
      csv <- atlas_table(cfg, "contributor_score_by_nucleus.csv")
      if (file.exists(csv)) {
        atlas_msg("repairing contributor score checkpoint from CSV")
        scores <- data.table::fread(csv)
        atlas_save_rds(scores, cp)
      }
    }
    return(scores)
  }

  atlas_msg("scoring contributors")
  contrib <- atlas_contributors(cfg)
  supervised_rows <- rows[supervised == TRUE]
  row_ids <- supervised_rows$row_id
  score_files <- character()
  prep_rows <- list()

  for (i in seq_len(nrow(contrib))) {
    meta <- contrib[i]
    out_file <- atlas_table(cfg, file.path("contributor_scores", paste0(meta$stem, ".csv")))
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(out_file)) {
      score_files <- c(score_files, out_file)
      next
    }
    block <- atlas_load_contributor_block(cfg, meta, row_ids, cfg$max_block_cols)
    x <- block$x
    prep_rows[[length(prep_rows) + 1L]] <- data.table::data.table(
      contributor_id = meta$stem,
      original_component_layout = meta$component_layout,
      scored_columns = ncol(x),
      preprocessing = block$note
    )
    res <- list()
    rix <- 1L
    for (view in c("within_1p9j", "between_720", "pooled")) {
      idx_view <- which(
        if (view == "within_1p9j") supervised_rows$dataset_id == "1p9j"
        else if (view == "between_720") supervised_rows$dataset_id == "720_static"
        else rep(TRUE, nrow(supervised_rows))
      )
      dtv <- supervised_rows[idx_view]
      xv <- x[idx_view, , drop = FALSE]
      for (granularity in c("per_residue_type", "pooled_over_residue")) {
        this_granularity <- granularity
        cell_vec <- if (granularity == "per_residue_type") dtv$residue_nucleus_id else dtv$pooled_nucleus_id
        sup <- support[dataset_scope == view & granularity == this_granularity]
        split_idx <- split(seq_along(cell_vec), cell_vec, drop = TRUE)
        for (cell_id in sup$cell_id) {
          this_cell_id <- cell_id
          srow <- sup[cell_id == this_cell_id]
          if (!length(split_idx[[cell_id]]) || srow$support_gate %in% c("singleton_shrinkage_only", "shrinkage_only", "descriptive_only")) {
            res[[rix]] <- data.table::data.table(
              dataset_scope = view, granularity = granularity, cell_id = cell_id,
              residue_scope = srow$residue_scope, nucleus_id = srow$nucleus_id,
              nucleus_label = srow$nucleus_label, contributor_id = meta$stem,
              mechanism = meta$mechanism, units = meta$units, component_layout = meta$component_layout,
              family_id = meta$family_id, legal_permutation_block = meta$legal_permutation_block,
              support_gate = srow$support_gate, n_rows = srow$n_rows,
              n_proteins = srow$n_proteins, n_eff = srow$n_eff,
              n_model_rows = 0L, n_fold_units = 0L, R_oof = NA_real_,
              R2_oof = NA_real_, R2_oof_floor = NA_real_,
              dft_floor_mse = srow$dft_floor_mse, permutation_null_p = NA_real_,
              fold_scheme = "support_gated_descriptive_only"
            )
            rix <- rix + 1L
            next
          }
          jj <- split_idx[[cell_id]]
          baseline_factor <- if (view == "pooled") dtv$dataset_id[jj] else NULL
          sc <- atlas_crossfit_block(
            y = dtv$sigma_iso[jj], x = xv[jj, , drop = FALSE],
            fold = dtv$fold_id[jj], floor_resid = dtv$dft_floor_resid[jj],
            baseline_factor = baseline_factor, lambda = cfg$ridge_lambda,
            n_perm = cfg$null_sign_permutations, min_rows = cfg$min_model_rows,
            min_units = cfg$min_model_units
          )
          res[[rix]] <- data.table::data.table(
            dataset_scope = view, granularity = granularity, cell_id = cell_id,
            residue_scope = srow$residue_scope, nucleus_id = srow$nucleus_id,
            nucleus_label = srow$nucleus_label, contributor_id = meta$stem,
            mechanism = meta$mechanism, units = meta$units, component_layout = meta$component_layout,
            family_id = meta$family_id, legal_permutation_block = meta$legal_permutation_block,
            support_gate = srow$support_gate, n_rows = srow$n_rows,
            n_proteins = srow$n_proteins, n_eff = srow$n_eff,
            n_model_rows = sc$n_model_rows, n_fold_units = sc$n_fold_units,
            R_oof = sc$R_oof, R2_oof = sc$R2_oof, R2_oof_floor = sc$R2_oof_floor,
            dft_floor_mse = sc$dft_floor_mse, permutation_null_p = sc$permutation_null_p,
            fold_scheme = sc$fold_scheme
          )
          rix <- rix + 1L
        }
      }
    }
    one <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
    atlas_write_csv(one, out_file)
    score_files <- c(score_files, out_file)
    rm(block, x, one)
    gc()
  }
  scores <- data.table::rbindlist(lapply(score_files, data.table::fread), use.names = TRUE, fill = TRUE)
  scores[, q_value := stats::p.adjust(permutation_null_p, method = "BH"), by = .(dataset_scope, granularity, mechanism)]
  atlas_write_csv(scores, atlas_table(cfg, "contributor_score_by_nucleus.csv"))
  if (length(prep_rows)) {
    prep <- data.table::rbindlist(prep_rows, use.names = TRUE, fill = TRUE)
  } else {
    prep <- data.table::rbindlist(lapply(score_files, function(f) NULL), fill = TRUE)
  }
  if (nrow(prep)) atlas_write_csv(prep, atlas_table(cfg, "contributor_block_preprocessing.csv"))
  atlas_save_rds(scores, cp)
  scores
}

atlas_score_t2 <- function(rows, support, cfg) {
  cp <- atlas_checkpoint(cfg, "tensor_T2_equivariant_scores")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }
  atlas_msg("scoring multivariate T2 contributors")
  contrib <- atlas_contributors(cfg)
  tensor <- contrib[component_layout %in% c("T2_5", "tensor9_raw")]
  supervised_rows <- rows[supervised == TRUE]
  row_ids <- supervised_rows$row_id
  y_all <- atlas_as_matrix(RcppCNPy::npyLoad(file.path(cfg$input_dir, "targets", "orca_total_T2.npy")))[row_ids + 1L, , drop = FALSE]
  res <- list()
  rix <- 1L
  for (i in seq_len(nrow(tensor))) {
    meta <- tensor[i]
    raw <- RcppCNPy::npyLoad(file.path(cfg$input_dir, meta$sidecar))
    x <- atlas_as_matrix(raw)[row_ids + 1L, , drop = FALSE]
    if (meta$component_layout == "tensor9_raw" && ncol(x) >= 9L) x <- x[, 5:9, drop = FALSE]
    x <- atlas_reduce_block(x, meta, min(cfg$max_block_cols, 5L))$x
    for (view in c("within_1p9j", "between_720", "pooled")) {
      idx_view <- which(
        if (view == "within_1p9j") supervised_rows$dataset_id == "1p9j"
        else if (view == "between_720") supervised_rows$dataset_id == "720_static"
        else rep(TRUE, nrow(supervised_rows))
      )
      dtv <- supervised_rows[idx_view]
      xv <- x[idx_view, , drop = FALSE]
      yv <- y_all[idx_view, , drop = FALSE]
      for (granularity in c("per_residue_type", "pooled_over_residue")) {
        this_granularity <- granularity
        cell_vec <- if (granularity == "per_residue_type") dtv$residue_nucleus_id else dtv$pooled_nucleus_id
        sup <- support[dataset_scope == view & granularity == this_granularity]
        split_idx <- split(seq_along(cell_vec), cell_vec, drop = TRUE)
        for (cell_id in sup$cell_id) {
          this_cell_id <- cell_id
          srow <- sup[cell_id == this_cell_id]
          jj <- split_idx[[cell_id]]
          if (!length(jj) || srow$support_gate %in% c("singleton_shrinkage_only", "shrinkage_only", "descriptive_only")) {
            t2 <- list(n_model_rows = 0L, n_fold_units = 0L, vector_R2_oof = NA_real_,
                       angular_cosine = NA_real_, norm_recovery_slope = NA_real_, norm_rmse = NA_real_)
          } else {
            t2 <- atlas_crossfit_t2(yv[jj, , drop = FALSE], xv[jj, , drop = FALSE], dtv$fold_id[jj],
                                    cfg$ridge_lambda, cfg$min_model_rows, cfg$min_model_units)
          }
          res[[rix]] <- data.table::data.table(
            dataset_scope = view, granularity = granularity, cell_id = cell_id,
            residue_scope = srow$residue_scope, nucleus_id = srow$nucleus_id,
            nucleus_label = srow$nucleus_label, tensor_contributor_block = meta$stem,
            mechanism = meta$mechanism, tensor_frame = "lab", valid_for_T2_model = TRUE,
            support_gate = srow$support_gate, n_rows = srow$n_rows, n_proteins = srow$n_proteins,
            n_eff = srow$n_eff, n_model_rows = t2$n_model_rows, n_fold_units = t2$n_fold_units,
            vector_R2_oof = t2$vector_R2_oof, angular_cosine = t2$angular_cosine,
            norm_recovery_slope = t2$norm_recovery_slope, norm_rmse = t2$norm_rmse
          )
          rix <- rix + 1L
        }
      }
    }
    rm(raw, x)
    gc()
  }
  out <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  atlas_write_csv(out, atlas_table(cfg, "tensor_T2_equivariant_scores.csv"))
  atlas_save_rds(out, cp)
  out
}

atlas_residualize_by_factor <- function(z, f) {
  if (length(unique(f[!is.na(f)])) < 2L) return(as.numeric(z - mean(z, na.rm = TRUE)))
  stats::resid(stats::lm(z ~ factor(f), na.action = stats::na.exclude))
}

atlas_relaimpo_family_lmg <- function(rows, support, cfg) {
  cp <- atlas_checkpoint(cfg, "contributor_family_relaimpo_lmg")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }

  atlas_msg("running relaimpo LMG for small manifest-governed families")
  contrib <- atlas_contributors(cfg)
  family_members <- contrib[
    family_id != "" & legal_permutation_block != "",
    .(family_size = .N),
    by = .(family_id, legal_permutation_block)
  ][family_size >= 2L & family_size <= 6L]
  if (!nrow(family_members)) {
    out <- data.table::data.table()
    atlas_write_csv(out, atlas_table(cfg, "contributor_family_relaimpo_lmg.csv"))
    atlas_save_rds(out, cp)
    return(out)
  }

  supervised_rows <- rows[supervised == TRUE]
  row_ids <- supervised_rows$row_id
  out <- list()
  rix <- 1L

  for (fg in seq_len(nrow(family_members))) {
    fam <- family_members[fg]
    members <- contrib[
      family_id == fam$family_id &
        legal_permutation_block == fam$legal_permutation_block
    ]
    if (nrow(members) < 2L) next
    atlas_msg(
      "relaimpo family ", fam$family_id, " / ", fam$legal_permutation_block,
      " with ", nrow(members), " contributors"
    )
    x_list <- list()
    name_map <- data.table::data.table()
    for (mi in seq_len(nrow(members))) {
      meta <- members[mi]
      blk <- atlas_load_contributor_block(cfg, meta, row_ids, 1L)
      v <- as.numeric(blk$x[, 1L])
      safe <- make.names(meta$stem, unique = TRUE)
      x_list[[safe]] <- v
      name_map <- data.table::rbindlist(list(
        name_map,
        data.table::data.table(variable = safe, contributor_id = meta$stem)
      ), use.names = TRUE, fill = TRUE)
      rm(blk)
      gc()
    }
    x_dt <- data.table::as.data.table(x_list)

    for (view in c("within_1p9j", "between_720", "pooled")) {
      idx_view <- which(
        if (view == "within_1p9j") supervised_rows$dataset_id == "1p9j"
        else if (view == "between_720") supervised_rows$dataset_id == "720_static"
        else rep(TRUE, nrow(supervised_rows))
      )
      dtv <- supervised_rows[idx_view]
      xv <- x_dt[idx_view]
      for (granularity in c("per_residue_type", "pooled_over_residue")) {
        this_granularity <- granularity
        cell_vec <- if (granularity == "per_residue_type") dtv$residue_nucleus_id else dtv$pooled_nucleus_id
        sup <- support[dataset_scope == view & granularity == this_granularity]
        split_idx <- split(seq_along(cell_vec), cell_vec, drop = TRUE)
        for (cell_id in sup$cell_id) {
          this_cell_id <- cell_id
          srow <- sup[cell_id == this_cell_id]
          jj <- split_idx[[this_cell_id]]
          emit_na <- function(status) {
            data.table::data.table(
              dataset_scope = view, granularity = granularity, cell_id = this_cell_id,
              residue_scope = srow$residue_scope, nucleus_id = srow$nucleus_id,
              nucleus_label = srow$nucleus_label, family_id = fam$family_id,
              legal_permutation_block = fam$legal_permutation_block,
              contributor_id = members$stem, n_rows = srow$n_rows,
              n_proteins = srow$n_proteins, n_eff = srow$n_eff,
              n_model_rows = 0L, n_fold_units = 0L, family_model_R2 = NA_real_,
              lmg_R2 = NA_real_, first_R2 = NA_real_, last_unique_R2 = NA_real_,
              lmg_share_of_family_model = NA_real_, relaimpo_status = status,
              partition_claim_gate = "linear_model_empirical_variance_not_physical_independence"
            )
          }
          if (!length(jj) || srow$support_gate %in% c("singleton_shrinkage_only", "shrinkage_only", "descriptive_only")) {
            out[[rix]] <- emit_na("support_gated_descriptive_only")
            rix <- rix + 1L
            next
          }
          x_cell <- as.data.frame(xv[jj])
          y <- dtv$sigma_iso[jj]
          fold <- dtv$fold_id[jj]
          ok <- is.finite(y) & !is.na(fold) & stats::complete.cases(x_cell)
          if (sum(ok) < cfg$relaimpo_min_rows || data.table::uniqueN(fold[ok]) < cfg$relaimpo_min_units) {
            out[[rix]] <- emit_na("insufficient_support")
            rix <- rix + 1L
            next
          }
          x_cell <- x_cell[ok, , drop = FALSE]
          y <- y[ok]
          fold <- fold[ok]
          if (view == "pooled") {
            y <- atlas_residualize_by_factor(y, dtv$dataset_id[jj][ok])
            x_cell[] <- lapply(x_cell, atlas_residualize_by_factor, f = dtv$dataset_id[jj][ok])
          } else {
            y <- as.numeric(y - mean(y, na.rm = TRUE))
          }
          x_cell[] <- lapply(x_cell, function(z) {
            z <- as.numeric(z)
            z <- z - mean(z, na.rm = TRUE)
            s <- stats::sd(z, na.rm = TRUE)
            if (!is.finite(s) || s == 0) return(rep(NA_real_, length(z)))
            z / s
          })
          nonconstant <- vapply(x_cell, function(z) all(is.finite(z)) && stats::sd(z) > 0, logical(1))
          x_cell <- x_cell[, nonconstant, drop = FALSE]
          if (ncol(x_cell) < 2L) {
            out[[rix]] <- emit_na("fewer_than_two_nonconstant_family_members")
            rix <- rix + 1L
            next
          }
          df <- data.frame(y = as.numeric(y), x_cell, check.names = FALSE)
          fit <- tryCatch(stats::lm(y ~ ., data = df), error = identity)
          rel <- if (inherits(fit, "error")) fit else tryCatch(
            relaimpo::calc.relimp(fit, type = c("lmg", "last", "first"), rela = FALSE),
            error = identity
          )
          if (inherits(rel, "error")) {
            out[[rix]] <- emit_na(paste("relaimpo_error", conditionMessage(rel)))
            rix <- rix + 1L
            next
          }
          lmg <- rel@lmg
          last <- rel@last
          first <- rel@first
          rel_dt <- data.table::data.table(
            variable = names(lmg),
            lmg_R2 = as.numeric(lmg),
            first_R2 = as.numeric(first[names(lmg)]),
            last_unique_R2 = as.numeric(last[names(lmg)])
          )
          rel_dt <- merge(rel_dt, name_map, by = "variable", all.x = TRUE, sort = FALSE)
          rel_dt[, lmg_share_of_family_model := lmg_R2 / max(rel@R2, .Machine$double.eps)]
          out[[rix]] <- data.table::data.table(
            dataset_scope = view, granularity = granularity, cell_id = this_cell_id,
            residue_scope = srow$residue_scope, nucleus_id = srow$nucleus_id,
            nucleus_label = srow$nucleus_label, family_id = fam$family_id,
            legal_permutation_block = fam$legal_permutation_block,
            contributor_id = rel_dt$contributor_id, n_rows = srow$n_rows,
            n_proteins = srow$n_proteins, n_eff = srow$n_eff,
            n_model_rows = nrow(df), n_fold_units = data.table::uniqueN(fold),
            family_model_R2 = as.numeric(rel@R2),
            lmg_R2 = rel_dt$lmg_R2, first_R2 = rel_dt$first_R2,
            last_unique_R2 = rel_dt$last_unique_R2,
            lmg_share_of_family_model = rel_dt$lmg_share_of_family_model,
            relaimpo_status = "ok",
            partition_claim_gate = "relaimpo_lmg_empirical_variance_not_physical_independence"
          )
          rix <- rix + 1L
        }
      }
    }
    rm(x_dt, x_list)
    gc()
  }

  out <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  atlas_write_csv(out, atlas_table(cfg, "contributor_family_relaimpo_lmg.csv"))
  atlas_save_rds(out, cp)
  out
}

atlas_independence_partition <- function(scores, rows = NULL, support = NULL, cfg) {
  cp <- atlas_checkpoint(cfg, "contributor_independence_partition")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }
  atlas_msg("building linear-independence partition table from legal manifest blocks")
  part <- copy(scores)
  for (nm in c("family_id", "legal_permutation_block", "contributor_id")) {
    if (!nm %in% names(part)) part[, (nm) := ""]
    part[is.na(get(nm)), (nm) := ""]
  }
  part[, contributor_block := data.table::fifelse(
    legal_permutation_block != "", legal_permutation_block,
    data.table::fifelse(family_id != "", family_id, contributor_id)
  )]
  part[, block_R2 := max(R2_oof_floor, na.rm = TRUE), by = .(dataset_scope, granularity, cell_id, contributor_block)]
  part[!is.finite(block_R2), block_R2 := NA_real_]
  part[, shared_family_R2 := fifelse(
    family_id == "", 0,
    max(R2_oof_floor[family_id == family_id[1L]], na.rm = TRUE)
  ), by = .(dataset_scope, granularity, cell_id, family_id)]
  part[!is.finite(shared_family_R2), shared_family_R2 := NA_real_]
  part[, unique_R2 := pmax(0, R2_oof_floor - pmax(0, shared_family_R2 - R2_oof_floor, na.rm = TRUE), na.rm = TRUE)]
  part[, shared_R2 := pmax(0, R2_oof_floor - unique_R2)]
  part[, partition_method := data.table::fifelse(
    family_id %in% c("ring_current_kernel_collinear", "ring_current_bs_total_parts", "ring_current_hm_total_parts"),
    "family_manifest_linear_commonality_candidate",
    "drop_one_block_oof_delta_proxy"
  )]
  out <- part[, .(
    block_R2 = max(block_R2, na.rm = TRUE),
    unique_R2 = max(unique_R2, na.rm = TRUE),
    shared_R2 = max(shared_R2, na.rm = TRUE),
    delta_R2_drop_block = max(R2_oof_floor, na.rm = TRUE),
    permutation_null_p = min(permutation_null_p, na.rm = TRUE),
    q_value = min(q_value, na.rm = TRUE),
    support_gate = support_gate[which.max(abs(R_oof %||% 0))][1L],
    n_rows = max(n_rows, na.rm = TRUE),
    n_proteins = max(n_proteins, na.rm = TRUE),
    n_eff = max(n_eff, na.rm = TRUE),
    partition_method = partition_method[1L],
    independence_claim_gate = data.table::fifelse(
      partition_method[1L] == "drop_one_block_oof_delta_proxy",
      "empirical_sample_variance_not_physical_independence",
      "linear_layer_only_family_manifest_governed"
    )
  ), by = .(dataset_scope, granularity, cell_id, residue_scope, nucleus_id, nucleus_label,
            contributor_block, family_id, legal_permutation_block)]
  for (nm in c("block_R2", "unique_R2", "shared_R2", "delta_R2_drop_block", "permutation_null_p", "q_value")) {
    out[!is.finite(get(nm)), (nm) := NA_real_]
  }
  atlas_write_csv(out, atlas_table(cfg, "contributor_independence_partition.csv"))
  if (!is.null(rows) && !is.null(support)) {
    atlas_relaimpo_family_lmg(rows, support, cfg)
  }
  atlas_save_rds(out, cp)
  out
}

atlas_geometry_summaries <- function(rows, scores, cfg) {
  cp <- atlas_checkpoint(cfg, "geometry_summaries")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }
  atlas_msg("building dihedral and secondary-structure summaries")
  dt <- rows[supervised == TRUE]
  bin_angle <- function(x, bins = 36L) {
    cut(x, breaks = seq(-pi, pi, length.out = bins + 1L), include.lowest = TRUE, labels = FALSE)
  }
  dt[, `:=`(
    phi_bin36 = bin_angle(phi),
    psi_bin36 = bin_angle(psi),
    omega_bin36 = bin_angle(omega),
    chi1_angle = atan2(chi1_sin, chi1_cos),
    chi2_angle = atan2(chi2_sin, chi2_cos),
    chi3_angle = atan2(chi3_sin, chi3_cos),
    chi4_angle = atan2(chi4_sin, chi4_cos)
  )]
  for (ch in paste0("chi", 1:4, "_angle")) dt[, paste0(ch, "_bin36") := bin_angle(get(ch))]

  out <- list()
  k <- 1L
  for (view in c("within_1p9j", "between_720", "pooled")) {
    dtv <- atlas_view_subset(dt, view)
    for (granularity in c("per_residue_type", "pooled_over_residue")) {
      cell <- if (granularity == "per_residue_type") "residue_nucleus_id" else "pooled_nucleus_id"
      residue <- if (granularity == "per_residue_type") dtv$residue_scope else rep("ALL", nrow(dtv))
      tmp <- copy(dtv)
      tmp[, `:=`(
        dataset_scope_out = view, granularity_out = granularity,
        cell_id = get(cell), residue_scope_out = residue
      )]
      tmp[, dihedral_out := "phi_psi"]
      rama <- tmp[is.finite(phi) & is.finite(psi), .(
        n_rows = .N, sigma_mean = mean(sigma_iso), sigma_sd = stats::sd(sigma_iso),
        dft_floor_mse = mean(dft_floor_resid^2)
      ), by = .(dataset_scope = dataset_scope_out, granularity = granularity_out, dihedral = dihedral_out,
                cell_id, residue_scope = residue_scope_out, nucleus_id, nucleus_label,
                phi_bin36, psi_bin36, ss3, ss8)]
      data.table::setnames(rama, c("phi_bin36", "psi_bin36"), c("x_bin", "y_bin"))
      out[[k]] <- rama
      k <- k + 1L
      for (dihedral in c("omega", paste0("chi", 1:4))) {
        bcol <- paste0(dihedral, if (dihedral == "omega") "_bin36" else "_angle_bin36")
        exists_col <- paste0(dihedral, "_exists")
        if (dihedral == "omega") {
          ok <- is.finite(tmp[[bcol]])
        } else {
          ok <- tmp[[exists_col]] == 1 & is.finite(tmp[[bcol]])
        }
        tmp[, dihedral_out := dihedral]
        one <- tmp[ok, .(
          n_rows = .N, sigma_mean = mean(sigma_iso), sigma_sd = stats::sd(sigma_iso),
          dft_floor_mse = mean(dft_floor_resid^2)
        ), by = .(dataset_scope = dataset_scope_out, granularity = granularity_out, dihedral = dihedral_out,
                  cell_id, residue_scope = residue_scope_out, nucleus_id, nucleus_label,
                  x_bin = get(bcol), ss3, ss8, rotamer_label)]
        one[, y_bin := NA_integer_]
        out[[k]] <- one
        k <- k + 1L
      }
    }
  }
  geom <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  atlas_write_csv(geom, atlas_table(cfg, "dihedral_sigma_bins.csv"))

  ss <- dt[, .(
    n_rows = .N, n_proteins = data.table::uniqueN(protein_id),
    sigma_mean = mean(sigma_iso), sigma_sd = stats::sd(sigma_iso)
  ), by = .(dataset_id, residue_scope, nucleus_id, nucleus_label, ss3, ss8)]
  atlas_write_csv(ss, atlas_table(cfg, "secondary_structure_breakdowns.csv"))

  circ_one <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) {
      return(list(mean_direction_rad = NA_real_, rho = NA_real_, angular_sd = NA_real_))
    }
    cx <- circular::circular(x, units = "radians", modulo = "asis", template = "none")
    list(
      mean_direction_rad = as.numeric(circular::mean.circular(cx)),
      rho = as.numeric(circular::rho.circular(cx)),
      angular_sd = as.numeric(circular::sd.circular(cx))
    )
  }
  directional <- list()
  dk <- 1L
  angle_specs <- data.table::data.table(
    dihedral = c("phi", "psi", "omega", paste0("chi", 1:4)),
    angle_col = c("phi", "psi", "omega", paste0("chi", 1:4, "_angle")),
    exists_col = c(NA_character_, NA_character_, NA_character_, paste0("chi", 1:4, "_exists"))
  )
  for (view in c("within_1p9j", "between_720", "pooled")) {
    dtv <- atlas_view_subset(dt, view)
    for (granularity in c("per_residue_type", "pooled_over_residue")) {
      cell <- if (granularity == "per_residue_type") "residue_nucleus_id" else "pooled_nucleus_id"
      residue <- if (granularity == "per_residue_type") dtv$residue_scope else rep("ALL", nrow(dtv))
      tmp <- copy(dtv)
      tmp[, `:=`(
        dataset_scope_out = view, granularity_out = granularity,
        cell_id = get(cell), residue_scope_out = residue
      )]
      for (ai in seq_len(nrow(angle_specs))) {
        spec <- angle_specs[ai]
        ok <- is.finite(tmp[[spec$angle_col]])
        if (!is.na(spec$exists_col)) ok <- ok & tmp[[spec$exists_col]] == 1
        if (!any(ok)) next
        tmp[, dihedral_out := spec$dihedral]
        one <- tmp[ok, {
          cs <- circ_one(get(spec$angle_col))
          .(
            n_rows = .N,
            sigma_mean = mean(sigma_iso),
            sigma_sd = stats::sd(sigma_iso),
            mean_direction_rad = cs$mean_direction_rad,
            circular_rho = cs$rho,
            circular_sd = cs$angular_sd
          )
        }, by = .(
          dataset_scope = dataset_scope_out, granularity = granularity_out, dihedral = dihedral_out,
          cell_id, residue_scope = residue_scope_out, nucleus_id, nucleus_label,
          ss3, ss8, rama_region_label, rotamer_label
        )]
        one[, region_estimation_status := "fold_support_summary_from_emitted_region_labels"]
        directional[[dk]] <- one
        dk <- dk + 1L
      }
    }
  }
  directional <- data.table::rbindlist(directional, use.names = TRUE, fill = TRUE)
  atlas_write_csv(directional, atlas_table(cfg, "dihedral_directional_basin_support.csv"))

  atlas_save_rds(list(dihedral = geom, ss = ss, directional = directional), cp)
  list(dihedral = geom, ss = ss, directional = directional)
}

atlas_render_figures <- function(support, scores, geom, cfg) {
  cp <- atlas_checkpoint(cfg, "figure_manifest")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }
  if (!cfg$render_figures) {
    fm <- data.table::data.table()
    atlas_save_rds(fm, cp)
    return(fm)
  }
  atlas_msg("rendering atlas figures")
  fig_manifest <- list()
  i <- 1L
  theme_atlas <- ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      legend.position = "right"
    )

  heat <- scores[is.finite(R2_oof_floor), .(
    R2_oof_floor = max(R2_oof_floor, na.rm = TRUE),
    R_oof = R_oof[which.max(R2_oof_floor)][1L],
    support_gate = support_gate[which.max(R2_oof_floor)][1L]
  ), by = .(dataset_scope, granularity, nucleus_label, contributor_id, mechanism)]
  for (view in unique(heat$dataset_scope)) {
    for (granularity in unique(heat$granularity)) {
      this_granularity <- granularity
      dat <- heat[dataset_scope == view & granularity == this_granularity]
      if (!nrow(dat)) next
      top_contrib <- dat[, .(mx = max(R2_oof_floor, na.rm = TRUE)), by = contributor_id][order(-mx)][seq_len(min(.N, 60))]
      top_nuc <- dat[, .(mx = max(R2_oof_floor, na.rm = TRUE)), by = nucleus_label][order(-mx)][seq_len(min(.N, 80))]
      dat <- dat[contributor_id %in% top_contrib$contributor_id & nucleus_label %in% top_nuc$nucleus_label]
      path <- file.path(cfg$figure_dir, "contributor_matrix", paste0(view, "__", granularity, "__R2_matrix.png"))
      p <- ggplot2::ggplot(dat, ggplot2::aes(contributor_id, nucleus_label, fill = pmax(0, R2_oof_floor))) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(mechanism ~ ., scales = "free_y", space = "free_y") +
        viridis::scale_fill_viridis(option = "C", name = "OOF R2 floor") +
        ggplot2::labs(title = paste("Contributor signal matrix", view, granularity), x = NULL, y = NULL) +
        theme_atlas +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
      ragg::agg_png(path, width = 14, height = 10, units = "in", res = 180)
      print(p)
      grDevices::dev.off()
      fig_manifest[[i]] <- data.table::data.table(artifact_type = "figure", figure_kind = "contributor_matrix",
                                                   dataset_scope = view, granularity = granularity, path = path)
      i <- i + 1L
    }
  }

  if (cfg$render_cell_figures && nrow(geom$dihedral)) {
    candidates <- geom$dihedral[, .(n = sum(n_rows)), by = .(dataset_scope, granularity, dihedral, cell_id, nucleus_label)][n >= 20]
    for (j in seq_len(nrow(candidates))) {
      row <- candidates[j]
      dat <- geom$dihedral[
        dataset_scope == row$dataset_scope & granularity == row$granularity &
          dihedral == row$dihedral & cell_id == row$cell_id
      ]
      if (!nrow(dat)) next
      safe_cell <- gsub("[^A-Za-z0-9_.-]+", "_", row$cell_id)
      dir <- file.path(cfg$figure_dir, "dihedral", row$dataset_scope, row$granularity, row$dihedral)
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      path <- file.path(dir, paste0(safe_cell, ".png"))
      if (row$dihedral == "phi_psi" && all(c("x_bin", "y_bin") %in% names(dat))) {
        dat2 <- dat[!is.na(y_bin) & is.finite(sigma_mean)]
        if (data.table::uniqueN(dat2$x_bin) >= 4L && data.table::uniqueN(dat2$y_bin) >= 4L && nrow(dat2) >= 24L) {
          p <- ggplot2::ggplot(dat2, ggplot2::aes(x_bin, y_bin, z = sigma_mean)) +
            ggplot2::geom_contour_filled(bins = 10, na.rm = TRUE) +
            viridis::scale_fill_viridis(discrete = TRUE, option = "D", name = "mean sigma") +
            ggplot2::coord_equal() +
            ggplot2::labs(title = paste(row$nucleus_label, row$dataset_scope, row$granularity, "phi/psi"),
                          x = "phi bin", y = "psi bin") +
            theme_atlas
        } else {
          p <- ggplot2::ggplot(dat2, ggplot2::aes(x_bin, y_bin, fill = sigma_mean)) +
            ggplot2::geom_tile() +
            viridis::scale_fill_viridis(option = "D", name = "mean sigma") +
            ggplot2::coord_equal() +
            ggplot2::labs(title = paste(row$nucleus_label, row$dataset_scope, row$granularity, "phi/psi"),
                          x = "phi bin", y = "psi bin") +
            theme_atlas
        }
      } else {
        p <- ggplot2::ggplot(dat, ggplot2::aes(x_bin, sigma_mean, colour = ss3)) +
          ggplot2::geom_line(ggplot2::aes(group = ss3), linewidth = 0.4, na.rm = TRUE) +
          ggplot2::geom_point(size = 0.7, na.rm = TRUE) +
          viridis::scale_colour_viridis(discrete = TRUE, option = "E", name = "SS3") +
          ggplot2::labs(title = paste(row$nucleus_label, row$dataset_scope, row$granularity, row$dihedral),
                        x = paste(row$dihedral, "bin"), y = "mean sigma") +
          theme_atlas
      }
      ragg::agg_png(path, width = 6, height = 4, units = "in", res = 160)
      print(p)
      grDevices::dev.off()
      fig_manifest[[i]] <- data.table::data.table(
        artifact_type = "figure", figure_kind = "dihedral_sigma",
        dataset_scope = row$dataset_scope, granularity = row$granularity,
        dihedral = row$dihedral, cell_id = row$cell_id, nucleus_label = row$nucleus_label,
        path = path
      )
      i <- i + 1L
    }
  }

  fm <- data.table::rbindlist(fig_manifest, use.names = TRUE, fill = TRUE)
  atlas_write_csv(fm, atlas_table(cfg, "figure_manifest.csv"))
  atlas_save_rds(fm, cp)
  fm
}

atlas_mixed_and_gam_ladder <- function(rows, cfg) {
  cp <- atlas_checkpoint(cfg, "mixed_gam_variance_ladder")
  if (file.exists(cp)) {
    atlas_msg("loading checkpoint ", cp)
    return(readRDS(cp))
  }
  atlas_msg("fitting R mixed/GAM ladder summaries")
  dt <- rows[supervised == TRUE]
  dt[, element_factor := factor(element)]
  dt[, nucleus_factor := factor(pooled_nucleus_id)]
  dt[, protein_factor := factor(protein_id)]
  dt[, fold_factor := factor(fold_id)]
  dt[, dataset_factor := factor(dataset_id)]
  out <- list()
  k <- 1L
  for (view in c("within_1p9j", "between_720", "pooled")) {
    d <- atlas_view_subset(dt, view)
    if (nrow(d) < 100L) next
    d <- droplevels(d[is.finite(phi) & is.finite(psi)])
    fixed_vec <- character()
    if (view == "pooled" && nlevels(d$dataset_factor) > 1L) fixed_vec <- c(fixed_vec, "dataset_factor")
    if (nlevels(d$element_factor) > 1L) fixed_vec <- c(fixed_vec, "element_factor")
    if (!length(fixed_vec)) fixed_vec <- "1"
    nucleus_smooth <- if (nlevels(d$nucleus_factor) > 1L) "s(nucleus_factor, bs = 're')" else NULL
    unit_factor <- if (nlevels(d$protein_factor) > 1L) {
      "protein_factor"
    } else if (nlevels(d$fold_factor) > 1L) {
      "fold_factor"
    } else {
      NA_character_
    }
    unit_smooth <- if (!is.na(unit_factor)) paste0("s(", unit_factor, ", bs = 're')") else NULL
    rhs0 <- paste(c(fixed_vec, nucleus_smooth), collapse = " + ")
    rhs1 <- paste(c(fixed_vec, nucleus_smooth, unit_smooth), collapse = " + ")
    rhs2 <- paste(c(fixed_vec, nucleus_smooth, unit_smooth, "te(phi, psi, bs = c('cc', 'cc'), k = c(20, 20))"), collapse = " + ")
    f0 <- stats::as.formula(paste("sigma_iso ~", rhs0))
    f1 <- stats::as.formula(paste("sigma_iso ~", rhs1))
    f2 <- stats::as.formula(paste(
      "sigma_iso ~", rhs2
    ))
    fits <- list(
      L_identity = tryCatch(mgcv::bam(f0, data = d, discrete = TRUE, method = "fREML"), error = identity),
      L_protein = tryCatch(mgcv::bam(f1, data = d, discrete = TRUE, method = "fREML"), error = identity),
      L_dihedral = tryCatch(
        mgcv::bam(
          f2, data = d, discrete = TRUE, method = "fREML",
          knots = list(phi = c(-pi, pi), psi = c(-pi, pi))
        ),
        error = identity
      )
    )
    base_var <- stats::var(d$sigma_iso)
    for (nm in names(fits)) {
      fit <- fits[[nm]]
      if (inherits(fit, "error")) {
        out[[k]] <- data.table::data.table(dataset_scope = view, ladder_step = nm, n_rows = nrow(d),
                                           deviance_explained = NA_real_, adjusted_R2 = NA_real_,
                                           response_variance = base_var, model_engine = "mgcv::bam(discrete=TRUE)",
                                           fit_status = paste("error", conditionMessage(fit)))
      } else {
        pred <- as.numeric(stats::predict(fit, newdata = d))
        out[[k]] <- data.table::data.table(dataset_scope = view, ladder_step = nm, n_rows = nrow(d),
                                           deviance_explained = summary(fit)$dev.expl,
                                           adjusted_R2 = 1 - sum((d$sigma_iso - pred)^2) / max(sum((d$sigma_iso - mean(d$sigma_iso))^2), .Machine$double.eps),
                                           response_variance = base_var, model_engine = "mgcv::bam(discrete=TRUE)",
                                           fit_status = "ok")
      }
      k <- k + 1L
    }

    random_terms <- c(
      if (nlevels(d$nucleus_factor) > 1L) "(1|nucleus_factor)" else NULL,
      if (!is.na(unit_factor)) paste0("(1|", unit_factor, ")") else NULL
    )
    lmer_formula <- stats::as.formula(paste(
      "sigma_iso ~", paste(c(fixed_vec, random_terms), collapse = " + ")
    ))
    lmer_fit <- tryCatch(
      lme4::lmer(
        lmer_formula, data = d, REML = TRUE,
        control = lme4::lmerControl(calc.derivs = FALSE, check.nobs.vs.nRE = "ignore")
      ),
      error = identity
    )
    if (inherits(lmer_fit, "error")) {
      out[[k]] <- data.table::data.table(
        dataset_scope = view, ladder_step = "L_lme4_identity_unit",
        n_rows = nrow(d), model_engine = "lme4::lmer", response_variance = base_var,
        fixed_variance = NA_real_, random_variance = NA_real_, residual_variance = NA_real_,
        marginal_R2 = NA_real_, conditional_R2 = NA_real_,
        fit_status = paste("error", conditionMessage(lmer_fit))
      )
    } else {
      fixed_pred <- as.numeric(stats::model.matrix(lmer_fit) %*% lme4::fixef(lmer_fit))
      vc <- as.data.frame(lme4::VarCorr(lmer_fit))
      random_var <- sum(vc$vcov, na.rm = TRUE)
      resid_var <- sigma(lmer_fit)^2
      fixed_var <- stats::var(fixed_pred)
      denom <- fixed_var + random_var + resid_var
      out[[k]] <- data.table::data.table(
        dataset_scope = view, ladder_step = "L_lme4_identity_unit",
        n_rows = nrow(d), model_engine = "lme4::lmer", response_variance = base_var,
        fixed_variance = fixed_var, random_variance = random_var,
        residual_variance = resid_var,
        marginal_R2 = fixed_var / max(denom, .Machine$double.eps),
        conditional_R2 = (fixed_var + random_var) / max(denom, .Machine$double.eps),
        fit_status = "ok"
      )
    }
    k <- k + 1L
  }
  ladder <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  atlas_write_csv(ladder, atlas_table(cfg, "variance_decomposition_ladder.csv"))
  atlas_save_rds(ladder, cp)
  ladder
}

atlas_build_index <- function(cfg, support, scores, partition, t2, geom, figures, ladder) {
  atlas_msg("building master index")
  relaimpo_path <- atlas_table(cfg, "contributor_family_relaimpo_lmg.csv")
  relaimpo_rows <- if (file.exists(relaimpo_path)) nrow(data.table::fread(relaimpo_path, select = "dataset_scope")) else 0L
  tables <- data.table::data.table(
    artifact_type = "table",
    table_name = c(
      "analysis_cell_support", "contributor_score_by_nucleus",
      "contributor_independence_partition", "tensor_T2_equivariant_scores",
      "dihedral_sigma_bins", "secondary_structure_breakdowns",
      "dihedral_directional_basin_support",
      "variance_decomposition_ladder", "contributor_family_relaimpo_lmg"
    ),
    path = file.path(cfg$table_dir, c(
      "analysis_cell_support.csv", "contributor_score_by_nucleus.csv",
      "contributor_independence_partition.csv", "tensor_T2_equivariant_scores.csv",
      "dihedral_sigma_bins.csv", "secondary_structure_breakdowns.csv",
      "dihedral_directional_basin_support.csv",
      "variance_decomposition_ladder.csv", "contributor_family_relaimpo_lmg.csv"
    )),
    rows = c(nrow(support), nrow(scores), nrow(partition), nrow(t2),
             nrow(geom$dihedral), nrow(geom$ss), nrow(geom$directional),
             nrow(ladder), relaimpo_rows)
  )
  cell_manifest <- scores[, .(
    score_rows = .N,
    best_R2_oof_floor = suppressWarnings(max(R2_oof_floor, na.rm = TRUE)),
    best_abs_R = suppressWarnings(max(abs(R_oof), na.rm = TRUE)),
    min_null_p = suppressWarnings(min(permutation_null_p, na.rm = TRUE)),
    support_gate = support_gate[1L],
    n_rows = n_rows[1L],
    n_proteins = n_proteins[1L],
    n_eff = n_eff[1L],
    score_table = file.path(cfg$table_dir, "contributor_score_by_nucleus.csv")
  ), by = .(dataset_scope, granularity, cell_id, residue_scope, nucleus_id, nucleus_label)]
  for (nm in c("best_R2_oof_floor", "best_abs_R", "min_null_p")) {
    cell_manifest[!is.finite(get(nm)), (nm) := NA_real_]
  }
  atlas_write_csv(cell_manifest, atlas_table(cfg, "cell_manifest.csv"))
  index <- data.table::rbindlist(list(tables, figures), use.names = TRUE, fill = TRUE)
  atlas_write_csv(index, file.path(cfg$output_dir, "master_index.csv"))
  disk <- system2("du", c("-sh", cfg$output_dir), stdout = TRUE)
  writeLines(c(
    paste("input_dir:", cfg$input_dir),
    paste("output_dir:", cfg$output_dir),
    paste("built_at:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste("support_rows:", nrow(support)),
    paste("score_rows:", nrow(scores)),
    paste("partition_rows:", nrow(partition)),
    paste("t2_rows:", nrow(t2)),
    paste("figure_rows:", nrow(figures)),
    paste("output_du:", paste(disk, collapse = " "))
  ), con = file.path(cfg$output_dir, "ATLAS_BUILD_SUMMARY.txt"))
  index
}
