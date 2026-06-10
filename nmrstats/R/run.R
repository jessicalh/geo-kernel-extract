audit_adjustments <- function(rows, nfold) {
  blocks <- feature_blocks(rows)
  excluded_efg <- attr(blocks, "excluded_bare_efg_kernel") %||% character()
  atom_uids <- data.table::uniqueN(paste(rows$protein_id, rows$atom_index, sep = ":"))
  adj <- list(
    manifest_row_count = "real schema_audit.json uses row_count, not expected_rows; loader accepts either and validates the real count",
    support_table = "support_table.csv is global; per-stratum support is computed from emitted rows in R",
    grouped_cv = sprintf("static leave-protein-out is grouped into %d full-scope folds; all rows are tested once, but this is not 720 separate LOPO fits", nfold),
    protein_atom_re = sprintf("s(protein_atom) omitted from bam formulas: %d protein:atom levels, most static levels are singletons; protein_id/residue/equivalence REs retained", atom_uids),
    gaulss_bam = "installed mgcv::bam reports 'general families not supported by bam' for gaulss; per-element heteroscedastic fits use full-scope two-stage bam mean + log-residual-scale models",
    region_round_trip = "region_def_id is the emitted row_design_v1 placeholder; regions are estimated and written, but not used as varying row tags until the C++ re-emit follow-on",
    spatial_moran = "row_design has no coordinate columns; reported Moran's I is a within-protein sequence-adjacency proxy",
    sign_controls = "column_irrep_schema declares no sign_flip_legal proxy columns, so sign-control is reported as non-applicable"
  )
  if (length(excluded_efg)) {
    adj$bare_efg_units <- sprintf(
      "excluded bare EFG kernel columns from ppm target models by unit gate: %s",
      paste(excluded_efg, collapse = ", ")
    )
  }
  if (!"mc" %in% names(blocks)) adj$mc_block <- "mc_lit_T0/mc_lit_absT2 are all missing in the real rows, so the McConnell block is omitted"
  adj
}

run_all <- function(dir, out, nfold = 3L, nthreads = 1L,
                    load_t2 = TRUE, save_fits = FALSE) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rows <- load_substrate(dir, load_t2 = load_t2, load_embedding = FALSE)
  if (!nrow(rows)) abort("no rows loaded")
  message(sprintf("run_all: loaded %d rows", nrow(rows)))

  support <- support_per_stratum(rows)
  folds <- make_folds(rows, k = nfold)
  folds_dt <- fold_summary(rows, folds)
  regions <- with_fold_each_regions(rows, folds)
  freeze_regions(regions, file.path(out, "region_spec.json"))

  ladder <- run_ladder(rows, folds, nthreads = nthreads)
  element_fits <- fit_per_element(rows, nthreads = nthreads, save_dir = if (save_fits) file.path(out, "fits") else NULL)
  joint <- joint_views(ladder$predictions)
  imp <- importance(rows, folds, predictions = ladder$predictions, nthreads = nthreads)
  nulls <- run_nulls(rows, folds, nthreads = nthreads)
  cal <- calibration_coverage(ladder$predictions)
  spatial <- within_protein_morans_i(ladder$predictions, rows)

  results <- list(
    dirs = normalizePath(dir, mustWork = TRUE),
    out = normalizePath(out, mustWork = FALSE),
    n_rows = nrow(rows),
    nfold = nfold,
    nthreads = nthreads,
    support_per_stratum = support,
    folds = folds_dt,
    regions = regions,
    region_summary = region_summary(regions),
    ladder = ladder,
    element_fits = element_fits,
    joint_views = joint,
    importance = imp,
    nulls = nulls,
    calibration = cal,
    spatial_moran = spatial,
    adjustments = audit_adjustments(rows, nfold)
  )
  results$protocol <- check_protocol(results)
  save_results(out, results)
  invisible(results)
}
