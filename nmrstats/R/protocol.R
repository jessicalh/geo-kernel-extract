read_protocol <- function(path = system.file("protocol", "primary_claims.json", package = "nmrstats")) {
  jsonlite::read_json(path, simplifyVector = TRUE)
}

check_protocol <- function(results) {
  deviations <- data.table::data.table(
    item = character(),
    status = character(),
    detail = character()
  )
  add <- function(item, status, detail) {
    deviations <<- rbind(deviations, data.table::data.table(item = item, status = status, detail = detail))
  }
  if (!is.null(results$adjustments)) {
    for (i in seq_along(results$adjustments)) add(names(results$adjustments)[i], "adjusted", results$adjustments[[i]])
  }
  if (is.null(results$ladder$average) || !nrow(results$ladder$average)) add("ladder", "missing", "no ladder average table")
  if (is.null(results$importance$drop_one_block$average)) add("importance", "missing", "no drop-one-block table")
  if (is.null(results$nulls)) add("nulls", "missing", "null controls absent")
  if (!nrow(deviations)) add("protocol", "ok", "no deviations recorded")
  deviations
}

save_results <- function(out, results) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_dt(results$support_per_stratum, file.path(out, "support_per_stratum.csv"))
  write_dt(results$folds, file.path(out, "folds.csv"))
  write_dt(results$region_summary, file.path(out, "region_summary.csv"))
  write_dt(results$ladder$average, file.path(out, "ladder_average.csv"))
  write_dt(results$ladder$by_fold, file.path(out, "ladder_by_fold.csv"))
  write_dt(results$ladder$timings, file.path(out, "ladder_timings.csv"))
  write_dt(results$joint_views, file.path(out, "joint_views.csv"))
  write_dt(results$element_fits, file.path(out, "element_gaulss_fits.csv"))
  write_dt(results$importance$relaimpo, file.path(out, "importance_relaimpo.csv"))
  write_dt(results$importance$drop_one_block$average, file.path(out, "importance_drop_one_block.csv"))
  write_dt(results$importance$drop_one_block$by_fold, file.path(out, "importance_drop_one_block_by_fold.csv"))
  write_dt(results$importance$correspondence, file.path(out, "importance_correspondence.csv"))
  write_dt(null_average_table(results$nulls), file.path(out, "nulls_average.csv"))
  write_dt(results$calibration$coverage, file.path(out, "calibration_coverage.csv"))
  write_dt(results$calibration$scores, file.path(out, "calibration_scores.csv"))
  write_dt(results$spatial_moran, file.path(out, "spatial_moran_sequence_proxy.csv"))
  write_dt(results$protocol, file.path(out, "protocol_deviations.csv"))
  saveRDS(results, file.path(out, "nmrstats_results.rds"))
  invisible(out)
}
