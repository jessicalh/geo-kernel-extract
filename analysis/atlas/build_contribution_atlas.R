#!/usr/bin/env Rscript

all_args <- commandArgs(FALSE)
file_arg <- all_args[grepl("^--file=", all_args)]
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  script_dir <- file.path(getwd(), "analysis", "atlas")
}
source(file.path(script_dir, "R", "atlas_lib.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
stage_arg <- args[grepl("^--stage=", args)]
stages <- if (length(stage_arg)) {
  strsplit(sub("^--stage=", "", stage_arg[[1]]), ",", fixed = TRUE)[[1]]
} else {
  c("support", "score", "t2", "partition", "geometry", "ladder", "figures", "index")
}
stage_requested <- function(stage, checkpoint = NULL, table = NULL) {
  stage %in% stages ||
    (!is.null(checkpoint) && file.exists(checkpoint)) ||
    (!is.null(table) && file.exists(table))
}

cfg <- atlas_default_config()
atlas_prepare_dirs(cfg)
sink(file.path(cfg$log_dir, paste0("build_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".log")), split = TRUE)
on.exit(sink(), add = TRUE)

atlas_msg("atlas build start")
atlas_msg("input: ", cfg$input_dir)
atlas_msg("output: ", cfg$output_dir)
atlas_require_packages()

rows <- atlas_read_rows(cfg)
support <- if ("support" %in% stages || file.exists(atlas_checkpoint(cfg, "analysis_cell_support"))) {
  atlas_make_support(rows, cfg)
} else {
  NULL
}

if (is.null(support)) stop("support stage is required before downstream stages", call. = FALSE)
if (identical(sort(stages), "support")) {
  atlas_msg("support-only stage complete")
  quit(save = "no", status = 0)
}

scores <- if (stage_requested("score", atlas_checkpoint(cfg, "contributor_score_by_nucleus"), atlas_table(cfg, "contributor_score_by_nucleus.csv"))) {
  atlas_score_contributors(rows, support, cfg)
} else {
  NULL
}

t2 <- if (stage_requested("t2", atlas_checkpoint(cfg, "tensor_T2_equivariant_scores"), atlas_table(cfg, "tensor_T2_equivariant_scores.csv"))) {
  atlas_score_t2(rows, support, cfg)
} else {
  NULL
}

partition <- if (stage_requested("partition", atlas_checkpoint(cfg, "contributor_independence_partition"), atlas_table(cfg, "contributor_independence_partition.csv"))) {
  if (is.null(scores)) scores <- data.table::fread(atlas_table(cfg, "contributor_score_by_nucleus.csv"))
  atlas_independence_partition(scores, rows, support, cfg)
} else {
  NULL
}

geom <- if (stage_requested("geometry", atlas_checkpoint(cfg, "geometry_summaries"), atlas_table(cfg, "dihedral_sigma_bins.csv"))) {
  atlas_geometry_summaries(rows, scores, cfg)
} else {
  NULL
}

ladder <- if (stage_requested("ladder", atlas_checkpoint(cfg, "mixed_gam_variance_ladder"), atlas_table(cfg, "variance_decomposition_ladder.csv"))) {
  atlas_mixed_and_gam_ladder(rows, cfg)
} else {
  NULL
}

figures <- if (stage_requested("figures", atlas_checkpoint(cfg, "figure_manifest"), atlas_table(cfg, "figure_manifest.csv"))) {
  if (is.null(scores)) scores <- data.table::fread(atlas_table(cfg, "contributor_score_by_nucleus.csv"))
  if (is.null(geom)) {
    geom <- list(
      dihedral = data.table::fread(atlas_table(cfg, "dihedral_sigma_bins.csv")),
      ss = data.table::fread(atlas_table(cfg, "secondary_structure_breakdowns.csv")),
      directional = data.table::fread(atlas_table(cfg, "dihedral_directional_basin_support.csv"))
    )
  }
  atlas_render_figures(support, scores, geom, cfg)
} else {
  NULL
}

if ("index" %in% stages) {
  if (is.null(scores)) scores <- data.table::fread(atlas_table(cfg, "contributor_score_by_nucleus.csv"))
  if (is.null(t2)) t2 <- data.table::fread(atlas_table(cfg, "tensor_T2_equivariant_scores.csv"))
  if (is.null(partition)) partition <- data.table::fread(atlas_table(cfg, "contributor_independence_partition.csv"))
  if (is.null(geom)) {
    geom <- list(
      dihedral = data.table::fread(atlas_table(cfg, "dihedral_sigma_bins.csv")),
      ss = data.table::fread(atlas_table(cfg, "secondary_structure_breakdowns.csv")),
      directional = data.table::fread(atlas_table(cfg, "dihedral_directional_basin_support.csv"))
    )
  }
  if (is.null(ladder)) ladder <- data.table::fread(atlas_table(cfg, "variance_decomposition_ladder.csv"))
  if (is.null(figures)) figures <- data.table::fread(atlas_table(cfg, "figure_manifest.csv"))
  atlas_build_index(cfg, support, scores, partition, t2, geom, figures, ladder)
}

atlas_msg("atlas build complete")
