#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

required <- c(
  "arrow", "data.table", "jsonlite", "RcppCNPy", "mgcv", "lme4",
  "relaimpo", "circular", "ggplot2", "patchwork", "ggrepel", "viridis",
  "ragg", "dplyr"
)

optional <- c(
  "ggtext", "ggdist", "ggridges", "ggnewscale", "scico", "MetBrewer",
  "sandwich", "clubSandwich", "fst", "qs", "metR"
)

pkgs <- unique(c(required, optional))
installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)
if (length(missing)) {
  install.packages(missing, Ncpus = max(1L, parallel::detectCores(logical = FALSE) - 1L))
}

still_missing <- setdiff(required, rownames(installed.packages()))
if (length(still_missing)) {
  stop("Required R packages still missing: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

optional_missing <- setdiff(optional, rownames(installed.packages()))
if (length(optional_missing)) {
  message("Optional graphics/SE packages not installed: ", paste(optional_missing, collapse = ", "))
}

message("R atlas package check complete.")
