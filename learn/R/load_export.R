# load_export.R — Load the Python-exported kernel matrix into R.
#
# Loads NPZ files via reticulate, returns a list with all data ready
# for analysis.  Source this once, then use the returned object.
#
# Usage:
#   source("load_export.R")
#   d <- load_export("../src/output/secondary/r_export")
#   # d$K_raw    — (N, K, 5) raw kernels
#   # d$K_norm   — (N, K, 5) per-protein normalized
#   # d$target   — (N, 5) DFT delta T2
#   # d$scalars  — (N, S) scalar features
#   # d$meta     — data.frame: protein_id, element, ring_dist, ...
#   # d$knames   — data.frame: kernel index, name, group
#   # d$slayout  — data.frame: scalar block name, offset, width
#   # d$cfg      — list: n_atoms, n_kernels, ridge_lambda, groups, ...

library(reticulate)
np <- import("numpy")

load_export <- function(export_dir) {
  dir <- normalizePath(export_dir, mustWork = TRUE)

  load_npz <- function(name) {
    f <- np$load(file.path(dir, name), allow_pickle = TRUE)
    f$f[["data"]]
  }

  K_raw  <- load_npz("kernels_raw.npz")
  K_norm <- load_npz("kernels_norm.npz")
  target <- load_npz("target.npz")
  scalars <- load_npz("scalars.npz")
  kscales <- load_npz("kernel_scales.npz")

  meta   <- read.csv(file.path(dir, "metadata.csv"), stringsAsFactors = FALSE)
  pmeta  <- read.csv(file.path(dir, "protein_meta.csv"), stringsAsFactors = FALSE)
  knames <- read.csv(file.path(dir, "kernel_names.csv"), stringsAsFactors = FALSE)
  slayout <- read.csv(file.path(dir, "scalar_layout.csv"), stringsAsFactors = FALSE)
  cfg    <- jsonlite::fromJSON(file.path(dir, "config.json"))

  cat(sprintf("Loaded: %d atoms, %d proteins, %d kernels, %d scalars\n",
              nrow(meta), nrow(pmeta), ncol(K_raw[1,,]), ncol(scalars)))

  list(
    K_raw = K_raw,
    K_norm = K_norm,
    target = target,
    scalars = scalars,
    kscales = kscales,
    meta = meta,
    pmeta = pmeta,
    knames = knames,
    slayout = slayout,
    cfg = cfg
  )
}
