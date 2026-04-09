#!/usr/bin/env Rscript
# workbench.R — Analytical experimental grid on the kernel matrix.
#
# Replaces train-and-check-in with: one export, all configurations
# tested analytically via ridge regression.  No GPU, no epochs, no
# batch size.  Pure linear algebra on the kernel matrix.
#
# Produces:
#   figures/wb_kernel_grid.pdf        — ridge R² for every kernel group combo
#   figures/wb_forward_selection.pdf  — greedy forward selection curve
#   figures/wb_eigenspectrum.pdf      — effective dimensionality
#   figures/wb_interaction_grid.pdf   — scalar × kernel interaction R² grid
#   figures/wb_stratum_grid.pdf       — ridge R² by stratum for top configs
#
# Usage:
#   cd learn/R
#   Rscript workbench.R ../src/output/secondary/r_export

args <- commandArgs(trailingOnly = TRUE)
source("common.R")

if (length(args) < 1) {
  stop("Usage: Rscript workbench.R <export_dir>")
}
export_dir <- args[1]

source("load_export.R")
d <- load_export(export_dir)

# ── Helpers ──────────────────────────────────────────────────────────

ridge_r2 <- function(X, y, lambda = 0.01) {
  # X: (N, p), y: (N, 5)
  # Returns R² of ridge fit
  XtX <- crossprod(X) + lambda * diag(ncol(X))
  w <- solve(XtX, crossprod(X, y))
  pred <- X %*% w
  ss_res <- sum((y - pred)^2)
  ss_tot <- sum((y - matrix(colMeans(y), nrow(y), ncol(y), byrow = TRUE))^2)
  1 - ss_res / ss_tot
}

flatten_kernels <- function(K, cols = NULL) {
  # K: (N, n_k, 5) → (N, n_k * 5) or subset
  n <- dim(K)[1]
  if (is.null(cols)) {
    matrix(K, nrow = n)
  } else {
    matrix(K[, cols, , drop = FALSE], nrow = n)
  }
}

lambda <- d$cfg$ridge_lambda

cat("\n══════════════════════════════════════════════════════════\n")
cat("  ANALYTICAL WORKBENCH\n")
cat(sprintf("  %d atoms, %d kernels, %d scalars, λ=%.4f\n",
            d$cfg$n_atoms, d$cfg$n_kernels, d$cfg$n_scalars, lambda))
cat("══════════════════════════════════════════════════════════\n\n")

# ── 1. Kernel group combinations ─────────────────────────────────────
# Test every subset of kernel groups
groups <- d$cfg$kernel_groups
group_names <- names(groups)

X_full <- flatten_kernels(d$K_norm)
y <- d$target
r2_full <- ridge_r2(X_full, y, lambda)
cat(sprintf("Full ridge (%d kernels): R² = %.4f\n\n", d$cfg$n_kernels, r2_full))

# Single groups
cat("Single kernel groups:\n")
single_rows <- list()
for (gn in group_names) {
  rng <- groups[[gn]]
  cols <- seq(rng[1] + 1, rng[2])  # R is 1-indexed
  X_g <- flatten_kernels(d$K_norm, cols)
  r2 <- ridge_r2(X_g, y, lambda)
  n_k <- length(cols)
  cat(sprintf("  %-12s  %3d kernels  R² = %.4f\n", gn, n_k, r2))
  single_rows[[length(single_rows) + 1]] <- data.frame(
    group = gn, n_kernels = n_k, r2 = r2, stringsAsFactors = FALSE)
}
single_df <- do.call(rbind, single_rows)

# All pairs
cat("\nKernel group pairs:\n")
pair_rows <- list()
for (i in seq_along(group_names)) {
  for (j in seq(i + 1, length(group_names))) {
    if (j > length(group_names)) next
    g1 <- group_names[i]; g2 <- group_names[j]
    cols <- c(seq(groups[[g1]][1] + 1, groups[[g1]][2]),
              seq(groups[[g2]][1] + 1, groups[[g2]][2]))
    X_p <- flatten_kernels(d$K_norm, cols)
    r2 <- ridge_r2(X_p, y, lambda)
    cat(sprintf("  %-12s + %-12s  R² = %.4f\n", g1, g2, r2))
    pair_rows[[length(pair_rows) + 1]] <- data.frame(
      group1 = g1, group2 = g2, r2 = r2, stringsAsFactors = FALSE)
  }
}
pair_df <- do.call(rbind, pair_rows)

# Plot
p1 <- ggplot(single_df, aes(x = reorder(group, -r2), y = r2, fill = group)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", r2)), vjust = -0.3, size = 3) +
  labs(title = "Ridge R² by kernel group (single)",
       x = NULL, y = "Ridge R²") +
  guides(fill = "none") +
  coord_cartesian(ylim = c(0, max(r2_full * 1.1, 0.5)))
save_fig(p1, "wb_kernel_groups", width = 8, height = 5)

# ── 2. Greedy forward selection (per kernel) ─────────────────────────
cat("\nGreedy forward selection:\n")
n_k <- d$cfg$n_kernels
remaining <- seq_len(n_k)
selected <- integer(0)
fwd_rows <- list()

for (step in seq_len(min(n_k, 25))) {
  best_k <- -1; best_r2 <- -Inf
  for (k in remaining) {
    trial <- c(selected, k)
    X_t <- flatten_kernels(d$K_norm, trial)
    r2 <- tryCatch(ridge_r2(X_t, y, lambda), error = function(e) -Inf)
    if (r2 > best_r2) { best_r2 <- r2; best_k <- k }
  }
  if (best_k < 0) break
  selected <- c(selected, best_k)
  remaining <- setdiff(remaining, best_k)
  kname <- d$knames$name[best_k]
  delta <- if (length(fwd_rows) > 0) best_r2 - fwd_rows[[length(fwd_rows)]]$r2 else best_r2
  cat(sprintf("  %2d. +%-25s  R² = %.4f  (+%.4f)\n", step, kname, best_r2, delta))
  fwd_rows[[length(fwd_rows) + 1]] <- data.frame(
    step = step, kernel = kname, r2 = best_r2, delta = delta,
    stringsAsFactors = FALSE)
}
fwd_df <- do.call(rbind, fwd_rows)

p2 <- ggplot(fwd_df, aes(x = step, y = r2)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_text(aes(label = kernel), angle = 45, hjust = -0.1, size = 2,
            data = fwd_df[fwd_df$step <= 10, ]) +
  labs(title = "Greedy forward kernel selection",
       x = "Step", y = "Cumulative ridge R²") +
  coord_cartesian(ylim = c(0, r2_full * 1.1))
save_fig(p2, "wb_forward_selection", width = 10, height = 6)

# ── 3. Eigenspectrum ─────────────────────────────────────────────────
cat("\nEigenspectrum:\n")
K_flat <- matrix(d$K_norm, nrow = dim(d$K_norm)[1])
K_cent <- scale(K_flat, center = TRUE, scale = FALSE)
# Use crossproduct for efficiency: (p x p) instead of (N x N)
cov_mat <- crossprod(K_cent) / nrow(K_cent)
eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
eig <- pmax(eig, 0)
cumvar <- cumsum(eig) / sum(eig)
dim90 <- which(cumvar >= 0.90)[1]
dim95 <- which(cumvar >= 0.95)[1]
cat(sprintf("  Effective dim: %d (90%%), %d (95%%)\n", dim90, dim95))

eigen_df <- data.frame(rank = seq_along(eig), eigenvalue = eig,
                        cumvar = cumvar)[1:min(30, length(eig)), ]

p3 <- ggplot(eigen_df, aes(x = rank, y = cumvar)) +
  geom_line(linewidth = 0.8, colour = "#0072B2") +
  geom_point(size = 1.5, colour = "#0072B2") +
  geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey50") +
  labs(title = sprintf("Kernel eigenspectrum (%d kernels, %d atoms)",
                       d$cfg$n_kernels, d$cfg$n_atoms),
       subtitle = sprintf("Effective dim: %d (90%%), %d (95%%)", dim90, dim95),
       x = "Eigenvalue rank", y = "Cumulative variance") +
  coord_cartesian(ylim = c(0, 1))
save_fig(p3, "wb_eigenspectrum", width = 8, height = 5)

# ── 4. Scalar interaction grid ───────────────────────────────────────
# Use top-10 from forward selection as base
cat("\nScalar interaction analysis (top-10 kernel base):\n")
top10 <- selected[1:min(10, length(selected))]
X_base <- flatten_kernels(d$K_norm, top10)
r2_base <- ridge_r2(X_base, y, lambda)
cat(sprintf("  Base (top-10): R² = %.4f\n", r2_base))

interact_rows <- list()
for (i in seq_len(nrow(d$slayout))) {
  sname <- d$slayout$name[i]
  soff  <- d$slayout$offset[i] + 1  # R 1-indexed
  swid  <- d$slayout$width[i]
  scols <- seq(soff, soff + swid - 1)

  S_group <- d$scalars[, scols, drop = FALSE]

  # Build interaction: each scalar × each base feature
  interactions <- do.call(cbind, lapply(seq_len(swid), function(si) {
    S_group[, si] * X_base
  }))
  X_int <- cbind(X_base, interactions)

  n_feat <- ncol(X_int)
  if (nrow(X_int) < n_feat * 2) {
    r2_int <- NA
  } else {
    r2_int <- tryCatch(ridge_r2(X_int, y, lambda), error = function(e) NA)
  }
  delta <- if (!is.na(r2_int)) r2_int - r2_base else NA

  cat(sprintf("  %-25s  width=%3d  R² = %.4f  Δ = %+.4f\n",
              sname, swid,
              ifelse(is.na(r2_int), 0, r2_int),
              ifelse(is.na(delta), 0, delta)))

  interact_rows[[length(interact_rows) + 1]] <- data.frame(
    group = sname, width = swid, r2 = r2_int, delta = delta,
    stringsAsFactors = FALSE)
}
interact_df <- do.call(rbind, interact_rows)
interact_df <- interact_df[!is.na(interact_df$delta), ]
interact_df <- interact_df[order(-interact_df$delta), ]

p4 <- ggplot(interact_df,
             aes(x = reorder(group, delta), y = delta,
                 fill = delta > 0.005)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey60"),
                    guide = "none") +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "grey30") +
  labs(title = "Scalar group interaction impact (R analytical)",
       subtitle = sprintf("Base: top-10 kernels, R² = %.4f", r2_base),
       x = NULL, y = expression(Delta * " R²"))
save_fig(p4, "wb_interaction_grid", width = 9, height = 7)

# ── 5. Per-stratum ridge ─────────────────────────────────────────────
cat("\nPer-stratum ridge:\n")
strata <- c("hie_only", "phe_only", "tyr_only", "trp_only", "no_hie", "all")

# Build stratum masks from ring_type_mask
assign_stratum <- function(mask) {
  if (bitwAnd(mask, 128) != 0 && bitwAnd(mask, bitwNot(128)) == 0) return("hie_only")
  if (bitwAnd(mask, 1) != 0 && bitwAnd(mask, bitwNot(1)) == 0) return("phe_only")
  if (bitwAnd(mask, 2) != 0 && bitwAnd(mask, bitwNot(2)) == 0) return("tyr_only")
  trp_bits <- bitwOr(bitwOr(4, 8), 16)
  if (bitwAnd(mask, trp_bits) != 0 && bitwAnd(mask, bitwNot(trp_bits)) == 0) return("trp_only")
  if (bitwAnd(mask, 128) == 0 && mask != 0) return("no_hie")
  return("all")
}

d$meta$stratum <- sapply(d$meta$ring_type_mask, assign_stratum)

strat_rows <- list()
for (s in strata) {
  if (s == "all") {
    mask <- rep(TRUE, nrow(d$meta))
  } else {
    mask <- d$meta$stratum == s
  }
  n <- sum(mask)
  if (n < 100) {
    cat(sprintf("  %-12s  n=%6d  (too few)\n", s, n))
    next
  }
  X_s <- flatten_kernels(d$K_norm)[mask, , drop = FALSE]
  y_s <- y[mask, , drop = FALSE]
  r2 <- tryCatch(ridge_r2(X_s, y_s, lambda), error = function(e) NA)
  cat(sprintf("  %-12s  n=%6d  R² = %.4f\n", s, n, ifelse(is.na(r2), 0, r2)))
  strat_rows[[length(strat_rows) + 1]] <- data.frame(
    stratum = s, n_atoms = n, r2 = r2, stringsAsFactors = FALSE)
}
strat_df <- do.call(rbind, strat_rows)

p5 <- ggplot(strat_df, aes(x = stratum, y = r2, fill = stratum)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f\nn=%d", r2, n_atoms)),
            vjust = -0.2, size = 2.5) +
  stratum_fill_scale() +
  labs(title = "Ridge R² by ring-type stratum (R analytical)",
       x = NULL, y = "Ridge R²") +
  guides(fill = "none") +
  coord_cartesian(ylim = c(0, max(strat_df$r2, na.rm = TRUE) * 1.2))
save_fig(p5, "wb_stratum_grid", width = 8, height = 5)

cat("\n══════════════════════════════════════════════════════════\n")
cat("  Done. Figures in output/secondary/figures/wb_*.pdf\n")
cat("══════════════════════════════════════════════════════════\n")
