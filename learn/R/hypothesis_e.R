#!/usr/bin/env Rscript
# hypothesis_e.R — Decompose ring_proximity to find the dominant coordinate.
#
# Tests which geometric coordinates within ring_proximity drive the
# +0.305 interaction delta.  Answers: is the path forward radial
# shaping, angular shaping, or azimuthal structure?
#
# Usage:
#   cd learn/R
#   Rscript hypothesis_e.R ../src/output/secondary/r_export

args <- commandArgs(trailingOnly = TRUE)
source("common.R")
if (length(args) < 1) stop("Usage: Rscript hypothesis_e.R <export_dir>")

source("load_export.R")
d <- load_export(args[1])

# ── Setup ────────────────────────────────────────────────────────────

ridge_r2 <- function(X, y, lambda = 0.01) {
  XtX <- crossprod(X) + lambda * diag(ncol(X))
  w <- solve(XtX, crossprod(X, y))
  pred <- X %*% w
  ss_res <- sum((y - pred)^2)
  ss_tot <- sum((y - matrix(colMeans(y), nrow(y), ncol(y), byrow = TRUE))^2)
  1 - ss_res / ss_tot
}

flatten_kernels <- function(K, cols = NULL) {
  n <- dim(K)[1]
  if (is.null(cols)) matrix(K, nrow = n)
  else matrix(K[, cols, , drop = FALSE], nrow = n)
}

lambda <- d$cfg$ridge_lambda
y <- d$target

# Forward-select top 10 (reuse from workbench logic)
n_k <- d$cfg$n_kernels
remaining <- seq_len(n_k)
selected <- integer(0)
for (step in seq_len(10)) {
  best_k <- -1; best_r2 <- -Inf
  for (k in remaining) {
    trial <- c(selected, k)
    X_t <- flatten_kernels(d$K_norm, trial)
    r2 <- tryCatch(ridge_r2(X_t, y, lambda), error = function(e) -Inf)
    if (r2 > best_r2) { best_r2 <- r2; best_k <- k }
  }
  selected <- c(selected, best_k)
  remaining <- setdiff(remaining, best_k)
}
X_base <- flatten_kernels(d$K_norm, selected)
r2_base <- ridge_r2(X_base, y, lambda)

cat(sprintf("\nBase (top-10 kernels): R² = %.4f\n", r2_base))
cat(sprintf("Top 10: %s\n\n", paste(d$knames$name[selected], collapse = ", ")))

# ── Ring proximity layout ────────────────────────────────────────────
# Per ring i (i=0..5), 10 columns each:
#   +0 mcconnell_factor, +1 exp_decay, +2 1/dist,
#   +3 z, +4 rho, +5 theta,
#   +6 disp_scalar, +7 disp_contacts,
#   +8 cos_phi, +9 sin_phi
# Column 60: n_pairs count

rp_offset <- 57 + 1  # R is 1-indexed, ring_proximity starts at scalar col 57
K_RINGS <- 6
COLS_PER_RING <- 13

# Helper: extract specific per-ring columns across all 6 rings
get_rp_cols <- function(within_ring_offsets) {
  cols <- integer(0)
  for (ri in 0:(K_RINGS - 1)) {
    for (off in within_ring_offsets) {
      cols <- c(cols, rp_offset + ri * COLS_PER_RING + off)
    }
  }
  cols
}

# ── Define subsets ───────────────────────────────────────────────────

subsets <- list(
  "1/dist only"            = get_rp_cols(2),
  "z + rho"                = get_rp_cols(c(3, 4)),
  "theta only"             = get_rp_cols(5),
  "z + rho + theta"        = get_rp_cols(c(3, 4, 5)),
  "cos_phi + sin_phi"      = get_rp_cols(c(8, 9)),
  "mcconnell_factor"       = get_rp_cols(0),
  "exp_decay"              = get_rp_cols(1),
  "disp_scalar + contacts" = get_rp_cols(c(6, 7)),
  "bs_T0 (signed)"         = get_rp_cols(10),
  "hm_T0 (signed)"         = get_rp_cols(11),
  "chi_T0 (signed)"        = get_rp_cols(12),
  "bs+hm+chi T0 (signed)"  = get_rp_cols(c(10, 11, 12)),
  "1/dist + z + rho"       = get_rp_cols(c(2, 3, 4)),
  "cylindrical (z,rho,theta,1/d)" = get_rp_cols(c(2, 3, 4, 5)),
  "cylindrical + phi"      = get_rp_cols(c(2, 3, 4, 5, 8, 9)),
  "signed T0 + phi"        = get_rp_cols(c(8, 9, 10, 11, 12)),
  "signed T0 + cylindrical" = get_rp_cols(c(2, 3, 4, 5, 10, 11, 12)),
  "signed T0 + cyl + phi"  = get_rp_cols(c(2, 3, 4, 5, 8, 9, 10, 11, 12)),
  "all geometry (no disp)"  = get_rp_cols(c(0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12)),
  "full ring_proximity"    = seq(rp_offset, rp_offset + 78)
)

# ── Test each subset ─────────────────────────────────────────────────

cat("══════════════════════════════════════════════════════════\n")
cat("  HYPOTHESIS E: Ring proximity decomposition\n")
cat("══════════════════════════════════════════════════════════\n\n")

results <- list()
for (name in names(subsets)) {
  scols <- subsets[[name]]
  S_sub <- d$scalars[, scols, drop = FALSE]
  n_scalars <- ncol(S_sub)

  # Build interaction: each scalar × each base kernel feature
  interactions <- do.call(cbind, lapply(seq_len(n_scalars), function(si) {
    S_sub[, si] * X_base
  }))
  X_int <- cbind(X_base, interactions)

  r2_int <- tryCatch(ridge_r2(X_int, y, lambda), error = function(e) NA)
  delta <- if (!is.na(r2_int)) r2_int - r2_base else NA

  cat(sprintf("  %-35s  %2d scalars  R² = %.4f  Δ = %+.4f\n",
              name, n_scalars,
              ifelse(is.na(r2_int), 0, r2_int),
              ifelse(is.na(delta), 0, delta)))

  results[[length(results) + 1]] <- data.frame(
    subset = name, n_scalars = n_scalars,
    r2 = r2_int, delta = delta, stringsAsFactors = FALSE)
}

res_df <- do.call(rbind, results)
res_df <- res_df[order(-res_df$delta), ]

# ── Plot ─────────────────────────────────────────────────────────────

p <- ggplot(res_df, aes(x = reorder(subset, delta), y = delta,
                         fill = delta > 0.1)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#56B4E9"),
                    guide = "none") +
  geom_text(aes(label = sprintf("%+.3f (%d)", delta, n_scalars)),
            hjust = -0.05, size = 2.5) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "grey30") +
  labs(title = "Hypothesis E: Which ring proximity coordinates matter?",
       subtitle = sprintf("Base: top-10 kernels, R² = %.4f.  Interaction delta shown.",
                          r2_base),
       x = NULL,
       y = expression(Delta * " R² (interaction with top-10 kernels)")) +
  theme(plot.title = element_text(size = 11))

save_fig(p, "wb_hypothesis_e", width = 10, height = 7)

cat("\n══════════════════════════════════════════════════════════\n")
cat("  Done. Figure: wb_hypothesis_e.pdf\n")
cat("══════════════════════════════════════════════════════════\n")
