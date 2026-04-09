#!/usr/bin/env Rscript
# ablation.R — Per-ring-type R^2 ablation figures.
#
# Reads: output/secondary/ablation/*.csv
# Writes: output/secondary/figures/ablation_*.pdf
#
# Usage:
#   cd learn/R
#   Rscript ablation.R [data_dir]

args <- commandArgs(trailingOnly = TRUE)
source("common.R")
if (length(args) >= 1) set_data_dir(args[1])

abl_dir <- file.path(DATA_DIR, "ablation")

# ── Load ─────────────────────────────────────────────────────────────
r2_strat  <- read_csv(file.path(abl_dir, "ring_type_r2.csv"), show_col_types = FALSE)
self_fit  <- read_csv(file.path(abl_dir, "self_fit_r2.csv"), show_col_types = FALSE)
counts    <- read_csv(file.path(abl_dir, "atom_counts.csv"), show_col_types = FALSE)

# ── 1. Per-stratum R² by kernel family (grouped bar) ────────────────
r2_long <- tidyr::pivot_longer(
  r2_strat, cols = starts_with("r2_"),
  names_to = "family", values_to = "r2",
  names_prefix = "r2_"
)
r2_long <- r2_long[!is.na(r2_long$r2), ]

p1 <- ggplot(r2_long, aes(x = stratum, y = r2, fill = family)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f", r2)),
            position = position_dodge(0.9), vjust = -0.3, size = 2.2) +
  labs(title = "Ridge R\u00b2 by stratum and kernel family",
       subtitle = paste("n =", paste(r2_strat$stratum, r2_strat$n_atoms,
                                      sep = ":", collapse = "  ")),
       x = NULL, y = "Ridge R\u00b2", fill = "Family") +
  coord_cartesian(ylim = c(min(0, min(r2_long$r2, na.rm = TRUE) - 0.05), 1))

save_fig(p1, "ablation_r2_by_stratum", width = 11, height = 6)

# ── 2. Self-fit R²: the money plot ──────────────────────────────────
self_fit <- self_fit[self_fit$n_atoms > 0, ]

p2 <- ggplot(self_fit, aes(x = reorder(ring_type, -r2_self_4kernel),
                            y = r2_self_4kernel, fill = ring_type)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", r2_self_4kernel)),
            vjust = -0.3, size = 3) +
  ring_fill_scale() +
  labs(title = "Self-fit R\u00b2: each ring type's own 4 kernels",
       subtitle = "How well does each ring type's geometry explain its own DFT delta?",
       x = NULL, y = "Ridge R\u00b2 (4 kernels)") +
  guides(fill = "none")

save_fig(p2, "ablation_self_fit", width = 8, height = 5)

# ── 3. Self-fit vs all-kernels R² (paired comparison) ───────────────
p3 <- ggplot(self_fit, aes(x = reorder(ring_type, -r2_all_kernels))) +
  geom_col(aes(y = r2_all_kernels, fill = "All 91 kernels"),
           alpha = 0.5, width = 0.7) +
  geom_col(aes(y = r2_self_4kernel, fill = "Self (4 kernels)"),
           alpha = 0.9, width = 0.4) +
  scale_fill_manual(values = c("All 91 kernels" = "grey60",
                                "Self (4 kernels)" = "#D55E00")) +
  labs(title = "Self-fit vs full kernel set by ring type",
       subtitle = "Gap = information from other calculators / other ring types",
       x = NULL, y = "Ridge R\u00b2", fill = NULL)

save_fig(p3, "ablation_self_vs_all", width = 9, height = 5)

# ── 4. Target magnitude by ring type ────────────────────────────────
p4 <- ggplot(self_fit, aes(x = reorder(ring_type, -target_mag_mean),
                            y = target_mag_mean, fill = ring_type)) +
  geom_col(alpha = 0.8) +
  ring_fill_scale() +
  labs(title = "DFT delta T2 magnitude by ring type",
       subtitle = "Larger = stronger aromatic signature in DFT",
       x = NULL, y = "Mean |T2| (ppm)") +
  guides(fill = "none")

save_fig(p4, "ablation_target_magnitude", width = 8, height = 5)

message("Done: 4 ablation figures.")
