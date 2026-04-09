#!/usr/bin/env Rscript
# scalar_ablation.R — Scalar group interaction analysis figures.
#
# Reads: output/secondary/scalar_ablation/*.csv
# Writes: output/secondary/figures/scalar_*.pdf
#
# Usage:
#   cd learn/R
#   Rscript scalar_ablation.R [data_dir]

args <- commandArgs(trailingOnly = TRUE)
source("common.R")
if (length(args) >= 1) set_data_dir(args[1])

sa_dir <- file.path(DATA_DIR, "scalar_ablation")

# ── Load ─────────────────────────────────────────────────────────────
interact <- read_csv(file.path(sa_dir, "interaction_r2.csv"), show_col_types = FALSE)
cumul    <- read_csv(file.path(sa_dir, "group_summary.csv"), show_col_types = FALSE)

eigen_file <- file.path(sa_dir, "eigenspectrum.csv")
has_eigen <- file.exists(eigen_file)
if (has_eigen) eigen <- read_csv(eigen_file, show_col_types = FALSE)

# ── 1. Per-group R² delta (horizontal bar chart) ────────────────────
interact <- interact[!is.na(interact$r2_delta), ]
interact$group_name <- factor(interact$group_name,
                               levels = interact$group_name[order(interact$r2_delta)])

p1 <- ggplot(interact, aes(x = group_name, y = r2_delta)) +
  geom_col(aes(fill = r2_delta > 0.005), alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey60"),
                    guide = "none") +
  geom_text(aes(label = sprintf("%+.4f", r2_delta)),
            hjust = ifelse(interact$r2_delta > 0, -0.1, 1.1), size = 2.5) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "grey30") +
  labs(title = "Scalar group interaction impact on ridge R\u00b2",
       subtitle = paste0("Base: top-10 kernels, R\u00b2 = ",
                         sprintf("%.4f", interact$r2_base[1])),
       x = NULL,
       y = expression(Delta * " R\u00b2 (with interactions - base)"))

save_fig(p1, "scalar_interaction_impact", width = 9, height = 7)

# ── 2. Cumulative R² as groups are added ─────────────────────────────
if (nrow(cumul) > 1) {
  cumul$step <- factor(cumul$step)

  p2 <- ggplot(cumul, aes(x = step, y = r2)) +
    geom_col(aes(fill = group), alpha = 0.8) +
    geom_text(aes(label = group), angle = 45, hjust = -0.1, size = 2.5) +
    geom_text(aes(label = sprintf("%.4f", r2)), vjust = -0.3, size = 2.5) +
    labs(title = "Cumulative ridge R\u00b2 as scalar groups are added",
         subtitle = "Groups added in order of individual impact",
         x = "Step", y = "Ridge R\u00b2") +
    guides(fill = "none") +
    coord_cartesian(ylim = c(0, max(cumul$r2, na.rm = TRUE) * 1.15))

  save_fig(p2, "scalar_cumulative_r2", width = 10, height = 5)
}

# ── 3. Eigenspectrum of interaction space ────────────────────────────
if (has_eigen) {
  p3 <- ggplot(eigen[eigen$rank <= 30, ],
               aes(x = rank, y = cumulative_variance)) +
    geom_line(linewidth = 0.8, colour = "#0072B2") +
    geom_point(size = 1.5, colour = "#0072B2") +
    geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey50") +
    annotate("text", x = 25, y = 0.91, label = "90%", colour = "grey40", size = 3) +
    annotate("text", x = 25, y = 0.96, label = "95%", colour = "grey40", size = 3) +
    labs(title = "Interaction feature space: eigenspectrum",
         subtitle = "How many effective dimensions when scalars modulate kernels?",
         x = "Eigenvalue rank", y = "Cumulative variance explained") +
    coord_cartesian(ylim = c(0, 1))

  save_fig(p3, "scalar_eigenspectrum", width = 8, height = 5)
}

# ── 4. Group width vs impact (scatter) ──────────────────────────────
p4 <- ggplot(interact, aes(x = group_width, y = r2_delta)) +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = interact[interact$r2_delta > 0.003 | interact$r2_delta < -0.003, ],
    aes(label = group_name), size = 2.5, max.overlaps = 20
  ) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  labs(title = "Scalar group size vs interaction impact",
       subtitle = "More features != more impact",
       x = "Number of scalars in group",
       y = expression(Delta * " R\u00b2"))

# ggrepel may not be installed — fall back gracefully
tryCatch(
  save_fig(p4, "scalar_width_vs_impact", width = 8, height = 6),
  error = function(e) {
    p4_fallback <- ggplot(interact, aes(x = group_width, y = r2_delta)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_text(data = interact[abs(interact$r2_delta) > 0.003, ],
                aes(label = group_name), size = 2.5, vjust = -1) +
      geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
      labs(title = "Scalar group size vs interaction impact",
           x = "Number of scalars in group",
           y = expression(Delta * " R\u00b2"))
    save_fig(p4_fallback, "scalar_width_vs_impact", width = 8, height = 6)
  }
)

message("Done: scalar ablation figures.")
