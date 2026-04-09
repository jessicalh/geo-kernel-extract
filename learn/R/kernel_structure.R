#!/usr/bin/env Rscript
# kernel_structure.R — Inter-kernel T2 structure figures.
#
# Reads: output/secondary/kernel_structure/*.csv
# Writes: output/secondary/figures/kernel_*.pdf
#
# Usage:
#   cd learn/R
#   Rscript kernel_structure.R [data_dir]

args <- commandArgs(trailingOnly = TRUE)
source("common.R")
if (length(args) >= 1) set_data_dir(args[1])

ks_dir <- file.path(DATA_DIR, "kernel_structure")

# ── Load ─────────────────────────────────────────────────────────────
cosmat <- read_csv(file.path(ks_dir, "cosine_matrix.csv"), show_col_types = FALSE)
eigen  <- read_csv(file.path(ks_dir, "eigenspectrum.csv"), show_col_types = FALSE)

# ── 1. Heatmap: kernel cosine similarity ("all" stratum) ────────────
cos_all <- cosmat[cosmat$stratum == "all", ]

p1 <- ggplot(cos_all, aes(x = kernel_i, y = kernel_j, fill = cosine_sim)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Cosine\nsimilarity") +
  labs(title = "Pairwise T2 cosine similarity (91 kernels)",
       subtitle = "All matched atoms",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        axis.text.y = element_text(size = 4))

save_fig(p1, "kernel_cosine_heatmap_all", width = 14, height = 12)

# ── 2. Eigenspectrum by stratum ──────────────────────────────────────
p2 <- ggplot(eigen[eigen$rank <= 20, ],
             aes(x = rank, y = cumulative_variance, colour = stratum)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey50") +
  annotate("text", x = 18, y = 0.91, label = "90%", colour = "grey40", size = 3) +
  annotate("text", x = 18, y = 0.96, label = "95%", colour = "grey40", size = 3) +
  stratum_fill_scale(aesthetics = "colour") +
  labs(title = "Eigenspectrum of kernel covariance by stratum",
       subtitle = "How many independent directions in the 91-kernel space?",
       x = "Eigenvalue rank", y = "Cumulative variance explained") +
  coord_cartesian(ylim = c(0, 1))

save_fig(p2, "kernel_eigenspectrum", width = 8, height = 5)

# ── 3. Side-by-side: HIE vs PHE heatmaps ────────────────────────────
cos_hie <- cosmat[cosmat$stratum == "hie_only", ]
cos_phe <- cosmat[cosmat$stratum == "phe_only", ]

if (nrow(cos_hie) > 0 && nrow(cos_phe) > 0) {
  both <- rbind(
    transform(cos_hie, panel = "HIE-only atoms"),
    transform(cos_phe, panel = "PHE-only atoms")
  )

  p3 <- ggplot(both, aes(x = kernel_i, y = kernel_j, fill = cosine_sim)) +
    geom_tile() +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Cosine\nsim") +
    facet_wrap(~panel, nrow = 1) +
    labs(title = "Kernel space structure: HIE vs PHE atoms",
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
          axis.text.y = element_text(size = 3))

  save_fig(p3, "kernel_hie_vs_phe", width = 20, height = 10)
}

message("Done: kernel structure figures.")
