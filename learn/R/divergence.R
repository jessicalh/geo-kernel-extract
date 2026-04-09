#!/usr/bin/env Rscript
# divergence.R — Classical vs MOPAC T2 divergence figures.
#
# Reads: output/secondary/divergence/*.csv
# Writes: output/secondary/figures/divergence_*.pdf
#
# Usage:
#   cd learn/R
#   Rscript divergence.R [data_dir]

args <- commandArgs(trailingOnly = TRUE)
source("common.R")
if (length(args) >= 1) set_data_dir(args[1])

div_dir <- file.path(DATA_DIR, "divergence")

# ── Load ─────────────────────────────────────────────────────────────
atom   <- read_csv(file.path(div_dir, "atom_divergence.csv"), show_col_types = FALSE)
by_elem <- read_csv(file.path(div_dir, "summary_by_element.csv"), show_col_types = FALSE)
by_rt   <- read_csv(file.path(div_dir, "summary_by_ring_type.csv"), show_col_types = FALSE)

# ── 1. Cosine similarity violin by element and pair ──────────────────
p1 <- ggplot(atom, aes(x = element, y = cosine_sim, fill = element)) +
  geom_violin(scale = "width", alpha = 0.7) +
  stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
  facet_wrap(~pair_name, nrow = 1) +
  scale_fill_manual(values = ELEMENT_COLOURS) +
  labs(title = "Classical vs MOPAC: T2 cosine similarity",
       subtitle = "1.0 = identical angular pattern, -1.0 = opposite",
       x = NULL, y = "Cosine similarity (5D)") +
  guides(fill = "none")

save_fig(p1, "divergence_cosine_by_element", width = 12, height = 5)

# ── 2. Divergence magnitude vs distance to ring ─────────────────────
p2 <- ggplot(atom[atom$pair_name == "coulomb_shielding", ],
             aes(x = ring_dist, y = diff_mag)) +
  geom_point(alpha = 0.05, size = 0.3, colour = "grey30") +
  geom_smooth(method = "loess", se = TRUE, colour = "#D55E00", linewidth = 0.8) +
  labs(title = "Coulomb shielding: classical-MOPAC T2 difference magnitude",
       x = "Distance to nearest removed ring (\u00c5)",
       y = "|T2(classical) - T2(MOPAC)| (ppm)") +
  coord_cartesian(xlim = c(0, 20))

save_fig(p2, "divergence_coulomb_vs_distance", width = 8, height = 5)

# ── 3. Mean divergence by ring type (bar chart) ─────────────────────
p3 <- ggplot(by_rt[by_rt$nearest_ring_type != "none", ],
             aes(x = reorder(nearest_ring_type, -diff_mag_mean),
                 y = diff_mag_mean, fill = nearest_ring_type)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = diff_mag_mean - diff_mag_std,
                     ymax = diff_mag_mean + diff_mag_std),
                width = 0.2) +
  facet_wrap(~pair_name, scales = "free_y") +
  ring_fill_scale() +
  labs(title = "Mean T2 divergence by nearest ring type",
       x = NULL, y = "Mean |T2 difference|") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p3, "divergence_by_ring_type", width = 10, height = 6)

# ── 4. Cosine similarity by ring type ───────────────────────────────
p4 <- ggplot(by_rt[by_rt$nearest_ring_type != "none", ],
             aes(x = nearest_ring_type, y = cosine_mean,
                 fill = nearest_ring_type)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = cosine_mean - cosine_std,
                     ymax = cosine_mean + cosine_std),
                width = 0.2) +
  facet_wrap(~pair_name) +
  ring_fill_scale() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "Classical-MOPAC T2 angular agreement by ring type",
       subtitle = "Higher = more similar angular pattern",
       x = NULL, y = "Mean cosine similarity") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p4, "divergence_cosine_by_ring_type", width = 10, height = 6)

message("Done: 4 divergence figures.")
