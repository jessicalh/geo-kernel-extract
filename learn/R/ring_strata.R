#!/usr/bin/env Rscript
# ring_strata.R — Ring-type stratified analysis figures.
#
# Reads: output/secondary/ring_strata/*.csv
# Writes: output/secondary/figures/strata_*.pdf
#
# Usage:
#   cd learn/R
#   Rscript ring_strata.R [data_dir]

args <- commandArgs(trailingOnly = TRUE)
source("common.R")
if (length(args) >= 1) set_data_dir(args[1])

strata_dir <- file.path(DATA_DIR, "ring_strata")

# ── Load ─────────────────────────────────────────────────────────────
ptags   <- read_csv(file.path(strata_dir, "protein_tags.csv"), show_col_types = FALSE)
atags   <- read_csv(file.path(strata_dir, "atom_tags.csv"), show_col_types = FALSE)
kstrat  <- read_csv(file.path(strata_dir, "kernel_by_stratum.csv"), show_col_types = FALSE)
ridge   <- read_csv(file.path(strata_dir, "ridge_by_stratum.csv"), show_col_types = FALSE)

# ── 1. Ridge R² by stratum (bar chart) ──────────────────────────────
ridge_long <- tidyr::pivot_longer(
  ridge, cols = starts_with("ridge_r2_"),
  names_to = "kernel_family", values_to = "r2",
  names_prefix = "ridge_r2_"
)

p1 <- ggplot(ridge_long[!is.na(ridge_long$r2), ],
             aes(x = stratum, y = r2, fill = kernel_family)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", r2)),
            position = position_dodge(0.9), vjust = -0.3, size = 2.5) +
  labs(title = "Ridge R\u00b2 by atom stratum and kernel family",
       subtitle = paste0("n_atoms: ",
                         paste(ridge$stratum, ridge$n_atoms, sep = "=",
                               collapse = ", ")),
       x = NULL, y = "Ridge R\u00b2", fill = "Kernel family") +
  coord_cartesian(ylim = c(0, 1))

save_fig(p1, "strata_ridge_r2", width = 10, height = 6)

# ── 2. Top 15 kernels by T2 magnitude, faceted by stratum ───────────
# Pick top 15 by mean magnitude across all strata
top_kernels <- kstrat |>
  dplyr::group_by(kernel_name) |>
  dplyr::summarise(overall = mean(t2_mag_mean), .groups = "drop") |>
  dplyr::slice_max(overall, n = 15) |>
  dplyr::pull(kernel_name)

kstrat_top <- kstrat[kstrat$kernel_name %in% top_kernels, ]

p2 <- ggplot(kstrat_top,
             aes(x = reorder(kernel_name, -t2_mag_mean),
                 y = t2_mag_mean, fill = stratum)) +
  geom_col(position = "dodge", alpha = 0.8) +
  stratum_fill_scale() +
  labs(title = "Top 15 kernels: T2 magnitude by stratum",
       x = NULL, y = "Mean T2 magnitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p2, "strata_kernel_magnitude", width = 12, height = 6)

# ── 3. Atom counts by stratum (with ring type labels) ───────────────
stratum_counts <- as.data.frame(table(atags$stratum))
names(stratum_counts) <- c("stratum", "n_atoms")

p3 <- ggplot(stratum_counts,
             aes(x = reorder(stratum, -n_atoms), y = n_atoms,
                 fill = stratum)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = n_atoms), vjust = -0.3, size = 3) +
  stratum_fill_scale() +
  labs(title = "Atom counts by ring-type stratum",
       x = NULL, y = "Number of matched atoms") +
  guides(fill = "none")

save_fig(p3, "strata_atom_counts", width = 8, height = 5)

# ── 4. Protein ring type composition ────────────────────────────────
has_cols <- grep("^has_", names(ptags), value = TRUE)
prot_long <- tidyr::pivot_longer(
  ptags, cols = all_of(has_cols),
  names_to = "ring_type", values_to = "present",
  names_prefix = "has_"
)
prot_long <- prot_long[prot_long$present == 1, ]

p4 <- ggplot(prot_long, aes(x = ring_type, fill = ring_type)) +
  geom_bar(alpha = 0.8) +
  ring_fill_scale() +
  labs(title = "Number of proteins containing each ring type",
       x = NULL, y = "Protein count") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p4, "strata_protein_composition", width = 8, height = 5)

message("Done: 4 ring strata figures.")
