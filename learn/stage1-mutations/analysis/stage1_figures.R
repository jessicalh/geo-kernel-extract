#!/usr/bin/env Rscript
# stage1_figures.R — Publication figures for Stage 1 thesis chapter.
#
# Reads JSON and CSV output from src/actual_physics/ scripts.
# Produces PDFs in stage1-mutations/figures/.
#
# Usage:
#   cd learn
#   Rscript stage1-mutations/analysis/stage1_figures.R
#
# Prerequisite: run_all.sh must have completed (all JSONs current).

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

# ── Paths ────────────────────────────────────────────────────────
LEARN    <- getwd()
OUT_DATA <- file.path(LEARN, "src", "output", "actual_physics")
FIG_DIR  <- file.path(LEARN, "stage1-mutations", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ELEMENT_COLOURS <- c(H = "#0072B2", C = "#D55E00", N = "#009E73", O = "#CC79A7")

save_pdf <- function(p, name, width = 8, height = 6) {
  path <- file.path(FIG_DIR, name)
  ggsave(path, p, width = width, height = height, device = cairo_pdf)
  message("Saved: ", path)
}

# ── 1. Cosine heatmap per element ────────────────────────────────
make_cosine_heatmaps <- function() {
  message("Figure: cosine heatmaps")

  plots <- list()
  for (el in c("H", "C", "N", "O")) {
    csv <- file.path(OUT_DATA, "full_space", paste0("cosine_matrix_", el, ".csv"))
    if (!file.exists(csv)) { message("  Skipping ", el); next }

    mat <- read_csv(csv, show_col_types = FALSE)
    kernel_names <- mat[[1]]
    mat <- mat[, -1]
    colnames(mat) <- kernel_names

    # Physics group classification
    classify <- function(n) {
      n <- as.character(n)
      if (grepl("^(BS_|HM_|RingSusc)", n)) return("Ring current")
      if (grepl("^Disp_", n)) return("Dispersion")
      if (grepl("^(PQ_|PQ$)", n)) return("Quadrupole")
      if (grepl("^MC_|^HBond", n)) return("Bond aniso")
      if (grepl("^MopacMC", n)) return("MOPAC bond")
      if (n %in% c("Coulomb_total", "EFG_bb", "EFG_aro")) return("ff14SB EFG")
      if (grepl("^Mopac(EFG|Coulomb)", n)) return("MOPAC EFG")
      if (grepl("^(APBS|Delta)", n)) return("Solvation")
      return("Other")
    }

    groups <- sapply(kernel_names, classify)

    # Long format for ggplot
    long <- mat %>%
      mutate(k1 = kernel_names, g1 = groups) %>%
      pivot_longer(cols = all_of(kernel_names), names_to = "k2", values_to = "cos") %>%
      mutate(g2 = groups[match(k2, kernel_names)],
             cos = as.numeric(cos))

    # Order by group
    group_order <- c("Ring current", "ff14SB EFG", "MOPAC EFG",
                     "Bond aniso", "MOPAC bond", "Quadrupole",
                     "Dispersion", "Solvation")
    long <- long %>%
      mutate(k1 = factor(k1, levels = kernel_names[order(match(groups, group_order))]),
             k2 = factor(k2, levels = kernel_names[order(match(groups, group_order))]))

    p <- ggplot(long, aes(x = k2, y = k1, fill = cos)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = 0.45, limits = c(0.25, 1.0),
                           name = "|cos|") +
      labs(title = paste0(el, ": inter-kernel T2 cosine similarity"),
           x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
            axis.text.y = element_text(size = 4))

    plots[[el]] <- p
    save_pdf(p, paste0("cosine_heatmap_", el, ".pdf"), width = 10, height = 9)
  }
}


# ── 2. Per-group R² bar chart: raw vs normalised ────────────────
make_group_r2_comparison <- function() {
  message("Figure: group R² raw vs normalised")

  fsa <- fromJSON(file.path(OUT_DATA, "full_space", "full_space_analysis.json"))

  groups <- c("ring_current", "ff14sb_efg", "mopac_efg", "bond_aniso",
              "mopac_bond", "quadrupole", "dispersion", "solvation")
  group_labels <- c("Ring current", "ff14SB EFG", "MOPAC EFG", "Bond aniso",
                    "MOPAC bond", "Quadrupole", "Dispersion", "Solvation")

  rows <- list()
  for (el in c("H", "C", "N", "O")) {
    for (i in seq_along(groups)) {
      g <- groups[i]
      r2_raw <- fsa[[el]]$per_group_r2$raw[[g]]
      r2_norm <- fsa[[el]]$per_group_r2$norm[[g]]
      if (is.null(r2_raw)) next
      rows <- c(rows, list(data.frame(
        element = el, group = group_labels[i],
        raw = r2_raw, norm = r2_norm,
        stringsAsFactors = FALSE
      )))
    }
  }
  df <- bind_rows(rows)

  df_long <- df %>%
    pivot_longer(cols = c(raw, norm), names_to = "space", values_to = "r2") %>%
    mutate(
      element = factor(element, levels = c("H", "C", "N", "O")),
      group = factor(group, levels = group_labels),
      space = factor(space, levels = c("raw", "norm"),
                     labels = c("Raw", "Normalised"))
    )

  p <- ggplot(df_long, aes(x = group, y = r2, fill = space)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~element, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c(Raw = "#999999", Normalised = "#0072B2"),
                      name = "Space") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = expression("Per-group ridge" ~ R^2 ~ ": raw vs normalised"),
         subtitle = "Per-protein normalisation reveals dispersion for O, helps N overall",
         x = NULL, y = expression(R^2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  save_pdf(p, "group_r2_raw_vs_norm.pdf", width = 10, height = 8)
}


# ── 3. Dimensionality PCA-ridge curves ──────────────────────────
make_dimensionality_curves <- function() {
  message("Figure: PCA-ridge curves")

  dim <- fromJSON(file.path(OUT_DATA, "dimensionality", "dimensionality.json"))

  rows <- list()
  for (el in c("H", "C", "N", "O")) {
    curve <- dim[[el]]$pca_ridge$curve
    for (i in seq_len(nrow(curve))) {
      rows <- c(rows, list(data.frame(
        element = el, k = curve[i, 1],
        train = curve[i, 2], val = curve[i, 3],
        stringsAsFactors = FALSE
      )))
    }
  }
  df <- bind_rows(rows) %>%
    mutate(element = factor(element, levels = c("H", "C", "N", "O")))

  p <- ggplot(df, aes(x = k)) +
    geom_line(aes(y = train, colour = "Train"), linewidth = 0.6) +
    geom_line(aes(y = val, colour = "Validation"), linewidth = 0.8) +
    facet_wrap(~element, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = c(Train = "#999999", Validation = "#0072B2"),
                        name = NULL) +
    labs(title = "PCA-then-ridge: predictive dimensionality after normalisation",
         subtitle = "Validation R² plateaus at H=20, C=6, N=3, O=12",
         x = "Number of PCA components (k)",
         y = expression("Validation" ~ R^2)) +
    theme_minimal()

  save_pdf(p, "pca_ridge_curves.pdf", width = 9, height = 7)
}


# ── 4. Forward selection: group contribution staircase ───────────
make_forward_selection <- function() {
  message("Figure: forward selection by group")

  fsa <- fromJSON(file.path(OUT_DATA, "full_space", "full_space_analysis.json"))

  group_colours <- c(
    ring_current = "#0072B2", ff14sb_efg = "#D55E00", mopac_efg = "#E69F00",
    bond_aniso = "#009E73", mopac_bond = "#56B4E9", quadrupole = "#CC79A7",
    dispersion = "#F0E442", solvation = "#999999"
  )

  rows <- list()
  for (el in c("H", "C", "N", "O")) {
    fwd <- fsa[[el]]$forward_selection$norm
    for (i in seq_along(fwd)) {
      step <- fwd[[i]]
      rows <- c(rows, list(data.frame(
        element = el, rank = step$rank,
        kernel = step$kernel, group = step$group,
        delta = step$delta, cumulative = step$cumulative,
        stringsAsFactors = FALSE
      )))
    }
  }
  df <- bind_rows(rows) %>%
    mutate(element = factor(element, levels = c("H", "C", "N", "O")),
           group = factor(group, levels = names(group_colours)))

  p <- ggplot(df %>% filter(rank <= 15),
              aes(x = rank, y = cumulative, colour = group)) +
    geom_step(linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap(~element, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = group_colours, name = "Physics group",
                        labels = c("Ring current", "ff14SB EFG", "MOPAC EFG",
                                   "Bond aniso", "MOPAC bond", "Quadrupole",
                                   "Dispersion", "Solvation")) +
    labs(title = "Forward kernel selection (normalised): which groups contribute",
         subtitle = "H = ring current dominated; N = all groups mixed; O = dispersion dominated",
         x = "Selection rank", y = expression("Cumulative" ~ R^2)) +
    theme_minimal()

  save_pdf(p, "forward_selection_norm.pdf", width = 10, height = 8)
}


# ── Run all ──────────────────────────────────────────────────────
message("Output directory: ", FIG_DIR)
message("---")

make_cosine_heatmaps()
make_group_r2_comparison()
make_dimensionality_curves()
make_forward_selection()

message("---")
message("Done: stage1_figures.R")
