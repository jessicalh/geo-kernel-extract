#!/usr/bin/env Rscript
# twenty_eight_realities.R — Publication-quality thesis figures for NMR
# shielding tensor physics: distance falloff, element-wise R-squared,
# forward selection, angular peaking, BS/HM convergence.
#
# Reads:
#   output/actual_physics/atom_tensor_data.csv
#   output/secondary/element_physics/element_group_r2.csv
#   output/secondary/element_physics/forward_{H,C,N,O}.csv
#
# Writes: PDFs to output directory (argv[1], default FIG_DIR from common.R)
#
# Usage:
#   cd learn/R
#   Rscript twenty_eight_realities.R [output_dir]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
source("common.R")

# ── Output directory ─────────────────────────────────────────────────
OUT_DIR <- if (length(args) >= 1) args[1] else FIG_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Data paths ───────────────────────────────────────────────────────
# Navigate relative to the R/ directory up to learn/src/output
LEARN_DIR   <- dirname(getwd())  # learn/
ATOM_CSV    <- file.path(LEARN_DIR, "src", "output", "actual_physics",
                         "atom_tensor_data.csv")
ELEM_R2_CSV <- file.path(LEARN_DIR, "src", "output", "secondary",
                         "element_physics", "element_group_r2.csv")
FORWARD_DIR <- file.path(LEARN_DIR, "src", "output", "secondary",
                         "element_physics")

# ── Physics group palette ────────────────────────────────────────────
GROUP_COLOURS <- c(
  ring_current = "#0072B2",   # blue
  efg          = "#D55E00",   # red/vermillion
  bond_aniso   = "#009E73",   # green
  quadrupole   = "#CC79A7",   # purple
  dispersion   = "#E69F00",   # orange
  other        = "#999999"    # grey
)

group_fill_scale <- function(...) {
  scale_fill_manual(values = GROUP_COLOURS, ...)
}

group_colour_scale <- function(...) {
  scale_colour_manual(values = GROUP_COLOURS, ...)
}

# ── Helper: classify kernel name into physics group ──────────────────
classify_kernel <- function(name) {
  name <- tolower(name)
  dplyr::case_when(
    grepl("^(rbfbs|bs_|hm_|ringsusc)", name)   ~ "ring_current",
    grepl("^(efg|mopacefg|deltaapbs)", name)    ~ "efg",
    grepl("^(bond|aniso)", name)                 ~ "bond_aniso",
    grepl("^pq", name)                           ~ "quadrupole",
    grepl("^(disp|dispchi)", name)               ~ "dispersion",
    TRUE                                         ~ "other"
  )
}

# ── Helper: save a figure to OUT_DIR ─────────────────────────────────
save_pdf <- function(p, name, width = 8, height = 6) {
  path <- file.path(OUT_DIR, name)
  ggsave(path, p, width = width, height = height, device = cairo_pdf)
  message("Saved: ", path)
  invisible(path)
}

# ── Helper: safe file read ───────────────────────────────────────────
safe_read <- function(path) {
  if (!file.exists(path)) {
    message("  Skipping: file not found: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}


# ====================================================================
# 1. distance_falloff.pdf — Three-panel log-log: BS, HM, PQ vs dist
# ====================================================================
make_distance_falloff <- function() {
  message("Figure 1: distance_falloff")
  atom <- safe_read(ATOM_CSV)
  if (is.null(atom)) return(invisible())

  panels <- list(
    list(col = "bs_mag",  label = "Biot-Savart |T2|",   slope = -3),
    list(col = "hm_mag",  label = "Haigh-Mallion |T2|",  slope = -3),
    list(col = "pq_mag",  label = "Point-Quadrupole |T2|", slope = -5)
  )

  plot_list <- list()

  for (pp in panels) {
    d <- atom %>%
      select(dist, mag = all_of(pp$col)) %>%
      filter(dist > 2, mag > 1e-15) %>%
      mutate(log_dist = log10(dist), log_mag = log10(mag))

    # Fit slope in log-log space
    fit <- lm(log_mag ~ log_dist, data = d)
    measured_slope <- round(coef(fit)[2], 2)

    # Binned medians (20 log-spaced bins)
    d <- d %>%
      mutate(dist_bin = cut(log_dist, breaks = 20)) %>%
      group_by(dist_bin) %>%
      mutate(bin_mid_x = median(log_dist),
             bin_mid_y = median(log_mag)) %>%
      ungroup()

    bin_medians <- d %>%
      distinct(dist_bin, .keep_all = TRUE) %>%
      select(bin_mid_x, bin_mid_y)

    # Reference line: passes through data centroid with theoretical slope
    cx <- mean(d$log_dist)
    cy <- mean(d$log_mag)
    ref_intercept <- cy - pp$slope * cx

    p <- ggplot(d, aes(x = log_dist, y = log_mag)) +
      geom_point(alpha = 0.05, size = 0.3, colour = "grey50") +
      geom_point(data = bin_medians, aes(x = bin_mid_x, y = bin_mid_y),
                 colour = ELEMENT_COLOURS["N"], size = 2.5) +
      geom_abline(slope = pp$slope, intercept = ref_intercept,
                  linetype = "dashed", colour = "grey30", linewidth = 0.6) +
      labs(title = paste0(pp$label, "  (slope = ", measured_slope,
                          ", theory = ", pp$slope, ")"),
           x = bquote(log[10] * "(distance / \u00c5)"),
           y = bquote(log[10] * "(|T2|)"))

    plot_list <- c(plot_list, list(p))
  }

  combined <- plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
    plot_annotation(title = "Distance dependence of ring-current tensor magnitudes",
                    subtitle = "Grey scatter, binned medians (blue), theoretical slope (dashed)")

  save_pdf(combined, "distance_falloff.pdf", width = 7, height = 10)
}


# ====================================================================
# 2. element_r2.pdf — Grouped bar: R-squared by element x physics group
# ====================================================================
make_element_r2 <- function() {
  message("Figure 2: element_r2")
  er2 <- safe_read(ELEM_R2_CSV)
  if (is.null(er2)) return(invisible())

  # Keep H, C, N, O only for the bars; get pooled "all"
  pooled_all <- er2 %>% filter(element == "all") %>% pull(all)

  er2_main <- er2 %>%
    filter(element %in% c("H", "C", "N", "O")) %>%
    select(element, ring_current, efg, bond_aniso, all) %>%
    pivot_longer(cols = c(ring_current, efg, bond_aniso, all),
                 names_to = "group", values_to = "r2") %>%
    mutate(
      group = factor(group,
                     levels = c("ring_current", "efg", "bond_aniso", "all"),
                     labels = c("Ring current", "EFG", "Bond aniso", "All")),
      element = factor(element, levels = c("H", "C", "N", "O"))
    )

  # Assign fill colours: use physics palette + black for "all"
  bar_colours <- c(
    "Ring current" = unname(GROUP_COLOURS["ring_current"]),
    "EFG"          = unname(GROUP_COLOURS["efg"]),
    "Bond aniso"   = unname(GROUP_COLOURS["bond_aniso"]),
    "All"          = "#333333"
  )

  p <- ggplot(er2_main, aes(x = element, y = r2, fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_hline(yintercept = pooled_all, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    annotate("text", x = 4.3, y = pooled_all + 0.02,
             label = paste0("pooled = ", round(pooled_all, 3)),
             colour = "grey40", size = 3, hjust = 1) +
    scale_fill_manual(values = bar_colours, name = "Physics group") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(title = expression(R^2 ~ "by element and kernel physics group"),
         subtitle = "Mutation-set calibration, per-element forward model",
         x = "Element", y = expression(R^2))

  save_pdf(p, "element_r2.pdf", width = 7, height = 5)
}


# ====================================================================
# 3. forward_selection.pdf — 2x2 panel: cumulative R-squared staircase
# ====================================================================
make_forward_selection <- function() {
  message("Figure 3: forward_selection")

  elements <- c("H", "C", "N", "O")
  all_fwd <- list()

  for (el in elements) {
    fname <- file.path(FORWARD_DIR, paste0("forward_", el, ".csv"))
    d <- safe_read(fname)
    if (is.null(d)) next
    d$element <- el
    d$group <- classify_kernel(d$kernel)
    all_fwd[[el]] <- d
  }

  if (length(all_fwd) == 0) {
    message("  No forward selection data found, skipping.")
    return(invisible())
  }

  fwd <- bind_rows(all_fwd) %>%
    mutate(element = factor(element, levels = elements),
           group = factor(group, levels = names(GROUP_COLOURS)))

  # Short kernel labels for first 5
  fwd <- fwd %>%
    mutate(short_label = ifelse(rank <= 5,
                                sub("^(RbfBs|Rbf|Mopac|DeltaAPBS_)", "", kernel),
                                NA_character_))

  p <- ggplot(fwd %>% filter(rank <= 15),
              aes(x = rank, y = r2, colour = group)) +
    geom_step(linewidth = 0.8, direction = "hv") +
    geom_point(size = 2) +
    geom_text(aes(label = short_label), size = 2.5, vjust = -1.0, hjust = 0.5,
              na.rm = TRUE, show.legend = FALSE, colour = "grey20") +
    facet_wrap(~element, scales = "free_y", ncol = 2) +
    group_colour_scale(name = "Physics group",
                       labels = c("Ring current", "EFG", "Bond aniso",
                                  "Quadrupole", "Dispersion", "Other")) +
    scale_x_continuous(breaks = seq(1, 15, by = 2)) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    labs(title = expression("Forward kernel selection: cumulative" ~ R^2),
         subtitle = "First 15 kernels, coloured by physics group",
         x = "Selection rank", y = expression("Cumulative" ~ R^2))

  save_pdf(p, "forward_selection.pdf", width = 8, height = 7)
}


# ====================================================================
# 4. angular_peaking.pdf — BS magnitude vs theta, with theoretical curve
# ====================================================================
make_angular_peaking <- function() {
  message("Figure 4: angular_peaking")
  atom <- safe_read(ATOM_CSV)
  if (is.null(atom)) return(invisible())

  d <- atom %>%
    filter(bs_mag > 1e-15, dist > 2) %>%
    select(theta, bs_mag, dist)

  # Bin theta into 20 equal bins from 0 to pi/2
  breaks <- seq(0, pi / 2, length.out = 21)
  d <- d %>%
    mutate(theta_bin = cut(theta, breaks = breaks, include.lowest = TRUE)) %>%
    filter(!is.na(theta_bin))

  # Compute bin midpoints and median BS magnitude (distance-normalized)
  # Normalize by dist^3 to remove distance effect and isolate angular part
  d <- d %>%
    mutate(bs_norm = bs_mag * dist^3)

  bin_stats <- d %>%
    group_by(theta_bin) %>%
    summarise(
      mid_theta = mean(theta),
      median_bs = median(bs_norm),
      q25 = quantile(bs_norm, 0.25),
      q75 = quantile(bs_norm, 0.75),
      n = n(),
      .groups = "drop"
    )

  # Theoretical curve: proportional to (1 - 3 cos^2 theta)^2
  theta_seq <- seq(0, pi / 2, length.out = 200)
  theory <- data.frame(
    theta = theta_seq,
    envelope = (1 - 3 * cos(theta_seq)^2)^2
  )
  # Scale theory to match the median magnitude range
  theory$envelope <- theory$envelope *
    max(bin_stats$median_bs) / max(theory$envelope)

  p <- ggplot() +
    geom_col(data = bin_stats,
             aes(x = mid_theta, y = median_bs),
             fill = GROUP_COLOURS["ring_current"], alpha = 0.7,
             width = diff(breaks)[1] * 0.85) +
    geom_errorbar(data = bin_stats,
                  aes(x = mid_theta, ymin = q25, ymax = q75),
                  width = 0.02, colour = "grey40", linewidth = 0.4) +
    geom_line(data = theory, aes(x = theta, y = envelope),
              linetype = "dashed", colour = "#D55E00", linewidth = 0.8) +
    scale_x_continuous(
      breaks = c(0, pi / 6, pi / 4, pi / 3, pi / 2),
      labels = c("0",
                 expression(pi / 6),
                 expression(pi / 4),
                 expression(pi / 3),
                 expression(pi / 2))) +
    labs(title = bquote("BS magnitude angular dependence (" *
                          r^3 * "-normalised)"),
         subtitle = bquote("Binned medians (IQR whiskers). Dashed: " *
                             (1 - 3 * cos^2 * theta)^2 ~ "envelope"),
         x = bquote(theta ~ "(polar angle from ring normal)"),
         y = bquote("|BS|" %*% r^3 ~ "(arb. units)"))

  save_pdf(p, "angular_peaking.pdf", width = 7, height = 5)
}


# ====================================================================
# 5. bs_hm_ratio.pdf — Histogram of HM/BS magnitude ratio
# ====================================================================
make_bs_hm_ratio <- function() {
  message("Figure 5: bs_hm_ratio")
  atom <- safe_read(ATOM_CSV)
  if (is.null(atom)) return(invisible())

  d <- atom %>%
    filter(bs_mag > 1e-15, hm_mag > 1e-15) %>%
    mutate(ratio = hm_mag / bs_mag)

  med  <- median(d$ratio)
  sdev <- sd(d$ratio)

  # Restrict to a sensible plotting range
  d_plot <- d %>% filter(ratio > 0.95, ratio < 1.10)

  p <- ggplot(d_plot, aes(x = ratio)) +
    geom_histogram(bins = 80, fill = GROUP_COLOURS["ring_current"],
                   alpha = 0.7, colour = NA) +
    geom_vline(xintercept = 1.0, linetype = "dashed", colour = "grey30") +
    geom_vline(xintercept = med, linetype = "solid",
               colour = "#D55E00", linewidth = 0.7) +
    annotate("text", x = med + 0.005, y = Inf, vjust = 2,
             label = sprintf("median = %.4f\nsd = %.4f", med, sdev),
             colour = "#D55E00", size = 3.5, hjust = 0) +
    labs(title = "Haigh-Mallion / Biot-Savart magnitude ratio",
         subtitle = "Expected near 1.0; tight distribution confirms ring-model consistency",
         x = "|HM| / |BS|",
         y = "Count")

  save_pdf(p, "bs_hm_ratio.pdf", width = 7, height = 4.5)
}


# ====================================================================
# 6. efg_cos_by_element.pdf — Violin: efg_target_cos by element
# ====================================================================
make_efg_cos_by_element <- function() {
  message("Figure 6: efg_cos_by_element")
  atom <- safe_read(ATOM_CSV)
  if (is.null(atom)) return(invisible())

  d <- atom %>%
    filter(element %in% c("H", "C", "N", "O")) %>%
    filter(!is.na(efg_target_cos)) %>%
    mutate(element = factor(element, levels = c("H", "C", "N", "O")))

  # Median per element for annotation
  el_meds <- d %>%
    group_by(element) %>%
    summarise(med = median(efg_target_cos), .groups = "drop")

  p <- ggplot(d, aes(x = element, y = efg_target_cos, fill = element)) +
    geom_violin(alpha = 0.6, colour = "grey50", linewidth = 0.3,
                scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 alpha = 0.8, linewidth = 0.4) +
    geom_text(data = el_meds,
              aes(x = element, y = med, label = sprintf("%.3f", med)),
              vjust = -0.8, size = 3, colour = "grey20") +
    scale_fill_manual(values = ELEMENT_COLOURS, guide = "none") +
    labs(title = bquote("EFG-target cosine similarity" ~
                          (Delta * R^2 ~ "on deletion") ~ "by element"),
         subtitle = "Negative values = kernel contributes. Tighter = more consistent.",
         x = "Element",
         y = "cos(EFG, target)")

  save_pdf(p, "efg_cos_by_element.pdf", width = 6, height = 5)
}


# ====================================================================
# 7. bs_hm_convergence.pdf — |bs_hm_cos| vs distance
# ====================================================================
make_bs_hm_convergence <- function() {
  message("Figure 7: bs_hm_convergence")
  atom <- safe_read(ATOM_CSV)
  if (is.null(atom)) return(invisible())

  d <- atom %>%
    filter(!is.na(bs_hm_cos), bs_mag > 1e-15, hm_mag > 1e-15) %>%
    mutate(abs_cos = abs(bs_hm_cos))

  # Bin by distance
  breaks <- seq(2, max(d$dist, na.rm = TRUE), length.out = 40)
  d <- d %>%
    filter(dist >= 2) %>%
    mutate(dist_bin = cut(dist, breaks = breaks, include.lowest = TRUE))

  bin_stats <- d %>%
    filter(!is.na(dist_bin)) %>%
    group_by(dist_bin) %>%
    summarise(
      mid_dist   = median(dist),
      median_cos = median(abs_cos),
      q10        = quantile(abs_cos, 0.10),
      q90        = quantile(abs_cos, 0.90),
      n          = n(),
      .groups = "drop"
    )

  p <- ggplot() +
    geom_point(data = d %>% sample_frac(0.05),
               aes(x = dist, y = abs_cos),
               alpha = 0.03, size = 0.3, colour = "grey60") +
    geom_ribbon(data = bin_stats,
                aes(x = mid_dist, ymin = q10, ymax = q90),
                fill = GROUP_COLOURS["ring_current"], alpha = 0.25) +
    geom_line(data = bin_stats,
              aes(x = mid_dist, y = median_cos),
              colour = GROUP_COLOURS["ring_current"], linewidth = 0.8) +
    geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey40") +
    coord_cartesian(ylim = c(0.990, 1.001)) +
    labs(title = "|cos(BS, HM)| vs distance: convergence to parallel",
         subtitle = "Band: 10th-90th percentile. Median converges to 1.0 at large r.",
         x = "Distance from ring centre (\u00c5)",
         y = "|cos(BS, HM)|")

  save_pdf(p, "bs_hm_convergence.pdf", width = 7, height = 5)
}


# ====================================================================
# Run all figures
# ====================================================================
message("Output directory: ", OUT_DIR)
message("---")

make_distance_falloff()
make_element_r2()
make_forward_selection()
make_angular_peaking()
make_bs_hm_ratio()
make_efg_cos_by_element()
make_bs_hm_convergence()

message("---")
message("Done: twenty_eight_realities.R")
