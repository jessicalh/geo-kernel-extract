#!/usr/bin/env Rscript

options(encoding = "UTF-8")
try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

capstone_dir <- "/tmp/rediscover-corrected-backbone-snapshot-1p9j/analysis"
mopac_dir <- "/tmp/rediscover-task1-mopac-field/analysis"

args <- commandArgs(trailingOnly = TRUE)
default_out <- file.path("/tmp", paste0("rediscover-rich-chart-suite-", format(Sys.time(), "%Y%m%d-%H%M%S")))
out_dir <- if (length(args) >= 1) args[[1]] else default_out
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  cap_var = file.path(capstone_dir, "variance_decomposition/variance_decomposition.csv"),
  cap_static = file.path(capstone_dir, "static_environment_calibration/static_environment_calibration.csv"),
  cap_ring = file.path(capstone_dir, "ring_literature_decirc/ring_literature_decirc.csv"),
  cap_mcconnell = file.path(capstone_dir, "mcconnell_literature_decirc/mcconnell_literature_decirc.csv"),
  cap_aimnet2 = file.path(capstone_dir, "aimnet2_ensemble/aimnet2_ceiling_ensemble.csv"),
  mopac_var = file.path(mopac_dir, "variance_decomposition/variance_decomposition.csv"),
  mopac_static = file.path(mopac_dir, "static_environment_calibration/static_environment_calibration.csv")
)

missing <- paths[!file.exists(paths)]
if (length(missing) > 0) {
  stop("Missing required emitted result CSV(s):\n", paste(missing, collapse = "\n"))
}

read_result <- function(path) {
  readr::read_csv(
    path,
    na = c("", "NA", "NaN", "nan"),
    locale = readr::locale(encoding = "UTF-8"),
    show_col_types = FALSE
  )
}

cap_var <- read_result(paths["cap_var"])
cap_static <- read_result(paths["cap_static"])
ring_decirc <- read_result(paths["cap_ring"])
mc_decirc <- read_result(paths["cap_mcconnell"])
aimnet2 <- read_result(paths["cap_aimnet2"])
mopac_var <- read_result(paths["mopac_var"])
mopac_static <- read_result(paths["mopac_static"])

gamma_sym <- "\u03B3"
sigma_sym <- "\u03C3"
sup2 <- "\u00B2"
glyph_decirc <- "\u25CF"
glyph_form <- "\u25D0"
glyph_cant <- "\u25CB"

strata_order <- c("N", "CA", "C", "O", "HN", "HA", "Aro-H")
backbone_order <- c("N", "CA", "C", "O", "HN", "HA")
mechanism_order <- c(
  "ring", "charge", "McConnell", "MOPAC-Coulomb-EFG",
  "MOPAC-Mc", "Buckingham", "APBS-EFG", "AIMNet2 ceiling"
)

theme_report <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "#d7d7d7"),
      axis.title = element_text(colour = "#202020"),
      axis.text = element_text(colour = "#202020"),
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(colour = "#3a3a3a", size = base_size),
      plot.caption = element_text(colour = "#555555", size = base_size - 2, hjust = 0),
      strip.text = element_text(face = "bold", colour = "#202020"),
      legend.position = "bottom"
    )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_neff <- function(x) {
  case_when(
    is.na(x) ~ "n/a",
    abs(x) >= 1000 ~ paste0(formatC(x / 1000, format = "f", digits = 1), "k"),
    TRUE ~ formatC(x, format = "f", digits = 1)
  )
}

wrap_caption <- function(x, width = 145) {
  paste(strwrap(x, width = width), collapse = "\n")
}

grade_from_bucket <- function(bucket) {
  case_when(
    is.na(bucket) ~ "can't-work-for-now",
    grepl("zero-circularity", bucket) ~ "de-circularised",
    grepl("form-recovered", bucket) ~ "form-recovered-scale-fitted",
    TRUE ~ "can't-work-for-now"
  )
}

glyph_for_grade <- function(grade) {
  case_when(
    grade == "de-circularised" ~ glyph_decirc,
    grade == "form-recovered-scale-fitted" ~ glyph_form,
    TRUE ~ glyph_cant
  )
}

save_pdf <- function(plot, stem, width, height) {
  path <- file.path(out_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  path
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

mc_law <- mc_decirc %>%
  filter(target == "T2", axis == "within", stratum %in% backbone_order) %>%
  transmute(stratum, law_grade = grade_from_bucket(verdict_bucket))

cap_t2_cells <- cap_var %>%
  filter(
    target == "T2",
    stratum %in% backbone_order,
    mechanism %in% c("ring_current_T2", "ff14sb_charge_EFG_T2", "bond_anisotropy_T2")
  ) %>%
  mutate(
    mechanism_display = recode(
      mechanism,
      ring_current_T2 = "ring",
      ff14sb_charge_EFG_T2 = "charge",
      bond_anisotropy_T2 = "McConnell"
    ),
    law_grade = case_when(
      mechanism_display == "ring" ~ "form-recovered-scale-fitted",
      mechanism_display == "charge" ~ "form-recovered-scale-fitted",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(mc_law, by = "stratum", suffix = c("", "_mc")) %>%
  mutate(
    law_grade = if_else(mechanism_display == "McConnell" & !is.na(law_grade_mc), law_grade_mc, law_grade),
    metric_value = within_frame_absT2_r,
    metric_code = "C rand |T2| r",
    metric_family = "capstone random within-frame |T2| r",
    n_eff = within_N_eff_median,
    source = "corrected capstone variance_decomposition.csv"
  ) %>%
  transmute(stratum, mechanism = mechanism_display, metric_value, metric_code, metric_family,
            n_eff, law_grade, source)

mopac_field_cells <- mopac_var %>%
  filter(
    target == "T2",
    split_strategy == "blocked",
    stratum %in% backbone_order,
    mechanism %in% c("apbs_efg_T2", "mopac_coulomb_shielding_T2", "mopac_mc_shielding_T2")
  ) %>%
  mutate(
    mechanism_display = recode(
      mechanism,
      apbs_efg_T2 = "APBS-EFG",
      mopac_coulomb_shielding_T2 = "MOPAC-Coulomb-EFG",
      mopac_mc_shielding_T2 = "MOPAC-Mc"
    ),
    metric_value = within_frame_absT2_r,
    metric_code = "M blocked |T2| r",
    metric_family = "MOPAC field blocked within-frame |T2| r",
    n_eff = within_N_eff_median,
    law_grade = "form-recovered-scale-fitted",
    source = "MOPAC field variance_decomposition.csv"
  ) %>%
  transmute(stratum, mechanism = mechanism_display, metric_value, metric_code, metric_family,
            n_eff, law_grade, source)

buckingham_cells <- cap_var %>%
  filter(
    target == "sigma_iso",
    stratum %in% backbone_order,
    mechanism == "apbs_efield_buckingham"
  ) %>%
  mutate(
    mechanism = "Buckingham",
    metric_value = within_frame_R2,
    metric_code = paste0("C rand ", sigma_sym, " R", sup2),
    metric_family = paste0("capstone random within-frame ", sigma_sym, " R", sup2),
    n_eff = within_N_eff_median,
    law_grade = "form-recovered-scale-fitted",
    source = "corrected capstone variance_decomposition.csv"
  ) %>%
  transmute(stratum, mechanism, metric_value, metric_code, metric_family,
            n_eff, law_grade, source)

aimnet2_cells <- aimnet2 %>%
  filter(target == "T2", feature_set == "ensemble_plus_aimnet2", stratum %in% backbone_order) %>%
  mutate(
    mechanism = "AIMNet2 ceiling",
    metric_value = within_frame_R2,
    metric_code = paste0("AI T2 R", sup2),
    metric_family = paste0("AIMNet2 ceiling within-frame R", sup2),
    n_eff = NA_real_,
    law_grade = "can't-work-for-now",
    source = "aimnet2_ceiling_ensemble.csv"
  ) %>%
  transmute(stratum, mechanism, metric_value, metric_code, metric_family,
            n_eff, law_grade, source)

ring_aro_cell <- ring_decirc %>%
  filter(target == "T2", axis == "within", ring_type == "all_valid") %>%
  slice(1) %>%
  mutate(
    stratum = "Aro-H",
    mechanism = "ring",
    metric_value = absT2_r,
    metric_code = "C de-circ |T2| r",
    metric_family = "aromatic-H ring de-circularised within |T2| r",
    n_eff = atom_signal_neff,
    law_grade = grade_from_bucket(literature_verdict_bucket),
    source = "ring_literature_decirc.csv"
  ) %>%
  transmute(stratum, mechanism, metric_value, metric_code, metric_family,
            n_eff, law_grade, source)

heat_cells <- bind_rows(
  cap_t2_cells, mopac_field_cells, buckingham_cells, aimnet2_cells, ring_aro_cell
) %>%
  mutate(
    recovery = clip01(metric_value),
    glyph = glyph_for_grade(law_grade),
    label = paste0(glyph, " ", fmt_num(metric_value), "\n", metric_code, "\nN_eff ", fmt_neff(n_eff))
  )

heat_full <- expand_grid(stratum = strata_order, mechanism = mechanism_order) %>%
  left_join(heat_cells, by = c("stratum", "mechanism")) %>%
  mutate(
    stratum = factor(stratum, levels = rev(strata_order)),
    mechanism = factor(mechanism, levels = mechanism_order),
    label = replace_na(label, ""),
    recovery = if_else(is.na(recovery), NA_real_, recovery),
    ceiling_cell = mechanism == "AIMNet2 ceiling" & !is.na(metric_value)
  )

heat_plot <- ggplot(heat_full, aes(x = mechanism, y = stratum)) +
  geom_tile(aes(fill = recovery), colour = "#f7f7f7", linewidth = 0.55) +
  geom_tile(
    data = filter(heat_full, ceiling_cell),
    fill = NA,
    colour = "#111111",
    linewidth = 1.2
  ) +
  geom_text(aes(label = label), size = 2.35, lineheight = 0.9, colour = "#151515") +
  scale_fill_gradientn(
    colours = c("#f3f3f3", "#d8e7b0", "#88c4b5", "#3d83a5", "#17324d"),
    values = scales::rescale(c(0, 0.15, 0.35, 0.6, 1)),
    limits = c(0, 1),
    na.value = "#eeeeee",
    name = "recovery\n(clipped 0-1)"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_equal() +
  labs(
    title = "Combined rediscover recovery grid",
    subtitle = paste0(
      "Rows are backbone strata plus aromatic H. Cells are annotated by metric family; ",
      "the AIMNet2 column is a ceiling band, not a physical law."
    ),
    x = NULL,
    y = NULL,
    caption = wrap_caption(paste0(
      glyph_decirc, " de-circularised; ", glyph_form, " form recovered/scale fitted; ",
      glyph_cant, " ceiling or not a usable law. ",
      "C rand = corrected capstone random within-frame; M blocked = MOPAC field blocked within-frame. ",
      "Do not compare |T2| r, ", sigma_sym, " R", sup2, ", and AIMNet2 R", sup2,
      " as identical effect sizes."
    ))
  ) +
  theme_report(10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

field_map <- tibble(
  mechanism = c("apbs_efg_T2", "mopac_coulomb_shielding_T2", "mopac_mc_shielding_T2"),
  mechanism_display = c("APBS-EFG", "MOPAC-Coulomb-EFG", "MOPAC-Mc")
)

field_base <- mopac_var %>%
  filter(split_strategy == "blocked", mechanism %in% field_map$mechanism, stratum %in% backbone_order) %>%
  inner_join(field_map, by = "mechanism") %>%
  mutate(
    stratum = factor(stratum, levels = backbone_order),
    mechanism_display = factor(mechanism_display, levels = field_map$mechanism_display),
    n_eff = within_N_eff_median,
    r_value = within_frame_absT2_r,
    r_for_ci = pmin(pmax(r_value, -0.999999), 0.999999),
    z = atanh(r_for_ci),
    z_se = 1 / sqrt(pmax(n_eff - 3, 1)),
    r_lo = tanh(z - 1.96 * z_se),
    r_hi = tanh(z + 1.96 * z_se)
  )

field_long <- bind_rows(
  field_base %>%
    transmute(
      stratum, mechanism_display, n_eff,
      metric = "blocked within-frame |T2| r",
      value = r_value, lo = r_lo, hi = r_hi,
      label = paste0("N_eff ", fmt_neff(n_eff))
    ),
  field_base %>%
    transmute(
      stratum, mechanism_display, n_eff,
      metric = paste0("blocked within-frame R", sup2),
      value = within_frame_R2, lo = NA_real_, hi = NA_real_,
      label = paste0("N_eff ", fmt_neff(n_eff))
    )
) %>%
  mutate(metric = factor(metric, levels = c("blocked within-frame |T2| r", paste0("blocked within-frame R", sup2))))

field_plot <- ggplot(field_long, aes(x = stratum, y = value, colour = mechanism_display, shape = mechanism_display)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#6a6a6a") +
  geom_errorbar(
    data = filter(field_long, !is.na(lo)),
    aes(ymin = lo, ymax = hi),
    width = 0.12,
    position = position_dodge(width = 0.55),
    linewidth = 0.55
  ) +
  geom_point(
    aes(size = n_eff),
    position = position_dodge(width = 0.55),
    stroke = 0.6
  ) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("APBS-EFG" = "#7a7a7a", "MOPAC-Coulomb-EFG" = "#1967a3", "MOPAC-Mc" = "#c46b27")) +
  scale_size_continuous(range = c(2.2, 4.8), labels = fmt_neff, name = "N_eff") +
  labs(
    title = "Field-leg recovery: APBS versus MOPAC sidecars",
    subtitle = paste0(
      "MOPAC field run only. Top-panel whiskers are Fisher-z 95% CI from emitted r and N_eff; ",
      "R", sup2, " CI is not emitted."
    ),
    x = "stratum",
    y = "emitted blocked within-frame metric",
    colour = "field leg",
    shape = "field leg",
    caption = "Correlate-not-match: MOPAC-Coulomb and MOPAC-Mc are sidecar correlations, not raw EFG equality claims."
  ) +
  theme_report(10)

ring_gamma <- ring_decirc %>%
  filter(target == "T2", axis %in% c("within", "between"), ring_type == "all_valid") %>%
  transmute(
    mechanism = "ring",
    stratum = "Aro-H",
    fit_axis = axis,
    gamma = gamma_literature_scaled,
    se = gamma_literature_scaled_se,
    law_grade = grade_from_bucket(literature_verdict_bucket),
    source_class = "de-circularised"
  )

mc_gamma <- mc_decirc %>%
  filter(target == "T2", axis %in% c("within", "between"), stratum %in% backbone_order) %>%
  transmute(
    mechanism = "McConnell",
    stratum,
    fit_axis = axis,
    gamma = gamma_lit,
    se = gamma_lit_se,
    law_grade = grade_from_bucket(verdict_bucket),
    source_class = "de-circularised"
  )

scale_gamma <- bind_rows(
  cap_static %>%
    filter(primary_coef_name == "gamma", mechanism %in% c("ff14sb_charge_EFG_T2", "apbs_charge_EFG_T2")) %>%
    mutate(mechanism_display = recode(
      mechanism,
      ff14sb_charge_EFG_T2 = "charge",
      apbs_charge_EFG_T2 = "APBS-EFG"
    )),
  mopac_static %>%
    filter(primary_coef_name == "gamma", mechanism %in% c("mopac_coulomb_shielding_T2", "mopac_mc_shielding_T2")) %>%
    mutate(mechanism_display = recode(
      mechanism,
      mopac_coulomb_shielding_T2 = "MOPAC-Coulomb-EFG",
      mopac_mc_shielding_T2 = "MOPAC-Mc"
    ))
) %>%
  filter(stratum %in% backbone_order) %>%
  transmute(
    mechanism = mechanism_display,
    stratum,
    fit_axis = "scale fit",
    gamma = primary_coef,
    se = primary_se,
    law_grade = "form-recovered-scale-fitted",
    source_class = "scale-fitted"
  )

gamma_df <- bind_rows(ring_gamma, mc_gamma, scale_gamma) %>%
  mutate(
    stratum = factor(stratum, levels = c(backbone_order, "Aro-H")),
    mechanism = factor(mechanism, levels = mechanism_order),
    fit_axis = factor(fit_axis, levels = c("within", "between", "scale fit")),
    ci_lo = gamma - 1.96 * se,
    ci_hi = gamma + 1.96 * se,
    glyph = glyph_for_grade(law_grade)
  )

gamma_breaks <- c(-1000, -300, -100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100, 300, 1000)

gamma_plot <- ggplot(gamma_df, aes(x = stratum, y = gamma, colour = law_grade, shape = fit_axis, group = fit_axis)) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", colour = "#333333") +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.16,
    linewidth = 0.45,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(size = 2.2, stroke = 0.7, position = position_dodge(width = 0.45)) +
  facet_wrap(~mechanism, ncol = 2, scales = "free_x") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10), breaks = gamma_breaks) +
  scale_colour_manual(
    values = c(
      "de-circularised" = "#1f6f5b",
      "form-recovered-scale-fitted" = "#b46a1f",
      "can't-work-for-now" = "#777777"
    )
  ) +
  labs(
    title = paste0(gamma_sym, " recovery against the literature reference"),
    subtitle = paste0(
      "Dashed line is ", gamma_sym, " = 1. Axis uses a signed pseudo-log scale so scale-fitted drift remains visible. ",
      "Buckingham is absent because no dimensionless ", gamma_sym, " is emitted."
    ),
    x = "stratum",
    y = paste0(gamma_sym, " with 95% CI from emitted SE"),
    colour = "law grade",
    shape = "fit axis",
    caption = "De-circ rows use literature-de-circularising CSVs; scale-fit rows use static calibration CSVs. No physics or trajectory access."
  ) +
  theme_report(10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

axis_cells <- bind_rows(
  cap_var %>%
    filter(
      target == "T2",
      stratum %in% backbone_order,
      mechanism %in% c("ring_current_T2", "ff14sb_charge_EFG_T2", "bond_anisotropy_T2", "apbs_efg_T2")
    ) %>%
    mutate(
      source = "corrected capstone",
      mechanism_display = recode(
        mechanism,
        ring_current_T2 = "ring",
        ff14sb_charge_EFG_T2 = "charge",
        bond_anisotropy_T2 = "McConnell",
        apbs_efg_T2 = "APBS-EFG"
      )
    ),
  mopac_var %>%
    filter(
      target == "T2",
      split_strategy == "blocked",
      stratum %in% backbone_order,
      mechanism %in% c("apbs_efg_T2", "mopac_coulomb_shielding_T2", "mopac_mc_shielding_T2")
    ) %>%
    mutate(
      source = "MOPAC field blocked",
      mechanism_display = recode(
        mechanism,
        apbs_efg_T2 = "APBS-EFG",
        mopac_coulomb_shielding_T2 = "MOPAC-Coulomb-EFG",
        mopac_mc_shielding_T2 = "MOPAC-Mc"
      )
    )
) %>%
  transmute(
    source,
    mechanism = mechanism_display,
    stratum,
    between_r = between_LOAO_absT2_r,
    within_r = within_frame_absT2_r,
    n_eff = within_N_eff_median
  ) %>%
  mutate(
    source = factor(source, levels = c("corrected capstone", "MOPAC field blocked")),
    mechanism = factor(mechanism, levels = mechanism_order),
    stratum = factor(stratum, levels = backbone_order)
  )

axis_plot <- ggplot(axis_cells, aes(x = between_r, y = within_r, colour = mechanism)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#7a7a7a") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#7a7a7a") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = "dotted", colour = "#909090") +
  geom_point(aes(size = n_eff), alpha = 0.88) +
  geom_text(aes(label = stratum), size = 2.4, nudge_y = 0.025, show.legend = FALSE) +
  facet_wrap(~source, ncol = 2) +
  scale_size_continuous(range = c(2.2, 5), labels = fmt_neff, name = "within N_eff") +
  labs(
    title = "Static-axis versus dynamic-axis recovery",
    subtitle = "Standalone sidecar correlations: between atom means on x; within-frame motion on y.",
    x = "between LOAO |T2| r",
    y = "within-frame |T2| r",
    colour = "mechanism",
    caption = "This view separates axis placement; it is not a joint model score."
  ) +
  theme_report(10)

aim_labels <- c(
  ensemble_physics_only = "physics only",
  ensemble_plus_crg = "+ CRG",
  ensemble_plus_embedding = "+ embedding",
  ensemble_plus_aimnet2_charge = "+ AIMNet2 charge",
  ensemble_plus_crg_embedding = "+ CRG + embedding",
  ensemble_plus_aimnet2 = "full AIMNet2"
)

aim_plot_df <- aimnet2 %>%
  filter(target == "T2", stratum %in% backbone_order, feature_set %in% names(aim_labels)) %>%
  mutate(
    stratum = factor(stratum, levels = backbone_order),
    feature_label = factor(recode(feature_set, !!!aim_labels), levels = aim_labels),
    label = paste0(
      "R", sup2, " ", fmt_num(within_frame_R2), "\n",
      "lift ", fmt_num(100 * within_lift_fraction_of_aimnet2_full, 0), "%"
    )
  )

aim_plot <- ggplot(aim_plot_df, aes(x = stratum, y = feature_label, fill = within_frame_R2)) +
  geom_tile(colour = "#f7f7f7", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.55, lineheight = 0.9, colour = "#111111") +
  scale_fill_gradientn(
    colours = c("#f4f4f4", "#d3e8c5", "#78c3bb", "#2e7fa5", "#0f334d"),
    values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.75)),
    limits = c(0, 0.75),
    name = paste0("within R", sup2)
  ) +
  labs(
    title = "AIMNet2 ceiling bands and feature lifts",
    subtitle = "Standalone physics row versus joint feature additions; lift is fraction of full AIMNet2 within-frame gain.",
    x = "stratum",
    y = NULL,
    caption = "Ceiling model context only; these cells are not physical-law recoveries."
  ) +
  theme_report(10) +
  theme(panel.grid = element_blank(), legend.position = "right")

chart_paths <- c(
  combined_set_heatmap = save_pdf(heat_plot, "combined_set_heatmap", 14.2, 8.4),
  field_leg_blocked_panel = save_pdf(field_plot, "field_leg_blocked_panel", 10.8, 7.6),
  decirc_gamma_panels = save_pdf(gamma_plot, "decirc_gamma_panels", 12.4, 9.8),
  axis_recovery_quadrants = save_pdf(axis_plot, "axis_recovery_quadrants", 11.2, 5.8),
  aimnet2_ceiling_lift = save_pdf(aim_plot, "aimnet2_ceiling_lift", 10.4, 6.2)
)

readr::write_csv(heat_cells, file.path(out_dir, "combined_set_heatmap_data.csv"), na = "")
readr::write_csv(field_long, file.path(out_dir, "field_leg_blocked_panel_data.csv"), na = "")
readr::write_csv(gamma_df, file.path(out_dir, "decirc_gamma_panels_data.csv"), na = "")
readr::write_csv(axis_cells, file.path(out_dir, "axis_recovery_quadrants_data.csv"), na = "")
readr::write_csv(aim_plot_df, file.path(out_dir, "aimnet2_ceiling_lift_data.csv"), na = "")

manifest <- tibble(
  chart = names(chart_paths),
  path = unname(chart_paths),
  shows = c(
    "Combined mechanism x stratum recovery grid with per-cell metric annotation and law-grade glyphs.",
    "APBS-EFG versus MOPAC-Coulomb-EFG versus MOPAC-Mc blocked within-frame field-leg recovery.",
    paste0(gamma_sym, " coefficient recovery with emitted-SE CIs against the ", gamma_sym, " = 1 reference."),
    "Between-static versus within-dynamic T2 recovery, separating axis placement from joint model claims.",
    paste0("AIMNet2 ceiling ", "R", sup2, " and feature lift bands for standalone-versus-joint context.")
  ),
  metric_note = c(
    paste0(
      "Mixed metrics are annotated per cell: capstone random within-frame |T2| r, ",
      "MOPAC blocked within-frame |T2| r, capstone ", sigma_sym, " R", sup2,
      ", AIMNet2 T2 R", sup2, ", and aromatic-H de-circ |T2| r."
    ),
    paste0("r panel has Fisher-z CI derived from emitted r and N_eff; R", sup2, " panel has no emitted CI."),
    "Buckingham omitted because no dimensionless gamma is emitted.",
    "Standalone sidecar correlations only; not a joint model comparison.",
    "Ceiling model context only; not a physical-law recovery."
  )
)
readr::write_csv(manifest, file.path(out_dir, "chart_manifest.csv"), na = "")

input_manifest <- tibble(input = names(paths), path = unname(paths))
readr::write_csv(input_manifest, file.path(out_dir, "input_manifest.csv"), na = "")

cat("wrote charts under ", out_dir, "\n", sep = "")
for (path in chart_paths) {
  cat(path, "\n", sep = "")
}
