# common.R — Shared theme, palettes, and helpers for secondary analysis.
#
# Source this before any analysis-specific R script:
#   source("common.R")

library(ggplot2)
library(patchwork)
library(scales)
library(readr)

# ── Paths ────────────────────────────────────────────────────────────
DATA_DIR <- file.path(dirname(getwd()), "src", "output", "secondary")
FIG_DIR  <- file.path(DATA_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Allow override from caller
set_data_dir <- function(path) {
  DATA_DIR <<- path
  FIG_DIR  <<- file.path(path, "figures")
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ── Ring type palette (colourblind-safe, 8 types) ────────────────────
RING_COLOURS <- c(
  PHE           = "#E69F00",  # orange
  TYR           = "#56B4E9",  # sky blue
  TRP_benzene   = "#009E73",  # bluish green
  TRP_pyrrole   = "#0072B2",  # blue
  TRP_perimeter = "#CC79A7",  # reddish purple
  HIS           = "#D55E00",  # vermillion
  HID           = "#F0E442",  # yellow
  HIE           = "#999999",  # grey
  none          = "#CCCCCC"   # light grey
)

ring_colour_scale <- function(...) {
  scale_colour_manual(values = RING_COLOURS, ...)
}

ring_fill_scale <- function(...) {
  scale_fill_manual(values = RING_COLOURS, ...)
}

# ── Stratum palette ──────────────────────────────────────────────────
STRATUM_COLOURS <- c(
  hie_only = "#999999",
  phe_only = "#E69F00",
  tyr_only = "#56B4E9",
  trp_only = "#009E73",
  no_hie   = "#0072B2",
  all      = "#333333"
)

stratum_fill_scale <- function(...) {
  scale_fill_manual(values = STRATUM_COLOURS, ...)
}

# ── Element palette ──────────────────────────────────────────────────
ELEMENT_COLOURS <- c(
  H = "#E69F00",
  C = "#333333",
  N = "#0072B2",
  O = "#D55E00",
  S = "#F0E442"
)

# ── Thesis theme ─────────────────────────────────────────────────────
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "grey70"),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text        = element_text(face = "bold", size = rel(0.9)),
      legend.position   = "bottom",
      plot.title        = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle     = element_text(colour = "grey40")
    )
}

# Set as default
theme_set(theme_thesis())

# ── Helpers ──────────────────────────────────────────────────────────
save_fig <- function(p, name, width = 8, height = 6) {
  path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  ggsave(path, p, width = width, height = height, device = cairo_pdf)
  message("Saved: ", path)
  invisible(path)
}
