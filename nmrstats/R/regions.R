angle_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(circular::mean.circular(circular::circular(x, units = "radians", modulo = "2pi")))
}

estimate_regions <- function(train_geom, fold_id = "full", grid_n = 72L, hdr = c(0.50, 0.80, 0.95)) {
  if (!all(c("phi", "psi") %in% names(train_geom))) abort("phi/psi columns missing for region estimation")
  finite <- train_geom[is.finite(phi) & is.finite(psi), ]
  breaks <- seq(-pi, pi, length.out = grid_n + 1L)
  if (nrow(finite)) {
    finite[, phi_bin := pmax(1L, pmin(grid_n, findInterval(phi, breaks, all.inside = TRUE)))]
    finite[, psi_bin := pmax(1L, pmin(grid_n, findInterval(psi, breaks, all.inside = TRUE)))]
    dens <- finite[, .N, by = .(phi_bin, psi_bin)][order(-N)]
    dens[, mass := N / sum(N)]
    dens[, cum_mass := cumsum(mass)]
    thresholds <- data.table::rbindlist(lapply(hdr, function(p) {
      inside <- dens[cum_mass <= p, N]
      data.table::data.table(hdr_mass = p, min_bin_n = if (length(inside)) min(inside) else max(dens$N))
    }))
  } else {
    dens <- data.table::data.table(phi_bin = integer(), psi_bin = integer(), N = integer(), mass = numeric(), cum_mass = numeric())
    thresholds <- data.table::data.table(hdr_mass = hdr, min_bin_n = NA_integer_)
  }
  rama <- train_geom[, .(
    n = .N,
    phi_mean = angle_mean(phi),
    psi_mean = angle_mean(psi),
    phi_finite = sum(is.finite(phi)),
    psi_finite = sum(is.finite(psi))
  ), by = .(rama_region)]
  chi <- data.table::rbindlist(lapply(1:4, function(i) {
    chi_col <- sprintf("chi%d", i)
    ex_col <- sprintf("chi%d_exists", i)
    if (!all(c(chi_col, ex_col) %in% names(train_geom))) {
      return(data.table::data.table(chi = i, n = 0L, mean_angle = NA_real_, rotamer_minus = 0L, rotamer_plus = 0L, rotamer_trans = 0L))
    }
    x <- train_geom[get(ex_col) == 1 & is.finite(get(chi_col)), get(chi_col)]
    data.table::data.table(
      chi = i,
      n = length(x),
      mean_angle = angle_mean(x),
      rotamer_minus = sum(x < -pi / 3),
      rotamer_plus = sum(x > pi / 3),
      rotamer_trans = sum(x >= -pi / 3 & x <= pi / 3)
    )
  }))
  list(
    schema_version = 1L,
    generated_by = "nmrstats::estimate_regions",
    training_split_id = fold_id,
    row_count = nrow(train_geom),
    finite_phi_psi = nrow(finite),
    angle_convention = "IUPAC signed radians",
    grid = list(n = grid_n, min = -pi, max = pi),
    hdr_thresholds = thresholds,
    rama_summary = rama,
    chi_summary = chi,
    note = "Core pass estimates regions only; region_def_id remains emitted placeholder until C++ re-emit round-trip."
  )
}

with_fold_each_regions <- function(rows, folds) {
  lapply(folds, function(f) estimate_regions(rows[f$train, ], fold_id = f$id))
}

freeze_regions <- function(regions, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(regions, path, auto_unbox = TRUE, pretty = TRUE, dataframe = "rows", na = "null")
  invisible(path)
}

region_summary <- function(regions) {
  data.table::rbindlist(lapply(regions, function(r) {
    data.table::data.table(
      training_split_id = r$training_split_id,
      row_count = r$row_count,
      finite_phi_psi = r$finite_phi_psi,
      finite_fraction = r$finite_phi_psi / r$row_count
    )
  }))
}
