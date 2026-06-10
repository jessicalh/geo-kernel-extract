factor_term <- function(dt, col) {
  if (!col %in% names(dt)) return(character())
  if (!is.factor(dt[[col]])) return(character())
  if (data.table::uniqueN(as.character(dt[[col]])) < 2L) return(character())
  col
}

re_term <- function(dt, col, max_levels = 5000L) {
  if (!col %in% names(dt)) return(character())
  n <- if (is.factor(dt[[col]])) data.table::uniqueN(as.character(dt[[col]])) else data.table::uniqueN(dt[[col]])
  if (n < 2L || n > max_levels) return(character())
  sprintf("s(%s, bs='re')", col)
}

linear_numeric_term <- function(dt, col) {
  if (!col %in% names(dt)) return(character())
  x <- dt[[col]]
  if (!is.numeric(x) && !is.integer(x)) return(character())
  if (data.table::uniqueN(x[is.finite(x)]) < 2L) return(character())
  col
}

smooth_numeric_term <- function(dt, col, k = 6L) {
  if (!col %in% names(dt)) return(character())
  x <- dt[[col]]
  if (!is.numeric(x) && !is.integer(x)) return(character())
  if (data.table::uniqueN(x[is.finite(x)]) < k) return(linear_numeric_term(dt, col))
  sprintf("s(%s, k=%d)", col, k)
}

element_k <- function(n) {
  if (n >= 100000L) list(k_phi = 8L, k_chi = 6L)
  else if (n >= 25000L) list(k_phi = 6L, k_chi = 5L)
  else list(k_phi = 5L, k_chi = 4L)
}

model_knots <- function() {
  list(
    phi_model = c(-pi, pi), psi_model = c(-pi, pi),
    chi1_model = c(-pi, pi), chi2_model = c(-pi, pi),
    chi3_model = c(-pi, pi), chi4_model = c(-pi, pi)
  )
}

build_model_terms <- function(dt, step = "L6", blocks = feature_blocks(dt), include_dataset_interactions = FALSE) {
  k <- element_k(nrow(dt))
  terms <- character()
  if (step %in% c("L1", "L2", "L3", "L4", "L5", "L6")) {
    terms <- c(terms, factor_term(dt, "dataset_id"), factor_term(dt, "iupac_role"))
  }
  if (step %in% c("L2", "L3", "L4", "L5", "L6")) {
    terms <- c(
      terms,
      re_term(dt, "residue_type"),
      re_term(dt, "equivalence_class"),
      re_term(dt, "protein_id", max_levels = 1500L)
    )
  }
  if (step %in% c("L3", "L4", "L5", "L6")) {
    if (all(c("phi_model", "psi_model", "phi_present") %in% names(dt)) && any(dt$phi_present == 1L)) {
      terms <- c(terms, sprintf("te(phi_model, psi_model, bs=c('cc','cc'), k=c(%d,%d), by=phi_present)", k$k_phi, k$k_phi), "phi_present")
    }
    for (i in 1:4) {
      chi <- sprintf("chi%d_model", i)
      ex <- sprintf("chi%d_exists", i)
      if (all(c(chi, ex) %in% names(dt)) && any(dt[[ex]] == 1L)) {
        terms <- c(terms, sprintf("s(%s, bs='cc', k=%d, by=%s)", chi, k$k_chi, ex), ex)
      }
    }
    terms <- c(
      terms,
      factor_term(dt, "prev_class"), factor_term(dt, "next_class"),
      factor_term(dt, "pre_proline"), factor_term(dt, "dssp_ss3"),
      factor_term(dt, "rama_region")
    )
  }
  if (step %in% c("L4", "L5", "L6")) {
    spatial <- spatial_features(dt)
    smooth_spatial <- intersect(c("nearest_ring_dist", "nearest_charge_dist", "nearest_bond_dist", "nearest_atom_dist", "ring_cyl_z", "ring_cyl_rho", "ring_angle_to_normal"), spatial)
    linear_spatial <- setdiff(spatial, smooth_spatial)
    terms <- c(terms, unlist(lapply(smooth_spatial, function(z) smooth_numeric_term(dt, z, k = 6L)), use.names = FALSE))
    terms <- c(terms, unlist(lapply(linear_spatial, function(z) linear_numeric_term(dt, z)), use.names = FALSE))
    terms <- c(terms, factor_term(dt, "region_def_id"))
  }
  if (step %in% c("L5", "L6")) {
    proxy <- unique(unlist(blocks, use.names = FALSE))
    proxy <- proxy[!grepl("bare_efg_kernel", proxy)]
    terms <- c(terms, unlist(lapply(proxy, function(z) linear_numeric_term(dt, z)), use.names = FALSE))
  }
  if (identical(step, "L6") && include_dataset_interactions && nlevels(dt$dataset_id) > 1L) {
    proxy <- unique(unlist(blocks, use.names = FALSE))
    proxy <- proxy[!grepl("bare_efg_kernel", proxy)]
    proxy <- proxy[vapply(proxy, function(z) length(linear_numeric_term(dt, z)) > 0, logical(1))]
    terms <- c(terms, sprintf("dataset_id:%s", proxy))
  }
  unique(terms[nzchar(terms)])
}

formula_for_step <- function(dt, response = "target_T0", step = "L6", blocks = feature_blocks(dt), include_dataset_interactions = FALSE) {
  terms <- build_model_terms(dt, step = step, blocks = blocks, include_dataset_interactions = include_dataset_interactions)
  stats::as.formula(paste(response, "~", if (length(terms)) paste(terms, collapse = " + ") else "1"))
}

fit_bam_model <- function(dt, formula, family = stats::gaussian(), nthreads = 1L) {
  mgcv::bam(
    formula,
    data = dt,
    family = family,
    method = "fREML",
    discrete = TRUE,
    nthreads = nthreads,
    knots = model_knots(),
    na.action = stats::na.fail
  )
}

predict_mu <- function(fit, newdata) {
  if (inherits(fit, "nmrstats_hetero_bam")) {
    return(as.numeric(stats::predict(fit$mean_fit, newdata = newdata, type = "response")))
  }
  p <- stats::predict(fit, newdata = newdata, type = "response")
  if (is.matrix(p)) p[, 1L] else as.numeric(p)
}

fit_element_two_stage <- function(dt, mu, sigma, nthreads, started, reason) {
  mean_fit <- fit_bam_model(dt, mu, family = stats::gaussian(), nthreads = nthreads)
  resid <- dt$target_T0 - as.numeric(stats::predict(mean_fit, newdata = dt, type = "response"))
  floor_var <- stats::quantile(resid^2, probs = 0.05, na.rm = TRUE)
  if (!is.finite(floor_var) || floor_var <= 0) floor_var <- .Machine$double.eps
  dt[, log_sigma_work := 0.5 * log(pmax(resid^2, floor_var))]
  sigma_rhs <- paste(deparse(sigma), collapse = " ")
  sigma_rhs <- sub("^~", "", sigma_rhs)
  scale_formula <- stats::as.formula(paste("log_sigma_work ~", sigma_rhs))
  scale_fit <- fit_bam_model(dt, scale_formula, family = stats::gaussian(), nthreads = nthreads)
  fit <- list(mean_fit = mean_fit, scale_fit = scale_fit, fallback_reason = reason)
  class(fit) <- "nmrstats_hetero_bam"
  attr(fit, "nmrstats_elapsed_sec") <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  attr(fit, "nmrstats_formula") <- list(mu = deparse(mu), sigma = deparse(scale_formula))
  fit
}

fit_element <- function(train_e, element = unique(train_e$element), nthreads = 1L, include_dataset_interactions = TRUE) {
  dt <- prepare_model_data(train_e)
  blocks <- feature_blocks(dt)
  mu <- formula_for_step(dt, response = "target_T0", step = "L6", blocks = blocks, include_dataset_interactions = include_dataset_interactions)
  sigma_terms <- c(factor_term(dt, "dataset_id"), factor_term(dt, "iupac_role"))
  sigma <- stats::as.formula(paste("~", if (length(sigma_terms)) paste(sigma_terms, collapse = " + ") else "1"))
  started <- Sys.time()
  fit <- tryCatch(
    mgcv::bam(
      list(mu, sigma),
      data = dt,
      family = mgcv::gaulss(),
      method = "fREML",
      discrete = TRUE,
      nthreads = nthreads,
      knots = model_knots(),
      na.action = stats::na.fail
    ),
    error = function(e) {
      if (grepl("general families not supported by bam", conditionMessage(e), fixed = TRUE)) {
        fit_element_two_stage(dt, mu, sigma, nthreads, started, conditionMessage(e))
      } else {
        stop(e)
      }
    }
  )
  if (!inherits(fit, "nmrstats_hetero_bam")) {
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    attr(fit, "nmrstats_elapsed_sec") <- elapsed
    attr(fit, "nmrstats_formula") <- list(mu = deparse(mu), sigma = deparse(sigma))
  }
  attr(fit, "nmrstats_element") <- element
  fit
}

summarize_fit <- function(fit, element, n) {
  if (inherits(fit, "nmrstats_hetero_bam")) {
    sm <- summary(fit$mean_fit)
    ss <- summary(fit$scale_fit)
    return(data.table::data.table(
      element = as.character(element),
      n = n,
      family = "gaussian_bam_mean_plus_log_resid_scale_bam",
      converged = TRUE,
      elapsed_sec = attr(fit, "nmrstats_elapsed_sec") %||% NA_real_,
      dev_expl = sm$dev.expl %||% NA_real_,
      r_sq = sm$r.sq %||% NA_real_,
      scale_dev_expl = ss$dev.expl %||% NA_real_,
      scale_r_sq = ss$r.sq %||% NA_real_,
      edf = sum(fit$mean_fit$edf) + sum(fit$scale_fit$edf),
      mu_formula = paste(attr(fit, "nmrstats_formula")$mu, collapse = " "),
      sigma_formula = paste(attr(fit, "nmrstats_formula")$sigma, collapse = " "),
      adjustment = fit$fallback_reason
    ))
  }
  s <- summary(fit)
  conv <- TRUE
  if (!is.null(fit$outer.info$conv)) conv <- isTRUE(fit$outer.info$conv == "full convergence")
  data.table::data.table(
    element = as.character(element),
    n = n,
    family = fit$family$family,
    converged = conv,
    elapsed_sec = attr(fit, "nmrstats_elapsed_sec") %||% NA_real_,
    dev_expl = s$dev.expl %||% NA_real_,
    r_sq = s$r.sq %||% NA_real_,
    edf = sum(fit$edf),
    mu_formula = paste(attr(fit, "nmrstats_formula")$mu, collapse = " "),
    sigma_formula = paste(attr(fit, "nmrstats_formula")$sigma, collapse = " ")
  )
}

fit_per_element <- function(rows, nthreads = 1L, save_dir = NULL) {
  elements <- sort(unique(rows$element))
  out <- vector("list", length(elements))
  names(out) <- elements
  for (el in elements) {
    message(sprintf("fit_per_element: element %s (%d rows)", el, sum(rows$element == el)))
    fit <- fit_element(rows[element == el, ], element = el, nthreads = nthreads)
    out[[as.character(el)]] <- summarize_fit(fit, el, sum(rows$element == el))
    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(fit, file.path(save_dir, sprintf("element_%s_gaulss.rds", el)))
    }
    rm(fit)
    gc()
  }
  data.table::rbindlist(out, fill = TRUE)
}
