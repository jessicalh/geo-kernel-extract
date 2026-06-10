normal_crps <- function(y, mu, sigma) {
  sigma <- pmax(sigma, .Machine$double.eps)
  z <- (y - mu) / sigma
  sigma * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}

calibration_coverage <- function(predictions, levels = c(0.50, 0.80, 0.90, 0.95)) {
  pred <- data.table::as.data.table(predictions)
  if ("step" %in% names(pred)) {
    keep <- pred[["step"]] == "L6"
    pred <- pred[keep, ]
  }
  if (!nrow(pred)) return(list(coverage = data.table::data.table(), scores = data.table::data.table()))
  pred[, resid := y - pred]
  sigmas <- pred[, .(sigma = stats::sd(resid)), by = element]
  pred <- merge(pred, sigmas, by = "element", all.x = TRUE)
  cov <- data.table::rbindlist(lapply(levels, function(level) {
    z <- stats::qnorm((1 + level) / 2)
    pred[, .(
      nominal = level,
      empirical = mean(abs(y - pred) <= z * sigma),
      n = .N,
      mean_width = mean(2 * z * sigma)
    ), by = element]
  }))
  scores <- pred[, .(
    n = .N,
    sigma = unique(sigma)[1],
    crps = mean(normal_crps(y, pred, sigma)),
    log_score = mean(stats::dnorm(y, mean = pred, sd = sigma, log = TRUE)),
    r2 = r2_score(y, pred)
  ), by = element]
  list(coverage = cov, scores = scores)
}
