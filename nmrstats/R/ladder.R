LADDER <- data.table::data.table(
  step = c("L0", "L1", "L2", "L3", "L4", "L5", "L6"),
  description = c(
    "in-fold posterior identity intercept baseline",
    "baseline + dataset/iupac role",
    "L1 + identity/protein shrinkage",
    "L2 + phi/psi, chi, sequence, secondary structure",
    "L3 + spatial propinquity and emitted region placeholder",
    "L4 + available proxy blocks",
    "L5 + dataset-by-proxy interactions where estimable"
  )
)

r2_score <- function(y, pred) {
  ok <- is.finite(y) & is.finite(pred)
  if (sum(ok) < 2L) return(NA_real_)
  denom <- sum((y[ok] - mean(y[ok]))^2)
  if (denom <= 0) return(NA_real_)
  1 - sum((y[ok] - pred[ok])^2) / denom
}

residual_icc <- function(resid, group) {
  ok <- is.finite(resid) & !is.na(group)
  if (sum(ok) < 3L || data.table::uniqueN(group[ok]) < 2L) return(NA_real_)
  dt <- data.table::data.table(resid = resid[ok], group = group[ok])
  means <- dt[, .(m = mean(resid), n = .N, v = stats::var(resid)), by = group]
  between <- stats::var(means$m)
  within <- stats::weighted.mean(means$v[is.finite(means$v)], means$n[is.finite(means$v)])
  if (!is.finite(between) || !is.finite(within) || between + within <= 0) return(NA_real_)
  between / (between + within)
}

metric_table <- function(pred_dt) {
  pred_dt[, .(
    n = .N,
    r2 = r2_score(y, pred),
    rmse = sqrt(mean((y - pred)^2)),
    mae = mean(abs(y - pred)),
    residual_icc_protein = residual_icc(y - pred, protein_id)
  ), by = .(step, element)]
}

fit_predict_step <- function(train_e, test_e, step, nthreads) {
  if (identical(step, "L0")) return(rep(0, nrow(test_e)))
  if ("protein_id" %in% names(train_e) && "protein_id" %in% names(test_e)) {
    train_proteins <- unique(as.character(train_e$protein_id))
    test_proteins <- unique(as.character(test_e$protein_id))
    if (length(setdiff(test_proteins, train_proteins))) {
      train_e <- data.table::copy(train_e)
      test_e <- data.table::copy(test_e)
      train_e[, protein_id := NULL]
      test_e[, protein_id := NULL]
    }
  }
  blocks <- feature_blocks(train_e)
  f <- formula_for_step(
    train_e,
    response = "target_work",
    step = step,
    blocks = blocks,
    include_dataset_interactions = identical(step, "L6")
  )
  if (!length(attr(stats::terms(f), "term.labels"))) {
    return(rep(mean(train_e$target_work), nrow(test_e)))
  }
  fit <- fit_bam_model(train_e, f, family = stats::gaussian(), nthreads = nthreads)
  as.numeric(predict_mu(fit, test_e))
}

run_ladder <- function(rows, folds, nthreads = 1L, steps = LADDER$step, keep_step_predictions = "L6") {
  dt <- prepare_model_data(rows)
  dt[, .row_id := seq_len(.N)]
  pred_parts <- list()
  metric_parts <- list()
  timings <- list()
  for (fold_i in seq_along(folds)) {
    fold <- folds[[fold_i]]
    message(sprintf("run_ladder: %s train=%d test=%d", fold$id, length(fold$train), length(fold$test)))
    train <- dt[fold$train, ]
    test <- dt[fold$test, ]
    base <- baseline_blup(train)
    train[, baseline_pred := predict(base, train)]
    test[, baseline_pred := predict(base, test)]
    train[, target_work := target_T0 - baseline_pred]
    elements <- sort(unique(as.character(test$element)))
    for (step in steps) {
      step_preds <- vector("list", length(elements))
      names(step_preds) <- elements
      started <- Sys.time()
      for (el in elements) {
        train_e <- train[as.character(element) == el, ]
        test_e <- test[as.character(element) == el, ]
        if (!nrow(train_e) || !nrow(test_e)) next
        pr <- fit_predict_step(train_e, test_e, step = step, nthreads = nthreads)
        step_preds[[el]] <- data.table::data.table(
          row_id = test_e$.row_id,
          fold = fold$id,
          step = step,
          element = as.character(test_e$element),
          dataset_id = as.character(test_e$dataset_id),
          protein_id = as.character(test_e$protein_id),
          y = test_e$target_T0,
          baseline = test_e$baseline_pred,
          pred = test_e$baseline_pred + pr
        )
      }
      pred_dt <- data.table::rbindlist(step_preds, fill = TRUE)
      metric_parts[[paste(fold$id, step, sep = "_")]] <- metric_table(pred_dt)[, fold := fold$id][]
      timings[[paste(fold$id, step, sep = "_")]] <- data.table::data.table(
        fold = fold$id, step = step,
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
      if (step %in% keep_step_predictions) {
        pred_parts[[paste(fold$id, step, sep = "_")]] <- pred_dt
      }
      rm(pred_dt)
      gc()
    }
  }
  metrics <- data.table::rbindlist(metric_parts, fill = TRUE)
  avg <- metrics[, .(
    n = sum(n),
    r2 = stats::weighted.mean(r2, n, na.rm = TRUE),
    rmse = stats::weighted.mean(rmse, n, na.rm = TRUE),
    mae = stats::weighted.mean(mae, n, na.rm = TRUE),
    residual_icc_protein = stats::weighted.mean(residual_icc_protein, n, na.rm = TRUE)
  ), by = .(step, element)]
  data.table::setorder(avg, element, step)
  avg[, delta_r2 := r2 - data.table::shift(r2), by = element]
  list(
    by_fold = metrics,
    average = avg,
    predictions = data.table::rbindlist(pred_parts, fill = TRUE),
    timings = data.table::rbindlist(timings, fill = TRUE),
    fold_summary = fold_summary(rows, folds)
  )
}
