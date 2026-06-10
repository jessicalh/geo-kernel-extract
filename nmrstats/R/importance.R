proxy_columns <- function(rows) {
  unique(unlist(feature_blocks(rows), use.names = FALSE))
}

linear_relaimpo <- function(rows) {
  dt <- prepare_model_data(rows)
  cols <- proxy_columns(dt)
  cols <- cols[!grepl("bare_efg_kernel", cols)]
  out <- list()
  for (el in sort(unique(dt$element))) {
    d <- dt[element == el, ]
    xcols <- cols[vapply(cols, function(z) z %in% names(d) && stats::sd(d[[z]]) > 0, logical(1))]
    if (length(xcols) < 2L) {
      out[[as.character(el)]] <- data.table::data.table(element = el, proxy = xcols, lmg = NA_real_, note = "fewer_than_two_proxy_columns")
      next
    }
    base <- baseline_blup(d)
    d[, target_work := target_T0 - predict(base, d)]
    f <- stats::as.formula(paste("target_work ~", paste(xcols, collapse = " + ")))
    fit <- stats::lm(f, data = d)
    ri <- relaimpo::calc.relimp(fit, type = "lmg", rela = TRUE)
    vals <- ri$lmg
    out[[as.character(el)]] <- data.table::data.table(
      element = el,
      proxy = names(vals),
      lmg = as.numeric(vals),
      note = "full-data linear residual diagnostic"
    )
  }
  data.table::rbindlist(out, fill = TRUE)
}

fit_predict_block_model <- function(train_e, test_e, blocks, nthreads) {
  f <- formula_for_step(train_e, response = "target_work", step = "L5", blocks = blocks)
  fit <- fit_bam_model(train_e, f, family = stats::gaussian(), nthreads = nthreads)
  as.numeric(predict_mu(fit, test_e))
}

drop_one_block_oof <- function(rows, folds, nthreads = 1L) {
  dt <- prepare_model_data(rows)
  dt[, .row_id := seq_len(.N)]
  out <- list()
  for (fold in folds) {
    message(sprintf("drop_one_block_oof: %s", fold$id))
    train <- dt[fold$train, ]
    test <- dt[fold$test, ]
    base <- baseline_blup(train)
    train[, baseline_pred := predict(base, train)]
    test[, baseline_pred := predict(base, test)]
    train[, target_work := target_T0 - baseline_pred]
    for (el in sort(unique(test$element))) {
      train_e <- train[element == el, ]
      test_e <- test[element == el, ]
      if (!nrow(train_e) || !nrow(test_e)) next
      blocks <- feature_blocks(train_e)
      if (!length(blocks)) next
      full_resid <- fit_predict_block_model(train_e, test_e, blocks, nthreads)
      y <- test_e$target_T0
      full_pred <- test_e$baseline_pred + full_resid
      full_r2 <- r2_score(y, full_pred)
      base_blocks <- list()
      base_r2 <- r2_score(y, test_e$baseline_pred + fit_predict_block_model(train_e, test_e, base_blocks, nthreads))
      out[[length(out) + 1L]] <- data.table::data.table(
        fold = fold$id, element = el, block = "all_proxy_vs_L4",
        r2_full = full_r2, r2_without = base_r2, delta_r2 = full_r2 - base_r2, n = nrow(test_e)
      )
      for (block in names(blocks)) {
        reduced <- blocks[names(blocks) != block]
        reduced_resid <- fit_predict_block_model(train_e, test_e, reduced, nthreads)
        reduced_r2 <- r2_score(y, test_e$baseline_pred + reduced_resid)
        out[[length(out) + 1L]] <- data.table::data.table(
          fold = fold$id, element = el, block = block,
          r2_full = full_r2, r2_without = reduced_r2, delta_r2 = full_r2 - reduced_r2, n = nrow(test_e)
        )
      }
    }
  }
  res <- data.table::rbindlist(out, fill = TRUE)
  avg <- res[, .(
    n = sum(n),
    r2_full = stats::weighted.mean(r2_full, n, na.rm = TRUE),
    r2_without = stats::weighted.mean(r2_without, n, na.rm = TRUE),
    delta_r2 = stats::weighted.mean(delta_r2, n, na.rm = TRUE)
  ), by = .(element, block)]
  list(by_fold = res, average = avg)
}

proxy_correspondence <- function(rows, predictions = NULL) {
  dt <- prepare_model_data(rows)
  cols_by_block <- feature_blocks(dt)
  if (!length(cols_by_block)) return(data.table::data.table())
  if (!is.null(predictions) && nrow(predictions)) {
    predictions <- data.table::as.data.table(predictions)
    err <- predictions[predictions[["step"]] == "L6", .(row_id, abs_resid = abs(y - pred))]
    dt[, row_id := seq_len(.N)]
    dt <- merge(dt, err, by = "row_id", all.x = FALSE, all.y = FALSE)
  } else {
    dt[, abs_resid := abs(target_T0 - mean(target_T0)), by = element]
  }
  region <- if ("rama_region" %in% names(dt)) "rama_region" else "region_def_id"
  out <- list()
  for (block in names(cols_by_block)) {
    cols <- cols_by_block[[block]]
    cols <- cols[!grepl("bare_efg_kernel", cols)]
    if (!length(cols)) next
    dt[, proxy_strength_tmp := rowMeans(abs(.SD), na.rm = TRUE), .SDcols = cols]
    tab <- dt[, .(
      proxy_strength = mean(proxy_strength_tmp, na.rm = TRUE),
      abs_resid = mean(abs_resid, na.rm = TRUE),
      n = .N
    ), by = c("element", region)]
    corr <- tab[, .(
      correspondence = suppressWarnings(stats::cor(proxy_strength, abs_resid, use = "complete.obs")),
      n_regions = .N
    ), by = element]
    corr[, block := block]
    out[[block]] <- corr
    dt[, proxy_strength_tmp := NULL]
  }
  data.table::rbindlist(out, fill = TRUE)
}

importance <- function(rows, folds, predictions = NULL, nthreads = 1L) {
  list(
    relaimpo = linear_relaimpo(rows),
    drop_one_block = drop_one_block_oof(rows, folds, nthreads = nthreads),
    correspondence = proxy_correspondence(rows, predictions)
  )
}
