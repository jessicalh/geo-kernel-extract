joint_views <- function(predictions) {
  pred <- data.table::as.data.table(predictions)
  if ("step" %in% names(pred)) {
    keep <- pred[["step"]] == "L6"
    pred <- pred[keep, ]
  }
  if (!nrow(pred)) return(data.table::data.table())
  by_dataset <- pred[, .(
    n = .N,
    r2 = r2_score(y, pred),
    rmse = sqrt(mean((y - pred)^2)),
    mae = mean(abs(y - pred))
  ), by = .(view = dataset_id, element)]
  combined <- pred[, .(
    n = .N,
    r2 = r2_score(y, pred),
    rmse = sqrt(mean((y - pred)^2)),
    mae = mean(abs(y - pred))
  ), by = .(element)]
  combined[, view := "combined_with_dataset_effect"]
  data.table::rbindlist(list(by_dataset, combined), fill = TRUE)
}

fit_joint <- function(rows, nthreads = 1L) {
  dt <- prepare_model_data(rows)
  out <- list()
  for (el in sort(unique(dt$element))) {
    message(sprintf("fit_joint: element %s", el))
    d <- dt[element == el, ]
    blocks <- feature_blocks(d)
    f <- formula_for_step(d, response = "target_T0", step = "L6", blocks = blocks, include_dataset_interactions = TRUE)
    started <- Sys.time()
    fit <- fit_bam_model(d, f, family = stats::gaussian(), nthreads = nthreads)
    s <- summary(fit)
    out[[as.character(el)]] <- data.table::data.table(
      element = el,
      n = nrow(d),
      elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      dev_expl = s$dev.expl %||% NA_real_,
      r_sq = s$r.sq %||% NA_real_,
      edf = sum(fit$edf),
      formula = paste(deparse(f), collapse = " ")
    )
    rm(fit)
    gc()
  }
  data.table::rbindlist(out, fill = TRUE)
}
