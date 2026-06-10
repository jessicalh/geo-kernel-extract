permute_cols_by <- function(dt, cols, by, seed = 1L) {
  cols <- intersect(cols, names(dt))
  if (!length(cols)) return(dt)
  set.seed(seed)
  by <- intersect(by, names(dt))
  if (!length(by)) by <- "element"
  dt[, (cols) := lapply(.SD, function(x) sample(x, length(x))), by = by, .SDcols = cols]
  dt
}

shuffle_target_work <- function(train, seed = 1L) {
  set.seed(seed)
  is_traj <- ("time_ps" %in% names(train) && any(is.finite(train$time_ps))) |
    ("pose_kind" %in% names(train) && any(grepl("trajectory", train$pose_kind)))
  if (is_traj && "atom_uid" %in% names(train)) {
    train[, target_work := sample(target_work, .N), by = atom_uid]
  } else {
    train[, target_work := sample(target_work, .N), by = .(dataset_id, element)]
  }
  train
}

score_oof_scenario <- function(rows, folds, scenario, transform_all = NULL, transform_train = NULL,
                               blocks_fun = feature_blocks, step = "L5", nthreads = 1L, seed = 1L) {
  dt <- prepare_model_data(rows)
  dt[, .row_id := seq_len(.N)]
  if (!is.null(transform_all)) dt <- transform_all(dt, seed = seed)
  out <- list()
  for (fold in folds) {
    message(sprintf("run_nulls[%s]: %s", scenario, fold$id))
    train <- dt[fold$train, ]
    test <- dt[fold$test, ]
    base <- baseline_blup(train)
    train[, baseline_pred := predict(base, train)]
    test[, baseline_pred := predict(base, test)]
    train[, target_work := target_T0 - baseline_pred]
    if (!is.null(transform_train)) train <- transform_train(train, seed = seed + match(fold$id, vapply(folds, `[[`, "", "id")))
    for (el in sort(unique(test$element))) {
      train_e <- train[element == el, ]
      test_e <- test[element == el, ]
      if (!nrow(train_e) || !nrow(test_e)) next
      blocks <- blocks_fun(train_e)
      f <- formula_for_step(train_e, response = "target_work", step = step, blocks = blocks)
      fit <- fit_bam_model(train_e, f, family = stats::gaussian(), nthreads = nthreads)
      pred <- test_e$baseline_pred + predict_mu(fit, test_e)
      out[[length(out) + 1L]] <- data.table::data.table(
        scenario = scenario,
        fold = fold$id,
        element = el,
        n = nrow(test_e),
        r2 = r2_score(test_e$target_T0, pred),
        rmse = sqrt(mean((test_e$target_T0 - pred)^2)),
        residual_icc_protein = residual_icc(test_e$target_T0 - pred, test_e$protein_id)
      )
    }
  }
  res <- data.table::rbindlist(out, fill = TRUE)
  avg <- res[, .(
    n = sum(n),
    r2 = stats::weighted.mean(r2, n, na.rm = TRUE),
    rmse = stats::weighted.mean(rmse, n, na.rm = TRUE),
    residual_icc_protein = stats::weighted.mean(residual_icc_protein, n, na.rm = TRUE)
  ), by = .(scenario, element)]
  list(by_fold = res, average = avg)
}

legal_sign_columns <- function(rows) {
  ir <- irrep_columns(rows)
  if (!nrow(ir) || !"sign_flip_legal" %in% names(ir)) return(character())
  flag <- as.logical(ir$sign_flip_legal)
  unique(ir[!is.na(flag) & flag, name])
}

run_sign_control <- function(rows, folds, nthreads = 1L) {
  cols <- intersect(legal_sign_columns(rows), proxy_columns(rows))
  if (!length(cols)) {
    return(list(
      by_fold = data.table::data.table(),
      average = data.table::data.table(
        scenario = "sign_flip",
        element = sort(unique(rows$element)),
        n = NA_integer_,
        r2 = NA_real_,
        rmse = NA_real_,
        residual_icc_protein = NA_real_,
        note = "no sign_flip_legal proxy columns declared by column_irrep_schema/null_spec"
      )
    ))
  }
  score_oof_scenario(
    rows, folds, scenario = "sign_flip",
    transform_all = function(dt, seed) {
      dt[, (cols) := lapply(.SD, function(x) -x), .SDcols = cols]
      dt
    },
    nthreads = nthreads
  )
}

run_nulls <- function(rows, folds, nthreads = 1L, seed = 1L) {
  proxy <- proxy_columns(rows)
  proxy <- proxy[!grepl("bare_efg_kernel", proxy)]
  list(
    target_shuffle = score_oof_scenario(
      rows, folds, scenario = "target_shuffle",
      transform_train = shuffle_target_work,
      nthreads = nthreads, seed = seed
    ),
    source_identity_shuffle = score_oof_scenario(
      rows, folds, scenario = "source_identity_shuffle",
      transform_all = function(dt, seed) permute_cols_by(dt, proxy, c("element", "residue_type", "rama_region"), seed = seed),
      nthreads = nthreads, seed = seed + 100L
    ),
    scrambled_physics = score_oof_scenario(
      rows, folds, scenario = "scrambled_physics",
      transform_all = function(dt, seed) {
        ring <- intersect(feature_blocks(dt)$ring %||% character(), names(dt))
        field <- intersect(feature_blocks(dt)$field_charge %||% character(), names(dt))
        dt <- permute_cols_by(dt, ring, c("element", "dataset_id"), seed = seed)
        dt <- permute_cols_by(dt, field, c("element", "dataset_id"), seed = seed + 1L)
        dt
      },
      nthreads = nthreads, seed = seed + 200L
    ),
    self_bonded_excluded = score_oof_scenario(
      rows, folds, scenario = "self_bonded_excluded",
      transform_all = function(dt, seed) {
        cols <- intersect(c("self_or_bonded_atom_count", "self_or_bonded_bond_count"), names(dt))
        if (length(cols)) dt[, (cols) := 0]
        dt
      },
      nthreads = nthreads, seed = seed + 300L
    ),
    sign_controls = run_sign_control(rows, folds, nthreads = nthreads)
  )
}

null_average_table <- function(nulls) {
  data.table::rbindlist(lapply(nulls, function(x) x$average), fill = TRUE)
}
