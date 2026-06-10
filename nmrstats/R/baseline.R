make_key <- function(dt, cols) {
  if (!length(cols)) return(rep("all", nrow(dt)))
  do.call(paste, c(lapply(cols, function(z) as.character(dt[[z]])), sep = "\r"))
}

shrink_map <- function(dt, cols, residual, lambda = 20) {
  group_key <- make_key(dt, cols)
  tab <- data.table::data.table(group_key = group_key, residual = residual)[, .(n = .N, mean = mean(residual)), by = group_key]
  tab[, value := (n * mean) / (n + lambda)]
  list(cols = cols, table = tab[, .(group_key, value)])
}

apply_shrink_map <- function(map, dt) {
  group_key <- make_key(dt, map$cols)
  out <- map$table$value[match(group_key, map$table$group_key)]
  out[is.na(out)] <- 0
  out
}

baseline_blup <- function(train, y = "target_T0", lambda = 20) {
  if (!y %in% names(train)) abort("baseline target %s missing", y)
  global <- mean(train[[y]])
  pred <- rep(global, nrow(train))
  specs <- list(
    c("element"),
    c("element", "iupac_role"),
    c("element", "residue_type"),
    c("element", "equivalence_class")
  )
  specs <- Filter(function(cols) all(cols %in% names(train)), specs)
  maps <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    m <- shrink_map(train, specs[[i]], train[[y]] - pred, lambda = lambda)
    maps[[i]] <- m
    pred <- pred + apply_shrink_map(m, train)
  }
  structure(
    list(global = global, maps = maps, y = y, train_keys = train$.row_id %||% seq_len(nrow(train))),
    class = "nmrstats_baseline"
  )
}

predict.nmrstats_baseline <- function(object, newdata, ...) {
  pred <- rep(object$global, nrow(newdata))
  for (m in object$maps) pred <- pred + apply_shrink_map(m, newdata)
  pred
}

resid_target <- function(rows, fold, y = "target_T0", lambda = 20) {
  train <- rows[fold$train, ]
  fit <- baseline_blup(train, y = y, lambda = lambda)
  pred <- predict(fit, rows[fold$test, ])
  data.table::data.table(
    row_id = fold$test,
    baseline = pred,
    resid = rows[[y]][fold$test] - pred
  )
}

with_fold <- function(rows, fold, fn, ...) {
  train <- rows[fold$train, ]
  test_ids <- fold$test
  if (".row_id" %in% names(train) && any(train$.row_id %in% test_ids)) {
    abort("fold leakage guard failed: training rows include held-out row ids")
  }
  fn(train, ...)
}
