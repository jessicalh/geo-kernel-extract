n_eff_ar1 <- function(x, block = NULL) {
  ok <- is.finite(x)
  if (!is.null(block)) ok <- ok & !is.na(block)
  x <- x[ok]
  if (!is.null(block)) block <- block[ok]
  if (length(x) < 3) return(length(x))
  if (is.null(block)) {
    rho <- suppressWarnings(stats::cor(utils::head(x, -1), utils::tail(x, -1)))
    n <- length(x)
  } else {
    pieces <- split(x, block)
    pairs <- lapply(pieces, function(z) {
      if (length(z) < 3) return(NULL)
      data.frame(a = utils::head(z, -1), b = utils::tail(z, -1))
    })
    pairs <- do.call(rbind, pairs)
    if (is.null(pairs) || nrow(pairs) < 3) return(length(x))
    rho <- suppressWarnings(stats::cor(pairs$a, pairs$b))
    n <- length(x)
  }
  if (!is.finite(rho)) rho <- 0
  rho <- max(min(rho, 0.99), -0.99)
  max(1, n * (1 - rho) / (1 + rho))
}

estimate_ar1_purge <- function(rows, max_purge = 25L) {
  if (!all(c("atom_index", "frame_slot", "target_T0") %in% names(rows))) return(0L)
  traj <- rows[is.finite(time_ps) | grepl("trajectory", pose_kind)]
  if (!nrow(traj)) return(0L)
  rhos <- traj[order(atom_index, frame_slot), {
    x <- target_T0
    if (.N < 8L) {
      list(rho = NA_real_)
    } else {
      list(rho = suppressWarnings(stats::cor(utils::head(x, -1), utils::tail(x, -1))))
    }
  }, by = atom_index]$rho
  rho <- stats::median(rhos[is.finite(rhos)], na.rm = TRUE)
  if (!is.finite(rho) || rho <= 0) return(1L)
  # AR1 is short-ranged here; use the e-folding horizon but cap so each block keeps training support.
  as.integer(max(1L, min(max_purge, ceiling(-1 / log(rho)))))
}

folds_leave_protein_out <- function(rows, k = NULL, seed = 1L) {
  proteins <- sort(unique(as.character(rows$protein_id)))
  if (is.null(k)) k <- length(proteins)
  k <- min(k, length(proteins))
  set.seed(seed)
  proteins <- sample(proteins)
  groups <- split(proteins, rep(seq_len(k), length.out = length(proteins)))
  lapply(seq_len(k), function(i) {
    test <- which(as.character(rows$protein_id) %in% groups[[i]])
    list(id = sprintf("protein_%02d", i), train = setdiff(seq_len(nrow(rows)), test), test = test, proteins = groups[[i]])
  })
}

folds_block_time <- function(rows, k = 3L, purge = NULL) {
  frames <- sort(unique(rows$frame_slot))
  if (!length(frames)) return(list())
  if (is.null(purge)) purge <- estimate_ar1_purge(rows)
  blocks <- split(frames, cut(seq_along(frames), breaks = k, labels = FALSE))
  lapply(seq_along(blocks), function(i) {
    test_frames <- blocks[[i]]
    lo <- min(test_frames) - purge
    hi <- max(test_frames) + purge
    test <- which(rows$frame_slot %in% test_frames)
    purged <- which(rows$frame_slot >= lo & rows$frame_slot <= hi)
    list(id = sprintf("time_%02d", i), train = setdiff(seq_len(nrow(rows)), purged), test = test, purge = purge)
  })
}

folds_leave_cluster <- function(rows, k = 3L, seed = 1L) {
  if (!"split_group_id" %in% names(rows)) {
    f <- folds_leave_protein_out(rows, k = k, seed = seed)
    attr(f, "sensitivity_flag") <- "split_group_id_absent_used_protein_id"
    return(f)
  }
  clusters <- sort(unique(as.character(rows$split_group_id)))
  set.seed(seed)
  clusters <- sample(clusters)
  groups <- split(clusters, rep(seq_len(min(k, length(clusters))), length.out = length(clusters)))
  lapply(seq_along(groups), function(i) {
    test <- which(as.character(rows$split_group_id) %in% groups[[i]])
    list(id = sprintf("cluster_%02d", i), train = setdiff(seq_len(nrow(rows)), test), test = test, clusters = groups[[i]])
  })
}

assert_disjoint_proteins <- function(rows, train, test) {
  tr <- unique(as.character(rows$protein_id[train]))
  te <- unique(as.character(rows$protein_id[test]))
  overlap <- intersect(tr, te)
  if (length(overlap)) abort("train/test protein overlap: %s", paste(head(overlap, 10), collapse = ", "))
  invisible(TRUE)
}

make_folds <- function(rows, k = 3L, seed = 1L, purge = NULL) {
  idx <- seq_len(nrow(rows))
  is_traj <- (("time_ps" %in% names(rows)) & is.finite(rows$time_ps)) |
    (("pose_kind" %in% names(rows)) & grepl("trajectory", rows$pose_kind))
  static_idx <- idx[!is_traj]
  traj_idx <- idx[is_traj]
  folds <- vector("list", k)
  static_parts <- if (length(static_idx)) folds_leave_protein_out(rows[static_idx, ], k = k, seed = seed) else vector("list", k)
  traj_parts <- if (length(traj_idx)) folds_block_time(rows[traj_idx, ], k = k, purge = purge) else vector("list", k)
  for (i in seq_len(k)) {
    static_test <- if (length(static_idx)) static_idx[static_parts[[i]]$test] else integer()
    static_train <- if (length(static_idx)) static_idx[static_parts[[i]]$train] else integer()
    traj_test <- if (length(traj_idx)) traj_idx[traj_parts[[i]]$test] else integer()
    traj_train <- if (length(traj_idx)) traj_idx[traj_parts[[i]]$train] else integer()
    folds[[i]] <- list(
      id = sprintf("fold_%02d", i),
      train = sort(c(static_train, traj_train)),
      test = sort(c(static_test, traj_test)),
      purge = if (length(traj_parts)) traj_parts[[i]]$purge else 0L
    )
  }
  covered <- sort(unique(unlist(lapply(folds, `[[`, "test"))))
  if (!identical(covered, idx)) abort("CV test folds do not cover every row exactly at least once")
  attr(folds, "design") <- "grouped leave-protein-out for static rows plus block-time trajectory folds"
  folds
}

fold_summary <- function(rows, folds) {
  data.table::rbindlist(lapply(folds, function(f) {
    data.table::data.table(
      fold = f$id,
      n_train = length(f$train),
      n_test = length(f$test),
      train_proteins = data.table::uniqueN(rows$protein_id[f$train]),
      test_proteins = data.table::uniqueN(rows$protein_id[f$test]),
      purge = f$purge %||% 0L
    )
  }))
}
