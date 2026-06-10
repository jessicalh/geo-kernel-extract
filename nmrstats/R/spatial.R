morans_i <- function(resid, coords = NULL, adjacency = NULL) {
  x <- resid[is.finite(resid)]
  if (length(x) < 3L) return(NA_real_)
  z <- x - mean(x)
  denom <- sum(z^2)
  if (denom <= 0) return(NA_real_)
  if (!is.null(adjacency)) {
    pairs <- adjacency
    pairs <- pairs[pairs[, 1] <= length(z) & pairs[, 2] <= length(z), , drop = FALSE]
    if (!nrow(pairs)) return(NA_real_)
    wsum <- nrow(pairs)
    return(length(z) / wsum * sum(z[pairs[, 1]] * z[pairs[, 2]]) / denom)
  }
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    ok <- apply(coords, 1, function(r) all(is.finite(r)))
    coords <- coords[ok, , drop = FALSE]
    z <- resid[ok] - mean(resid[ok])
    if (nrow(coords) < 3L) return(NA_real_)
    d <- as.matrix(stats::dist(coords))
    w <- 1 / pmax(d, .Machine$double.eps)
    diag(w) <- 0
    wsum <- sum(w)
    if (wsum <= 0) return(NA_real_)
    return(length(z) / wsum * as.numeric(t(z) %*% w %*% z) / sum(z^2))
  }
  NA_real_
}

within_protein_morans_i <- function(predictions, rows) {
  pred <- data.table::as.data.table(predictions)
  if ("step" %in% names(pred)) {
    keep <- pred[["step"]] == "L6"
    pred <- pred[keep, ]
  }
  if (!nrow(pred)) return(data.table::data.table())
  meta <- data.table::copy(data.table::as.data.table(rows)[, .(protein_id_meta = protein_id, residue_index, element_meta = element)])
  meta[, .row_id_tmp := seq_len(.N)]
  pred <- merge(pred, meta, by.x = "row_id", by.y = ".row_id_tmp", all.x = TRUE)
  pred[, resid := y - pred]
  by_res <- pred[, .(resid = mean(resid)), by = .(protein_id = protein_id_meta, element = element_meta, residue_index)]
  by_res <- by_res[is.finite(residue_index), ]
  out <- by_res[order(residue_index), {
    n <- .N
    if (n < 3L) {
      list(n_residues = n, moran_i = NA_real_)
    } else {
      adjacency <- cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)
      adjacency <- rbind(adjacency, adjacency[, 2:1])
      list(n_residues = n, moran_i = morans_i(resid, adjacency = adjacency))
    }
  }, by = .(protein_id, element)]
  out[, note := "sequence-adjacency Moran proxy; coordinate columns absent from row_design emit"]
  out[]
}
