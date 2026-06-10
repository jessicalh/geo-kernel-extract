`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

read_json_file <- function(path) {
  jsonlite::read_json(path, simplifyVector = FALSE)
}

load_substrate <- function(dir, load_t2 = TRUE, load_embedding = FALSE) {
  if (length(dir) > 1) {
    parts <- lapply(dir, load_substrate, load_t2 = load_t2, load_embedding = load_embedding)
    headers <- lapply(parts, names)
    if (!all(vapply(headers[-1], identical, logical(1), headers[[1]]))) {
      abort("row_design headers are not byte-column-compatible across inputs")
    }
    rows <- data.table::rbindlist(parts, use.names = TRUE, fill = FALSE)
    attrs <- list(
      manifests = lapply(parts, attr, "manifest"),
      support_global = data.table::rbindlist(lapply(parts, attr, "support_global"), fill = TRUE),
      null_spec = lapply(parts, attr, "null_spec"),
      provenance = lapply(parts, attr, "provenance"),
      irrep_schema = lapply(parts, attr, "irrep_schema"),
      region_placeholder = lapply(parts, attr, "region_placeholder"),
      target_T2 = if (load_t2) do.call(rbind, lapply(parts, attr, "target_T2")) else NULL
    )
    for (nm in names(attrs)) attr(rows, nm) <- attrs[[nm]]
    data.table::setkeyv(rows, intersect(c("dataset_id", "protein_id", "atom_index", "frame_slot"), names(rows)))
    rows[]
  } else {
    dir <- normalizePath(dir, mustWork = TRUE)
    required <- c(
      "row_design_rows.csv", "schema_audit.json", "column_provenance.json",
      "support_table.csv", "null_spec.json", "region_spec.json",
      "column_irrep_schema.json", "row_design_target_T2.npy"
    )
    missing <- required[!file.exists(file.path(dir, required))]
    if (length(missing)) abort("missing emit files in %s: %s", dir, paste(missing, collapse = ", "))

    schema <- read_json_file(file.path(dir, "schema_audit.json"))
    if (isTRUE(schema$fixture)) abort("refusing fixture-tagged table: %s", dir)
    if (!identical(schema$subsampling %||% "none", "none")) {
      abort("refusing subsampled table %s: subsampling=%s", dir, schema$subsampling)
    }

    rows <- data.table::fread(file.path(dir, "row_design_rows.csv"), showProgress = FALSE)
    expected <- schema$expected_rows %||% schema$row_count
    if (!is.null(expected) && nrow(rows) != as.integer(expected)) {
      abort("row count mismatch for %s: got %d expected %d", dir, nrow(rows), as.integer(expected))
    }
    if (!"target_T0" %in% names(rows)) abort("target_T0 missing from %s", dir)
    if (any(!is.finite(rows$target_T0))) abort("target_T0 contains non-finite values in %s", dir)

    if (!all(c("sin_phi", "cos_phi", "sin_psi", "cos_psi") %in% names(rows))) {
      abort("phi/psi sin-cos columns missing from %s", dir)
    }
    rows[, phi_from_sc := atan2(sin_phi, cos_phi)]
    rows[, psi_from_sc := atan2(sin_psi, cos_psi)]
    if (!"phi" %in% names(rows)) rows[, phi := phi_from_sc]
    if (!"psi" %in% names(rows)) rows[, psi := psi_from_sc]
    for (j in 1:4) {
      s <- sprintf("sin_chi%d", j)
      c <- sprintf("cos_chi%d", j)
      out <- sprintf("chi%d", j)
      if (all(c(s, c) %in% names(rows))) rows[, (out) := atan2(get(s), get(c))]
    }
    rows[, emit_dir := basename(dir)]

    support <- data.table::fread(file.path(dir, "support_table.csv"), showProgress = FALSE)
    target_t2 <- NULL
    if (load_t2) {
      target_t2 <- RcppCNPy::npyLoad(file.path(dir, "row_design_target_T2.npy"))
      if (nrow(target_t2) != nrow(rows)) {
        abort("target_T2 row count mismatch in %s: got %d expected %d", dir, nrow(target_t2), nrow(rows))
      }
    }
    if (load_embedding) {
      emb <- RcppCNPy::npyLoad(file.path(dir, "row_design_aimnet2_embedding.npy"))
      if (nrow(emb) != nrow(rows)) abort("embedding row count mismatch in %s", dir)
      attr(rows, "aimnet2_embedding") <- emb
    }

    attr(rows, "manifest") <- schema
    attr(rows, "support_global") <- support
    attr(rows, "null_spec") <- read_json_file(file.path(dir, "null_spec.json"))
    attr(rows, "provenance") <- read_json_file(file.path(dir, "column_provenance.json"))
    attr(rows, "irrep_schema") <- read_json_file(file.path(dir, "column_irrep_schema.json"))
    attr(rows, "region_placeholder") <- read_json_file(file.path(dir, "region_spec.json"))
    attr(rows, "target_T2") <- target_t2
    data.table::setkeyv(rows, intersect(c("dataset_id", "protein_id", "atom_index", "frame_slot"), names(rows)))
    rows[]
  }
}

schema_columns <- function(rows) {
  manifests <- attr(rows, "manifests") %||% list(attr(rows, "manifest"))
  cols <- manifests[[1]]$columns
  if (is.null(cols)) return(data.table::data.table(name = names(rows)))
  data.table::rbindlist(lapply(cols, data.table::as.data.table), fill = TRUE)
}

irrep_columns <- function(rows) {
  schemas <- attr(rows, "irrep_schema")
  if (is.null(schemas)) schemas <- list(attr(rows, "irrep_schema"))
  data.table::rbindlist(lapply(schemas, function(x) {
    if (is.null(x$columns)) data.table::data.table() else data.table::rbindlist(lapply(x$columns, data.table::as.data.table), fill = TRUE)
  }), fill = TRUE)
}

support_per_stratum <- function(rows) {
  req <- c("dataset_id", "element", "residue_type", "iupac_role", "equivalence_class", "rama_region", "dssp_ss3")
  by <- intersect(req, names(rows))
  if (!length(by)) abort("no support strata columns are present")
  rows[, atom_uid_tmp := paste(protein_id, atom_index, sep = ":")]
  rows[, residue_uid_tmp := paste(protein_id, residue_index, sep = ":")]
  out <- rows[, .(
    n = .N,
    n_atoms = data.table::uniqueN(atom_uid_tmp),
    n_residues = data.table::uniqueN(residue_uid_tmp),
    n_proteins = data.table::uniqueN(protein_id),
    n_frames = data.table::uniqueN(frame_slot),
    target_sd = stats::sd(target_T0)
  ), by = by]
  rows[, c("atom_uid_tmp", "residue_uid_tmp") := NULL]
  data.table::setorder(out, dataset_id, element, -n)
  out[]
}

feature_blocks <- function(rows, allow_bare_efg = FALSE) {
  present <- function(cols) cols[cols %in% names(rows) & vapply(cols, function(z) any(is.finite(rows[[z]])), logical(1))]
  ring <- present(c("ring_bs_T0", "ring_bs_absT2", "ring_hm_T0", "ring_hm_absT2", "ring_jb_T0", "ring_jb_absT2"))
  mc <- present(c("mc_lit_T0", "mc_lit_absT2"))
  field_charge <- present(c(
    "ff14sb_efield_mag", "mopac_efield_mag", "apbs_efield_mag",
    "aimnet2_charge", "ff14sb_charge", "mopac_charge", "eeq_charge", "eeq_cn",
    "sasa", "larsen_hbond_count", "larsen_hbond_absT2"
  ))
  bare_efg <- present(c("mopac_bare_efg_kernel_absT2", "apbs_bare_efg_kernel_absT2", "aimnet2_bare_efg_kernel_absT2"))
  blocks <- list(ring = ring, mc = mc, field_charge = field_charge)
  if (allow_bare_efg) blocks$bare_efg_kernel <- bare_efg
  blocks <- blocks[vapply(blocks, length, integer(1)) > 0]
  attr(blocks, "excluded_bare_efg_kernel") <- bare_efg
  blocks
}

spatial_features <- function(rows) {
  cols <- c(
    grep("^n_(atoms|rings|charges|bonds)_[0-9]+A$", names(rows), value = TRUE),
    "nearest_ring_dist", "nearest_charge_dist", "nearest_charge_sign",
    "nearest_bond_dist", "nearest_atom_dist", "ring_cyl_z", "ring_cyl_rho",
    "ring_angle_to_normal", "self_or_bonded_atom_count", "self_or_bonded_bond_count"
  )
  cols[cols %in% names(rows)]
}

factor_model_cols <- function(rows) {
  intersect(c(
    "dataset_id", "protein_id", "residue_type", "iupac_role",
    "backbone_role", "equivalence_class", "prev_class", "next_class",
    "pre_proline", "dssp_ss3", "region_def_id", "rama_region",
    "residue_class", "terminal_state"
  ), names(rows))
}

prepare_model_data <- function(rows, extra_numeric = character()) {
  dt <- data.table::copy(rows)
  for (col in factor_model_cols(dt)) {
    x <- dt[[col]]
    if (col %in% c("dataset_id", "protein_id", "region_def_id")) {
      x <- ifelse(is.na(x) | x == "", "missing", as.character(x))
    } else {
      x <- ifelse(is.na(x) | x == -1 | x == "", "missing", paste0("v", x))
    }
    data.table::set(dt, j = col, value = factor(x))
  }
  dt[, atom_uid := factor(paste(protein_id, atom_index, sep = ":"))]
  dt[, phi_present := as.integer(is.finite(phi) & is.finite(psi))]
  dt[, phi_model := ifelse(is.finite(phi), phi, 0)]
  dt[, psi_model := ifelse(is.finite(psi), psi, 0)]
  for (j in 1:4) {
    chi <- sprintf("chi%d", j)
    exists <- sprintf("chi%d_exists", j)
    model <- sprintf("chi%d_model", j)
    if (chi %in% names(dt)) {
      if (!exists %in% names(dt)) dt[, (exists) := as.integer(is.finite(get(chi)))]
      dt[, (model) := ifelse(is.finite(get(chi)) & get(exists) == 1, get(chi), 0)]
      dt[, (exists) := as.integer(!is.na(get(exists)) & get(exists) == 1)]
    }
  }
  numeric_cols <- unique(c(
    spatial_features(dt),
    unlist(feature_blocks(dt), use.names = FALSE),
    extra_numeric
  ))
  for (col in intersect(numeric_cols, names(dt))) {
    x <- dt[[col]]
    if (!is.numeric(x) && !is.integer(x)) next
    if (all(!is.finite(x))) {
      data.table::set(dt, j = col, value = rep(0, nrow(dt)))
    } else {
      med <- stats::median(x[is.finite(x)], na.rm = TRUE)
      x[!is.finite(x)] <- med
      data.table::set(dt, j = col, value = as.numeric(x))
    }
  }
  dt[]
}

write_dt <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path)
  invisible(path)
}
