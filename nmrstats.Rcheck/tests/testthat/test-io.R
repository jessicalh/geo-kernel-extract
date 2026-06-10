test_that("fixture-tagged emit is refused", {
  dir <- tempfile("emit")
  dir.create(dir)
  data.table::fwrite(data.table::data.table(
    dataset_id = "d", protein_id = "p", atom_index = 1L, frame_slot = 1L,
    sin_phi = 0, cos_phi = 1, sin_psi = 0, cos_psi = 1, target_T0 = 1
  ), file.path(dir, "row_design_rows.csv"))
  jsonlite::write_json(list(fixture = TRUE, row_count = 1, subsampling = "none"), file.path(dir, "schema_audit.json"), auto_unbox = TRUE)
  jsonlite::write_json(list(), file.path(dir, "column_provenance.json"), auto_unbox = TRUE)
  jsonlite::write_json(list(), file.path(dir, "null_spec.json"), auto_unbox = TRUE)
  jsonlite::write_json(list(), file.path(dir, "region_spec.json"), auto_unbox = TRUE)
  jsonlite::write_json(list(columns = list()), file.path(dir, "column_irrep_schema.json"), auto_unbox = TRUE)
  data.table::fwrite(data.table::data.table(metric = "rows", count = 1L), file.path(dir, "support_table.csv"))
  RcppCNPy::npySave(file.path(dir, "row_design_target_T2.npy"), matrix(0, nrow = 1, ncol = 5))
  expect_error(load_substrate(dir), "fixture")
})

test_that("bare EFG kernels are not default ppm features", {
  rows <- data.table::data.table(
    ring_bs_T0 = c(1, 2),
    mopac_bare_efg_kernel_absT2 = c(3, 4),
    apbs_bare_efg_kernel_absT2 = c(5, 6)
  )
  blocks <- feature_blocks(rows)
  expect_false(any(grepl("bare_efg_kernel", unlist(blocks, use.names = FALSE))))
  expect_true(length(attr(blocks, "excluded_bare_efg_kernel")) > 0)
})
