test_that("leave-protein folds are protein disjoint", {
  rows <- data.table::data.table(
    protein_id = rep(letters[1:6], each = 3),
    frame_slot = 1:18,
    target_T0 = rnorm(18)
  )
  folds <- folds_leave_protein_out(rows, k = 3, seed = 1)
  for (f in folds) expect_true(assert_disjoint_proteins(rows, f$train, f$test))
})

test_that("n_eff_ar1 is smaller for positive autocorrelation", {
  x <- as.numeric(stats::filter(rnorm(200), 0.8, method = "recursive"))
  expect_lt(n_eff_ar1(x), length(x))
  expect_gt(n_eff_ar1(x), 1)
})
