test_that("baseline uses training rows only", {
  rows <- data.table::data.table(
    .row_id = 1:6,
    element = c(1, 1, 1, 6, 6, 6),
    iupac_role = c(1, 1, 2, 1, 1, 2),
    residue_type = c(1, 1, 2, 1, 2, 2),
    equivalence_class = c(1, 1, 2, 1, 2, 2),
    target_T0 = c(1, 2, 1000, 4, 5, 1000)
  )
  fold <- list(train = c(1, 2, 4, 5), test = c(3, 6))
  fit <- baseline_blup(rows[fold$train, ])
  pred <- predict(fit, rows[fold$test, ])
  expect_true(all(pred < 100))
})
