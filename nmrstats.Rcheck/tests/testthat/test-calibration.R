test_that("normal calibration is close on a well specified synthetic case", {
  set.seed(1)
  pred <- data.table::data.table(
    row_id = 1:5000,
    step = "L6",
    element = 1,
    y = rnorm(5000),
    pred = 0
  )
  cal <- calibration_coverage(pred, levels = 0.90)
  expect_equal(cal$coverage$nominal, 0.90)
  expect_lt(abs(cal$coverage$empirical - 0.90), 0.03)
})
