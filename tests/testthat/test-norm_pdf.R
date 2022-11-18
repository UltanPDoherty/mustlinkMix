test_that("norm_pdf computes", {
  expect_equal(round(norm_pdf(rep(0, 2), rep(0, 2), diag(2)), digits = 7), 0.1591549)
})
