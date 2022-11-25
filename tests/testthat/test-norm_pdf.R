test_that("norm_pdf computes", {
  expect_equal(round(norm_pdf(matrix(rep(2, 4), nrow = 2), mu = rep(1, 2), sigma_inv = diag(2))),
               rep(-1, 2))
})
