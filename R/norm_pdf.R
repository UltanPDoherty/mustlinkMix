norm_pdf <- function(x, mu, sigma) {
  if (is.vector(x)) {
    x <- t(x)
  }

  x_n <- nrow(x)
  x_p <- ncol(x)

  if (x_n > 1) {
    x_centr <- x - matrix(1, x_n, 1) %*% mu
  } else {
    x_centr <- x - mu
  }

  norm_const <- 1 / (2 * pi)^(x_p / 2) * 1 / det(sigma)^(1 / 2)
  x_exp      <- exp(-1 / 2 * diag(x_centr %*% solve(sigma) %*% t(x_centr)))

  return(norm_const * x_exp)
}
