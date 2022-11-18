#' Compute the pdf of a multivariate Normal distribution.
#'
#' @param x A numeric vector giving the value of a single observation or a
#' matrix where each row is an observation.
#' @param mu A numeric vector giving the mean of the distribution.
#' @param sigma A matrix giving the covariance matrix of the distribution.
#'
#' @return A numeric vector with one entry per observation provided.
#' @export
#'
#' @examples
#' value <- norm_pdf(rep(0, 4), rep(0, 4), diag(4))
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
