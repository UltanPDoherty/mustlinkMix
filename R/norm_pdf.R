#' Compute the pdf of a multivariate Normal distribution.
#'
#' @param x A numeric vector giving the value of a single observation or a
#' matrix where each row is an observation.
#' @param mu A numeric vector giving the mean of the distribution.
#' @param sigma_inv A matrix giving the inverse of the covariance matrix of the distribution.
#'
#' @return A numeric vector with one entry per observation provided.
#' @export
#'
#' @examples
#' norm_pdf(matrix(rep(2, 4), nrow = 2), mu = rep(1, 2), sigma_inv = diag(2))
norm_pdf <- function(x, mu, sigma_inv) {
  x_n <- nrow(x)
  x_centr <- x - matrix(1, x_n, 1) %*% mu

  out <- c()
  for(i in 1:x_n) {
    out[i]      <- -1 / 2 * diag(t(x_centr[i, ]) %*% sigma_inv %*% x_centr[i, ])
  }

  return(out)
}

