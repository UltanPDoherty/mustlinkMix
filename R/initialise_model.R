#' Compute initial parameter values of a GMM.
#'
#' @param clust_num Number of clusters.
#' @param dimension Number of variables.
#'
#' @return A list.
#' @export
#'
#' @examples
#' initialise_model(3, 4)
initialise_model <- function(clust_num, dimension) {
  prop  <- rep(1 / clust_num, clust_num)
  mu    <- matrix(stats::rnorm(clust_num * dimension, 0, 1), ncol = dimension)
  sigma <- array(NA, c(clust_num, dimension, dimension))
  for (k in 1:clust_num) {
    sigma[k, , ] <- diag(dimension)
  }

  init <- list(prop = prop,
               mu = mu,
               sigma = sigma)
  return(init)
}
