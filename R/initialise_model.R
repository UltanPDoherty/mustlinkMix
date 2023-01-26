#' @title Compute initial parameter values of a GMM.
#'
#' @description
#' Initialise a Gaussian Mixture Model either with equal proportions, mean
#' vectors, and covariance matrices, or based on a k-Means clustering, or using
#' parameters computed by mclust.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param start Initialisation option. One of "central", "k-Means", or "mclust".
#' @param init_seed Seed.
#'
#' @return A list consisting of a mixing proportions vector, a matrix of component means,
#'         and an array containing component covariance matrices.
#' @export
#' @import mclust
#'
#' @examples
#' initialise_model(iris[, 1:4], 4)
initialise_model <- function(data, clust_num, start = "central", init_seed = NULL) {
  p <- ncol(data)

  sigma     <- array(NA, c(p, p, clust_num))

  start_options <- c("central", "k-Means", "mclust")
  stopifnot( "start must be one of \"central\", \"k-Means\", or \"mclust\"." =
               any(start == start_options))

  if (!is.null(init_seed)) {
    set.seed(init_seed)
  }

  if (start == "central") {
    prop      <- rep(1 / clust_num, clust_num)
    mu        <- matrix(rep(colMeans(data), clust_num),
                        nrow = clust_num, byrow = FALSE)
    for (k in 1:clust_num) {
      sigma[, , k] <- diag(p)
    }
  } else if (start == "k-Means") {
    km    <- stats::kmeans(x = data, centers = clust_num)
    prop  <- km$size / sum(km$size)
    mu    <- km$centers
    for (k in 1:clust_num) {
      sigma[, , k] <- stats::cov(data[km$cluster == k, ])
    }
  } else if (start == "mclust") {
    mc    <- mclust::Mclust(data, G = clust_num, modelNames = c("VVV"))
    prop  <- mc$parameters$pro
    mu    <- t(mc$parameters$mean)
    sigma[, , k] <- mc$parameters$variance$sigma[, , k]
  } else {
    stop("start must be one of \"central\", \"k-Means\", or \"mclust\".")
  }

  init <- list(prop = prop,
               mu = mu,
               sigma = sigma)
  return(init)
}
