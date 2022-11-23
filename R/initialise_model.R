#' Compute initial parameter values of a GMM.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param start Initialisation option.
#'
#' @return A list.
#' @export
#'
#' @examples
#' initialise_model(iris[, 1:4], 4)
initialise_model <- function(data, clust_num, start = "vanilla") {
  dimension <- ncol(data)

  if (start == "vanilla"){
    prop      <- rep(1 / clust_num, clust_num)
    mu        <- matrix(rep(colMeans(data), clust_num),
                        nrow = clust_num, byrow = FALSE)
    sigma     <- array(NA, c(clust_num, dimension, dimension))
    for (k in 1:clust_num) {
      sigma[k, , ] <- diag(dimension)
    }
  }

  if (start == "k-Means"){
    km    <- stats::kmeans(x = data, centers = clust_num)
    prop  <- km$size
    mu    <- km$centers
    sigma     <- array(NA, c(clust_num, dimension, dimension))
    for (k in 1:clust_num) {
      sigma[k, , ] <- stats::cov(data[km$cluster == k], )
    }
  }

  if (start == "mclust"){
    library(mclust)
    mc    <- mclust::Mclust(data, G = clust_num)
    prop  <- mc$parameters$pro
    mu    <- t(mc$parameters$mean)
    sigma     <- array(NA, c(clust_num, dimension, dimension))
    for (k in 1:clust_num) {
      sigma[k, , ] <- mc$parameters$variance$sigma[, , k]
    }
  }

  init <- list(prop = prop,
               mu = mu,
               sigma = sigma)
  return(init)
}
