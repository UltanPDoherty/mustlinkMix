#' Compute initial parameter values of a GMM.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param start Initialisation option.
#'
#' @return A list.
#' @export
#' @import mclust
#'
#' @examples
#' initialise_model(iris[, 1:4], 4)
initialise_model <- function(data, clust_num, start = "vanilla") {
  p <- ncol(data)

  sigma     <- array(NA, c(p, p, clust_num))

  if (start == "vanilla"){
    prop      <- rep(1 / clust_num, clust_num)
    mu        <- matrix(rep(colMeans(data), clust_num),
                        nrow = clust_num, byrow = FALSE)
    for (k in 1:clust_num) {
      sigma[, , k] <- diag(p)
    }
  }

  if (start == "k-Means"){
    km    <- stats::kmeans(x = data, centers = clust_num)
    prop  <- km$size / sum(km$size)
    mu    <- km$centers
    for (k in 1:clust_num) {
      sigma[, , k] <- stats::cov(data[km$cluster == k, ])
    }
  }

  if (start == "mclust"){
    #library(mclust)
    mc    <- mclust::Mclust(data, G = clust_num)
    prop  <- mc$parameters$pro
    mu    <- t(mc$parameters$mean)
    for (k in 1:clust_num) {
      sigma[, , k] <- mc$parameters$variance$sigma[, , k]
    }
  }

  init <- list(prop = prop,
               mu = mu,
               sigma = sigma)
  return(init)
}
