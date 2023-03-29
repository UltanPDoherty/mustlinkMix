#' @title Compute initial parameter values of a GMM.
#'
#' @description
#' Initialise a Gaussian Mixture Model either with equal proportions, mean
#' vectors, and covariance matrices, or based on a k-Means clustering, or using
#' parameters computed by mclust.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param labels Set of initial labels.
#' @param start Initialisation option. One of "central", "k-Means", or "mclust".
#' @param init_seed Seed.
#'
#' @return A list consisting of a mixing proportions vector, a matrix of
#'         component means, and an array containing component covariance
#'         matrices.
#' @export
#' @import mclust
#'
#' @examples
#' initialise_model(iris[, 1:4], 4)
initialise_model <- function(data, clust_num, labels = NULL,
                             start = "k-Means", init_seed = NULL) {
  var_num <- ncol(data)
  obs_num <- nrow(data)

  if (!is.null(labels)) {
    clust_num <- length(unique(labels))
  }

  init <- list(prop  = vector("numeric", length = clust_num),
               mu    = matrix(NA, nrow = clust_num, ncol = var_num),
               sigma = array(NA, c(var_num, var_num, clust_num)))

  if (!is.null(labels)) {
    names     <- sort(unique(labels))
    sizes     <- as.numeric(table(labels))
    init$prop <- sizes / obs_num
    init$mu   <- rowsum(data, group = labels) / sizes
    for (k in 1:clust_num) {
      init$sigma[, , k] <- stats::cov(data[labels == names[k], ])
    }
    return(init)
  }

  start_options <- c("central", "k-Means", "mclust")
  stopifnot("start must be one of \"central\", \"k-Means\", or \"mclust\"." =
               any(start == start_options))

  if (!is.null(init_seed)) {
    set.seed(init_seed)
  }

  if (start == "central") {
    init$prop      <- rep(1 / clust_num, clust_num)
    init$mu        <- matrix(rep(colMeans(data), clust_num),
                             nrow = clust_num, byrow = FALSE)
    for (k in 1:clust_num) {
      init$sigma[, , k] <- diag(var_num)
    }
  } else if (start == "k-Means") {
    km <- stats::kmeans(x = data, centers = clust_num)
    init$prop  <- km$size / sum(km$size)
    init$mu    <- km$centers
    for (k in 1:clust_num) {
      init$sigma[, , k] <- stats::cov(data[km$cluster == k, ])
    }
  } else if (start == "mclust") {
    mc <- mclust::Mclust(data, G = clust_num, modelNames = c("VVV"))
    init$prop  <- mc$parameters$pro
    init$mu    <- t(mc$parameters$mean)
    init$sigma <- mc$parameters$variance$sigma
  } else {
    stop("start must be one of \"central\", \"k-Means\", or \"mclust\".")
  }

  return(init)
}
