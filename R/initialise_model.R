#' @title Compute initial parameter values of a GMM.
#'
#' @description
#' Initialise a Gaussian Mixture Model either with equal proportions, mean
#' vectors, and covariance matrices, or based on a k-Means clustering, or using
#' parameters computed by mclust.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param constraint_labels Label indicating which constrained set an event
#' belongs to, 0 for unconstrained events.
#' @param init_labels Set of initial labels.
#' @param init_method Initialisation option. One of "Must-Link k-Means++",
#' "k-Means++", "k-Means", "Mclust", or "use_labels".
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
initialise_model <- function(data, clust_num, constraint_labels = NULL,
                             init_labels = NULL, init_seed = NULL,
                             init_method = c("Must-Link k-Means++",
                                             "k-Means++", "k-Means",
                                             "Mclust", "use_labels")) {

  init_method <- rlang::arg_match(init_method)

  var_num <- ncol(data)
  obs_num <- nrow(data)
  params <- list()

  set.seed(init_seed)

  if (init_method == "Must-Link k-Means++") {
    constr_num <- length(unique(constraint_labels[constraint_labels != 0]))
    kmpp <- ClusterR::KMeans_rcpp(data[constraint_labels == 0, ],
                                  clusters = clust_num - constr_num)

    combined_labels <- integer(obs_num)
    combined_labels <- constraint_labels
    combined_labels[constraint_labels == 0] <- kmpp$clusters + constr_num

    sizes        <- as.numeric(table(combined_labels))
    params$prop  <- sizes / sum(sizes)
    params$mu    <- rowsum(data, group = combined_labels) / sizes
    params$sigma <- array(NA, c(var_num, var_num, clust_num))
    for (k in 1:clust_num) {
      params$sigma[, , k] <- stats::cov(data[combined_labels == k, ])
    }
  }

  if (init_method == "k-Means++") {
    kmpp <- ClusterR::KMeans_rcpp(data, clusters = clust_num)
    sizes        <- as.numeric(table(kmpp$clusters))
    params$prop  <- sizes / sum(sizes)
    params$mu    <- kmpp$centroids
    params$sigma <- array(NA, c(var_num, var_num, clust_num))
    for (k in 1:clust_num) {
      params$sigma[, , k] <- stats::cov(data[kmpp$clusters == k, ])
    }
  }

  if (init_method == "k-Means") {
    km <- stats::kmeans(data, centers = clust_num)
    params$prop  <- km$size / sum(km$size)
    params$mu    <- km$centers
    params$sigma <- array(NA, c(var_num, var_num, clust_num))
    for (k in 1:clust_num) {
      params$sigma[, , k] <- stats::cov(data[km$cluster == k, ])
    }
  }

  if (init_method == "Mclust") {
    mc <- mclust::Mclust(data, G = clust_num, modelNames = c("VVV"))
    params$prop  <- mc$parameters$pro
    params$mu    <- t(mc$parameters$mean)
    params$sigma <- mc$parameters$variance$sigma
  }

  if (init_method == "use_labels") {
    sizes      <- as.numeric(table(init_labels))

    stopifnot(!is.null(init_labels),
              length(init_labels) == nrow(data),
              min(table(init_labels)) > 1)

    params$prop  <- sizes / obs_num
    params$mu    <- rowsum(data, group = init_labels) / sizes

    clust_num  <- length(unique(init_labels))
    params$sigma <- array(NA, c(var_num, var_num, clust_num))
    names      <- sort(unique(init_labels))
    for (k in 1:clust_num) {
      params$sigma[, , k] <- stats::cov(data[init_labels == names[k], ])
    }
  }

  return(params)
}
