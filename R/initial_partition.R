#' @title Construct an initial partition.
#'
#' @description
#' Partition the data.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param constraint_labels Label indicating which constrained set an event
#' belongs to, 0 for unconstrained events.
#' @param init_method Initialisation option. One of "Must-Link k-Means++",
#' "k-Means++", "k-Means", "Mclust", or "use_labels".
#' @param init_seed Seed.
#'
#' @return An integer vector.
#' @export
#' @import mclust
#'
#' @examples
#' initialise_model(iris[, 1:4], 4)
initial_partition <- function(data, clust_num, constraint_labels = NULL,
                             init_seed = NULL,
                             init_method = c("Must-Link k-Means++",
                                             "k-Means++", "k-Means")) {

  init_method <- rlang::arg_match(init_method)

  obs_num <- nrow(data)

  set.seed(init_seed)

  if (init_method == "k-Means") {
    partition <- stats::kmeans(data, centers = clust_num)$cluster
  } else if (init_method == "k-Means++") {
    partition <- ClusterR::KMeans_rcpp(data, clusters = clust_num)$clusters
  } else if (init_method == "Must-Link k-Means++") {
    constr_num <- length(unique(constraint_labels[constraint_labels != 0]))
    kmpp <- ClusterR::KMeans_rcpp(data[constraint_labels == 0, ],
                                  clusters = clust_num - constr_num)$clusters

    partition <- integer(obs_num)
    partition <- constraint_labels
    partition[partition == 0] <- kmpp + constr_num
  }

  return(partition)
}
