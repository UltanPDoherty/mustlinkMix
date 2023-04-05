#' @title Construct an initial partition.
#'
#' @description
#' Partition the data.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param linked_set_labels Label indicating which linked set an event
#' belongs to, 0 for unconstrained events.
#' @param init_method Initialisation option. One of "Must-Link k-Means++",
#' "Must-Link k-Means", "k-Means++", or "k-Means".
#' @param init_seed Seed.
#'
#' @return An integer vector.
#' @export
#' @import mclust
#'
#' @examples
#' initialise_model(iris[, 1:4], 4)
initial_partition <- function(data, clust_num, linked_set_labels = NULL,
                              init_seed = NULL,
                              init_method = c("k-Means++", "k-Means",
                                              "Must-Link k-Means++",
                                              "Must-Link k-Means")) {

  init_method <- rlang::arg_match(init_method)

  obs_num <- nrow(data)

  if (is.null(init_seed)) {
    init_seed <- as.numeric(Sys.time())
  }

  if (init_method == "k-Means") {
    set.seed(init_seed)
    partition <- stats::kmeans(data, centers = clust_num)$cluster
  } else if (init_method == "k-Means++") {
    partition <- ClusterR::KMeans_rcpp(data, clusters = clust_num,
                                       seed = init_seed)$clusters
  } else if (init_method == "Must-Link k-Means") {
    linked_num <- length(unique(linked_set_labels[linked_set_labels != 0]))
    set.seed(init_seed)
    km <- stats::kmeans(data[linked_set_labels == 0, ],
                        centers = clust_num - linked_num)$cluster

    partition <- integer(obs_num)
    partition <- linked_set_labels
    partition[partition == 0] <- km + linked_num
  } else if (init_method == "Must-Link k-Means++") {
    linked_num <- length(unique(linked_set_labels[linked_set_labels != 0]))
    kmpp <- ClusterR::KMeans_rcpp(data[linked_set_labels == 0, ],
                                  clusters = clust_num - linked_num,
                                  seed = init_seed)$clusters

    partition <- integer(obs_num)
    partition <- linked_set_labels
    partition[partition == 0] <- kmpp + linked_num
  }

  return(partition)
}
