#' @title Construct an initial partition.
#'
#' @description
#' Partition the data.
#'
#' @param data Matrix or dataframe.
#' @param clust_num Number of clusters.
#' @param linked_set_labels Label indicating which linked set an event
#' belongs to, 0 for unconstrained events.
#' @param init_method Initialisation option. One of "mlkmpp",
#' "mlkm", "kmpp", or "km".
#' @param init_seed Seed.
#'
#' @return An integer vector.
#' @export
initial_partition <- function(data, clust_num, linked_set_labels = NULL,
                              init_seed = NULL,
                              init_method = c("kmpp", "km",
                                              "mlkmpp",
                                              "mlkm")) {

  init_method <- rlang::arg_match(init_method)

  obs_num <- nrow(data)

  if (is.null(init_seed)) {
    init_seed <- as.numeric(Sys.time())
  }

  if (init_method == "km") {
    set.seed(init_seed)
    partition <- stats::kmeans(data, centers = clust_num)$cluster
  } else if (init_method == "kmpp") {
    partition <- ClusterR::KMeans_rcpp(data, clusters = clust_num,
                                       seed = init_seed)$clusters
  } else if (init_method %in% c("mlkm", "mlkmpp")) {
    linked_num <- length(unique(linked_set_labels[linked_set_labels != 0]))
    if (clust_num <= linked_num) {
      message(paste0("Cluster number is equal to number of constrained sets.",
                     "\n", "km is not implemented.",
                     "\n", "Constrained sets are used for initialisation.",
                     "\n"))
      partition <- linked_set_labels
    } else {
      set.seed(init_seed)
      if (init_method == "mlkm") {
        km <- stats::kmeans(data[linked_set_labels == 0, ],
                            centers = clust_num - linked_num)$cluster
      } else {
        km <- ClusterR::KMeans_rcpp(data[linked_set_labels == 0, ],
                                    clusters = clust_num - linked_num,
                                    seed = init_seed)$clusters
      }
      partition <- integer(obs_num)
      partition <- linked_set_labels
      partition[partition == 0] <- km + linked_num
    }
  }

  return(partition)
}
