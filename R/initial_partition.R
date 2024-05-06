#' @title Construct an initial partition.
#'
#' @description
#' Partition the data.
#'
#' @inheritParams mustlink
#' @param linked_set_labels Label indicating which linked set an event
#' belongs to, 0 for unconstrained events.
#'
#' @return An integer vector.
#' @export
initial_partition <- function(data, clust_num, linked_set_labels = NULL,
                              init_seed = NULL,
                              init_method = c(
                                "mlkmpp", "mlkm",
                                "kmpp", "km",
                                "old_mlkmpp", "old_mlkm"
                              )) {
  init_method <- rlang::arg_match(init_method)

  if (is.null(init_seed)) {
    init_seed <- as.numeric(Sys.time())
  }

  if (init_method == "km") {
    set.seed(init_seed)
    partition <- stats::kmeans(data, centers = clust_num)$cluster
  } else if (init_method == "kmpp") {
    partition <- ClusterR::KMeans_rcpp(data,
      clusters = clust_num,
      seed = init_seed
    )$clusters
  } else if (init_method == "mlkm") {
    partition <- mustlink_kmeans(
      data,
      linked_set_labels,
      clust_num,
      init_seed,
      plusplus = FALSE
    )
  } else if (init_method == "mlkmpp") {
    partition <- mustlink_kmeans(
      data,
      linked_set_labels,
      clust_num,
      init_seed,
      plusplus = TRUE
    )
  } else if (init_method == "old_mlkm") {
    partition <- old_mustlink_kmeans(
      data,
      linked_set_labels,
      clust_num,
      init_seed,
      plusplus = FALSE
    )
  } else if (init_method == "old_mlkmpp") {
    partition <- old_mustlink_kmeans(
      data,
      linked_set_labels,
      clust_num,
      init_seed,
      plusplus = TRUE
    )
  }

  return(partition)
}


mustlink_kmeans <- function(
    data,
    linked_set_labels,
    clust_num,
    init_seed,
    plusplus = FALSE) {
  linked_num <- length(unique(linked_set_labels)) - 1

  big_k <- 2 * clust_num
  if (plusplus) {
    big_kmeans <- ClusterR::KMeans_rcpp(
      data,
      clusters = big_k,
      seed = init_seed
    )
  } else {
    set.seed(init_seed)
    big_kmeans <- kmeans(
      data,
      centers = big_k,
      nstart = 3
    )
  }

  well <- rep(0, nrow(data))
  for (k in seq_len(big_k)) {
    big_kmeans_k <- big_kmeans$cluster == k
    for (p in seq_len(linked_num)) {
      linked_p <- linked_set_labels == p
      if (sum(linked_p & big_kmeans_k) > sum(big_kmeans_k) / 2) {
        well[big_kmeans_k] <- p
      }
    }
  }

  partition <- well
  partition[linked_set_labels != 0] <- linked_set_labels[linked_set_labels != 0]

  if (plusplus) {
    small_kmeans <- ClusterR::KMeans_rcpp(
      data[partition == 0, ],
      clusters = clust_num - linked_num,
      seed = init_seed
    )
  } else {
    set.seed(init_seed)
    small_kmeans <- kmeans(
      data[partition == 0, ],
      centers = clust_num - linked_num,
      nstart = 3
    )
  }

  partition[partition == 0] <- small_kmeans$cluster + linked_num

  return(partition)
}




old_mustlink_kmeans <- function(
    data,
    linked_set_labels,
    clust_num,
    init_seed,
    plusplus = FALSE) {
  linked_num <- length(unique(linked_set_labels[linked_set_labels != 0]))

  if (clust_num <= linked_num) {
    message(paste0(
      "Cluster number is less than or equal to number of constrained sets.\n",
      "km is not implemented.\n",
      "Constrained sets are used for initialisation.\n"
    ))
    partition <- linked_set_labels
  } else {
    set.seed(init_seed)
    if (plusplus) {
      km <- ClusterR::KMeans_rcpp(
        data[linked_set_labels == 0, ],
        clusters = clust_num - linked_num,
        seed = init_seed
      )$clusters
    } else {
      km <- stats::kmeans(
        data[linked_set_labels == 0, ],
        centers = clust_num - linked_num
      )$cluster
    }

    partition <- linked_set_labels
    partition[partition == 0] <- km + linked_num
  }

  return(partition)
}
