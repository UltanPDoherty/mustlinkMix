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
                                              "mlkmpp", "mlkm",
                                              "old_mlkmpp", "old_mlkm")) {

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
  plusplus = FALSE
) {
  linked_num <- length(unique(linked_set_labels)) - 1

  bigK <- 2 * clust_num
  if (plusplus) {
    bigKmeans <- ClusterR::KMeans_rcpp(
      data,
      clusters = bigK,
      seed = init_seed
    )
  } else {
    set.seed(init_seed)
    bigKmeans <- kmeans(
      data,
      centers = bigK,
      nstart = 3
    )
  }

  well <- rep(0, nrow(data))
  for (k in seq_len(bigK)) {
    bigKmeans_k <- bigKmeans$cluster == k
    for (p in seq_len(linked_num)) {
      linked_p <- linked_set_labels == p
      if (sum(linked_p & bigKmeans_k) > sum(bigKmeans_k) / 2) {
        well[bigKmeans_k] <- p
      }
    }
  }

  partition <- well
  partition[linked_set_labels != 0] <- linked_set_labels[linked_set_labels != 0]

  if (plusplus) {
    smallKmeans <- ClusterR::KMeans_rcpp(
      data[partition == 0, ],
      clusters = clust_num - linked_num,
      seed = init_seed
    )
  } else {
    set.seed(init_seed)
    smallKmeans <- kmeans(
      data[partition == 0, ],
      centers = clust_num - linked_num,
      nstart = 3
    )
  }

  partition[partition == 0] <- smallKmeans$cluster + linked_num

  return(partition)
}




old_mustlink_kmeans <- function(
  data,
  linked_set_labels,
  clust_num,
  init_seed,
  plusplus = FALSE
) {
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
