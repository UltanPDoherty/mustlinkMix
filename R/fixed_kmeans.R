#' Runs fixed_kmeans given constraints_common
#'
#' @param data .
#' @param constraints_common .
#' @param clust_num .
#' @param init_seed .
#'
#' @return label vector
#' @export
constrained_kmeans <- function(
    data, constraints_common, clust_num, init_seed = 123) {
  unconstrained_data <- data[constraints_common == 0, ]
  constrained_data <- data[constraints_common != 0, ]

  constrained_labels <- constraints_common[constraints_common != 0]

  constrained_sums <- rowsum(constrained_data, constrained_labels)
  constrained_counts <- as.numeric(table(constrained_labels))
  constraint_centres <- constrained_sums / constrained_counts

  unconstrained_labels <- fixed_kmeans(
    unconstrained_data, constraint_centres, clust_num, init_seed
  )

  labels <- rep(NA, nrow(data))
  labels[constraints_common == 0] <- unconstrained_labels
  labels[constraints_common != 0] <- constrained_labels

  return(labels)
}

#' k-Means with some fixed centres provided.
#'
#' @param data .
#' @param centres .
#' @param clust_num .
#' @param init_seed .
#' @param max_changes .
#'
#' @return label vector
#' @export
fixed_kmeans <- function(
    data, centres, clust_num, init_seed = 123, max_changes = 0) {
  obs_num <- nrow(data)
  fixed_num <- nrow(centres)

  extra_num <- clust_num - fixed_num

  if (extra_num > 0) {
    set.seed(init_seed)
    km0 <- stats::kmeans(data, centers = extra_num)
    added_centres <- km0$centers
    fixed_centres <- centres
    centres <- rbind(fixed_centres, added_centres)
  }

  label_changes <- Inf
  iter_count <- 0
  while (label_changes > max_changes) {
    dists <- matrix(nrow = obs_num, ncol = clust_num)
    for (k in seq_len(clust_num)) {
      dists[, k] <- apply(data, 1, function(x) sqrt(sum((x - centres[k, ])^2)))
    }

    new_labels <- apply(dists, 1, which.min)

    if (extra_num == 0) {
      label_changes <- 0
    } else if (iter_count == 0) {
      label_changes <- Inf
    } else {
      label_changes <- sum(label_vec != new_labels)
    }

    label_vec <- new_labels

    if (extra_num > 0) {
      for (k in seq_len(extra_num)) {
        added_centres[k, ] <- colMeans(data[label_vec == (k + fixed_num), ])
      }
      centres <- rbind(fixed_centres, added_centres)
    }

    iter_count <- iter_count + 1
    cat(paste0(
      "Iteration ", iter_count, ", Label Changes = ", label_changes, "\n"
    ))
  }

  return(label_vec)
}
