#' @title Compute initial parameter values of a GMM from a vector of partition
#' labels.
#'
#' @description
#' Initialise a Gaussian Mixture Model.
#'
#' @param data Matrix or dataframe.
#' @param init_labels Set of initial labels.
#'
#' @return A list consisting of a mixing proportions vector, a matrix of
#'         component means, and an array containing component covariance
#'         matrices.
#' @export
initial_parameters <- function(data, init_labels) {
  stopifnot(
    length(init_labels) == nrow(data),
    min(table(init_labels)) > 1
  )

  params <- list()

  sizes <- as.numeric(table(init_labels))
  names <- sort(unique(init_labels))

  params$prop <- sizes / sum(sizes)

  params$mu <- rowsum(data, group = init_labels) / sizes

  clust_num <- length(unique(init_labels))
  var_num <- ncol(data)
  params$sigma <- array(NA, c(var_num, var_num, clust_num))
  for (k in seq_len(clust_num)) {
    params$sigma[, , k] <- stats::cov(data[init_labels == names[k], ])
  }

  return(params)
}
