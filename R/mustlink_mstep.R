#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link
#' constraints.
#'
#' @inheritParams mustlink_em
#' @inheritParams mustlink_estep
#' @param postprob_events Expanded observation posterior probability matrix.
#' @param postprob_sets Posterior probability matrix for constrained sets.
#' @param sets_num Number of constrained sets including unconstrained events.
#'
#' @return List:
#' * prop
#' * mu
#' * sigma
#' @export
mustlink_mstep <- function(
    data,
    postprob_events,
    postprob_sets,
    sets_num = nrow(postprob_sets),
    clust_num = ncol(postprob_sets),
    event_num = nrow(data),
    var_num = ncol(data),
    drop_cluster = FALSE) {
  postprob_events_sums <- colSums(postprob_events)

  if (drop_cluster) {
    empty_clusters <- postprob_events_sums < 2
    if (any(empty_clusters)) {
      clust_num <- sum(!empty_clusters)
      postprob_events_sums <- postprob_events_sums[!empty_clusters]
      postprob_sets <- postprob_sets[, !empty_clusters, drop = FALSE]
      postprob_events <- postprob_events[, !empty_clusters, drop = FALSE]
    }
  }

  # block mixing proportions
  prop <- postprob_events_sums / event_num

  postprob_events_div <- sweep(
    x = postprob_events, MARGIN = 2,
    STATS = postprob_events_sums, FUN = "/"
  )

  # Mean vector
  mu <- t(postprob_events_div) %*% data

  # Covariance matrix
  data_mu <- array(dim = c(event_num, var_num, clust_num))
  sigma <- array(dim = c(var_num, var_num, clust_num))
  for (k in 1:clust_num) {
    data_mu[, , k] <-
      sqrt(postprob_events_div[, k]) * scale(data, mu[k, ], FALSE)
    sigma[, , k] <- crossprod(data_mu[, , k])
  }

  return(list(prop = prop, mu = mu, sigma = sigma))
}
