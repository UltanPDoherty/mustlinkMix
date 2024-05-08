#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement E-step of EM algorithm for GMM with positive / must-link
#' constraints.
#'
#' @inheritParams mustlink_em
#' @param constraints_info .
#' @param event_num Number of observations in the data set.
#' @param var_num Number of variables in the data set.
#'
#' @return List:
#' * loglike: log likelihood value
#' * postprob_sets: posterior probability matrix
#' * postprob_events: posterior probability matrix
#' @export
mustlink_estep <- function(
    data,
    constraints_info,
    params,
    event_num = nrow(data),
    var_num = ncol(params$mu),
    clust_num = nrow(params$mu)) {
  lpdf_constrained <- compute_lpdf_constrained(
    data = data,
    constraints_info = constraints_info,
    params = params,
    event_num = event_num,
    clust_num = clust_num
  )

  # for loop computes the block posterior probability matrix, postprob_sets
  sets_unnorm <- postprob_sets <- matrix(
    nrow = constraints_info$num,
    ncol = clust_num
  )
  log_maxes <- loglike_vec <- sets_unnorm_sums <- vector(
    mode = "numeric",
    length = constraints_info$num
  )

  for (l in 1:constraints_info$num) {
    # Add the log mixing proportions and then un-log this sum with exp.
    # Subtract lpdf_constrained row maxes to prevent exp mapping large values to
    # Inf.
    log_maxes[l] <- max(
      lpdf_constrained[l, ] + constraints_info$size[l] * log(params$prop)
    )
    sets_unnorm[l, ] <- exp(
      lpdf_constrained[l, ]
      + constraints_info$size[l] * log(params$prop)
        - log_maxes[l]
    )

    # Normalise rows of block_unnorm to obtain postprob_sets.
    sets_unnorm_sums[l] <- sum(sets_unnorm[l, ])
    postprob_sets[l, ] <- sets_unnorm[l, ] / sets_unnorm_sums[l]

    loglike_vec[l] <- log(sets_unnorm_sums[l]) + log_maxes[l]
  }

  loglike <- sum(loglike_vec)

  postprob_events <- postprob_sets[constraints_info$labels, ]

  return(list(
    loglike = loglike,
    postprob_sets = postprob_sets,
    postprob_events = postprob_events
  ))
}

compute_lpdf_constrained <- function(
    data,
    constraints_info,
    params,
    event_num = nrow(data),
    clust_num = nrow(params$mu)) {
  # lpdf_event is an event_num x clust_num matrix
  # it is the log pdf for each component evaluated at every point
  lpdf_event <- vapply(
    1:clust_num,
    FUN.VALUE = double(event_num),
    FUN = function(k) {
      mvtnorm::dmvnorm(
        data,
        log = TRUE,
        mean = params$mu[k, ],
        sigma = params$sigma[, , k]
      )
    }
  )

  constrained_sets <- seq_len(constraints_info$zone_num)
  unconstrained_sets <- seq(constraints_info$zone_num + 1, constraints_info$num)

  constrained_events <- constraints_info$labels %in% constrained_sets
  unconstrained_events <- constraints_info$labels %in% unconstrained_sets

  # lpdf_constrained is the sum of lpdf_event values within each block
  lpdf_constrained <- matrix(NA, nrow = constraints_info$num, ncol = clust_num)
  lpdf_constrained[unconstrained_sets, ] <- lpdf_event[unconstrained_events, ]
  lpdf_constrained[constrained_sets, ] <- rowsum(
    lpdf_event[constrained_events, , drop = FALSE],
    group = constraints_info$labels[constrained_events]
  )

  return(lpdf_constrained)
}
