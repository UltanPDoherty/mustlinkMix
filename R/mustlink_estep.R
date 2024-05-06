#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement E-step of EM algorithm for GMM
#' with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param params List containing prop, mu, sigma.
#' @param block Object from make_block.
#' @param event_num Number of observations in the dataset.
#' @param var_num Number of varibles in the dataset.
#' @param clust_num Number of clusters pre-specified.
#' @param model Model to be used. Either "vm" for Melnykov et al. or "ns" for
#' Shental et al.
#'
#' @return A list containing a log-likelihood value
#'         and a posterior probability matrix.
#' @export
mustlink_estep <- function(data, block, params,
                           event_num = nrow(data), var_num = ncol(params$mu),
                           clust_num = nrow(params$mu),
                           model = c("vm", "ns")) {
  model <- rlang::arg_match(model)

  lpdf_block <- compute_lpdf_block(
    data = data, block = block, params = params,
    event_num = event_num, clust_num = clust_num
  )

  # for loop computes the block posterior probability matrix, postprob_block
  block_unnorm <- postprob_block <- matrix(nrow = block$num, ncol = clust_num)
  log_maxes <- loglike_vec <- block_unnorm_sums <- vector(
    mode = "numeric",
    length = block$num
  )

  prop_exponent <- switch(model,
    vm = block$size,
    ns = rep(1, block$num)
  )
  for (l in 1:block$num) {
    # Add the log mixing proportions and then un-log this sum with exp.
    # Subtract lpdf_block row maxes to prevent exp mapping large values to Inf.
    log_maxes[l] <- max(lpdf_block[l, ] + prop_exponent[l] * log(params$prop))
    block_unnorm[l, ] <- exp(
      lpdf_block[l, ] + prop_exponent[l] * log(params$prop) - log_maxes[l]
    )

    # Normalise rows of block_unnorm to obtain postprob_block.
    block_unnorm_sums[l] <- sum(block_unnorm[l, ])
    postprob_block[l, ] <- block_unnorm[l, ] / block_unnorm_sums[l]

    loglike_vec[l] <- log(block_unnorm_sums[l]) + log_maxes[l]
  }

  loglike <- sum(loglike_vec)

  postprob_event <- postprob_block[block$labels, ]

  return(list(
    loglike = loglike,
    postprob_block = postprob_block,
    postprob_event = postprob_event
  ))
}




compute_lpdf_block <- function(data, block, params,
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

  linked_blocks <- seq_len(block$zone_num)
  unlinked_blocks <- seq(block$zone_num + 1, block$num)

  linked_events <- block$labels %in% linked_blocks
  unlinked_events <- block$labels %in% unlinked_blocks

  # lpdf_block is the sum of lpdf_event values within each block
  lpdf_block <- matrix(NA, nrow = block$num, ncol = clust_num)
  lpdf_block[unlinked_blocks, ] <- lpdf_event[unlinked_events, ]
  lpdf_block[linked_blocks, ] <- rowsum(
    lpdf_event[linked_events, , drop = FALSE],
    group = block$labels[linked_events]
  )

  return(lpdf_block)
}
