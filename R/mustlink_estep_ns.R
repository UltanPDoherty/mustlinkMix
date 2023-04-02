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
#'
#' @return A list containing a log-likelihood value
#'         and a posterior probability matrix.
#' @export
#'
#' @examples
#' block_labels1 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' block1 <- list(labels = block_labels1,
#'                num = length(unique(block_labels1)),
#'                size = as.numeric(table(block_labels1)))
#' params1 <- initialise_model(iris[, 1:4], clust_num = 3)
#' mustlink_estep_ns(as.matrix(iris[, 1:4]), block = block1, params = params1,
#'                   event_num = 150, var_num = 4, clust_num = 3)
mustlink_estep_ns <- function(data, block, params,
                           event_num = nrow(data), var_num = ncol(params$mu),
                           clust_num = nrow(params$mu)) {

  # lpdf_event is an event_num x clust_num matrix
  # it is the log pdf for each component evaluated at every point
  lpdf_event <- vapply(1:clust_num, FUN.VALUE = double(event_num),
                    FUN = function(k) {
                      mvtnorm::dmvnorm(data, log = TRUE,
                                       mean = params$mu[k, ],
                                       sigma = params$sigma[, , k])
                      }
                    )

  # singles_block is a logical of length block$num
  # it identifies which blocks are singletons,
  # i.e. unconstrained observations
  singles_block <- block$size == 1
  # singles_obs is a logical of length event_num
  # it identifies which observations correspond to the singleton blocks
  singles_obs   <- block$labels %in% (1:block$num)[singles_block]

  nonsingle_blocknum <- (1:block$num)[!singles_block]
  single_blocknum   <- block$labels[singles_obs]

  # lpdf_block is the sum of lpdf_event values within each block
  lpdf_block <- matrix(NA, nrow = block$num, ncol = clust_num)
  lpdf_block[single_blocknum, ] <- lpdf_event[singles_obs, ]
  lpdf_block[nonsingle_blocknum, ] <- rowsum(lpdf_event[!singles_obs, ],
                                             group = block$labels[!singles_obs])

  # for loop computes the block posterior probability matrix, postprob_block
  block_unnorm <- postprob_block <- matrix(nrow = block$num, ncol = clust_num)
  log_maxes <- loglike_vec <- block_unnorm_sums <- vector(mode = "numeric",
                                                     length = block$num)
  for (l in 1:block$num) {

    # Add the log mixing proportions and then un-log this sum with exp.
    # Subtract lpdf_block row maxes to prevent exp mapping large values to Inf.
    log_maxes[l]      <- max(lpdf_block[l, ])
    block_unnorm[l, ] <- exp(lpdf_block[l, ] - log_maxes[l] + log(params$prop))

    # Normalise rows of block_unnorm to obtain postprob_block.
    block_unnorm_sums[l] <- sum(block_unnorm[l, ])
    postprob_block[l, ] <- block_unnorm[l, ] / block_unnorm_sums[l]

    loglike_vec[l] <- log(block_unnorm_sums[l]) + log_maxes[l]
  }

  loglike  <- sum(loglike_vec)

  postprob_event <- postprob_block[block$labels, ]

  return(list(loglike = loglike,
              postprob_block = postprob_block,
              postprob_event = postprob_event))
}
