#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link
#' constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param postprob_event Expanded observation posterior probability matrix.
#' @param postprob_block Block posterior probability matrix.
#' @param block_num Number of blocks.
#' @param clust_num Number of clusters pre-specified.
#' @param event_num Number of observations in the dataset.
#' @param var_num Number of variables in the dataset.
#'
#' @return List containing prop, mu, sigma.
#' @export
#'
#' @examples
#' block_labels1 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' block1 <- list(labels = block_labels1,
#'                num = length(unique(block_labels1)),
#'                size = as.numeric(table(block_labels1)))
#' params1 <- initialise_model(iris[, 1:4], clust_num = 3)
#' e_out1  <- mustlink_estep_ns(as.matrix(iris[, 1:4]),
#'                              block = block1, params = params1,
#'                              event_num = 150, var_num = 4, clust_num = 3)
#' postprob_event1 <- e_out1$postprob_event
#' postprob_block1 <- e_out1$postprob_block
#' mustlink_mstep_ns(as.matrix(iris[, 1:4]),
#'                   postprob_event = postprob_event1, postprob_block = postprob_block1)
mustlink_mstep_ns <- function(data, postprob_event, postprob_block,
                              block_num = nrow(postprob_block),
                              clust_num = ncol(postprob_block),
                              event_num = nrow(data),
                              var_num = ncol(data)) {

  # Block mixing proportions
  prop <- colSums(postprob_block) / block_num

  postprob_event_sums <- colSums(postprob_event)
  postprob_event2 <- sweep(x = postprob_event, MARGIN = 2,
                           STATS = postprob_event_sums, FUN = "/")

  # Mean vector
  mu <- t(postprob_event2) %*% data

  # Covariance matrix
  data_mu <- array(dim = c(event_num, var_num, clust_num))
  sigma   <- array(dim = c(var_num, var_num, clust_num))
  for (k in 1:clust_num) {
    data_mu[, , k] <- sqrt(postprob_event2[, k]) * scale(data, center = mu[k, ],
                                                         scale = FALSE)
    sigma[, , k] <- crossprod(data_mu[, , k])
  }

  return(list(prop = prop,
              mu = mu,
              sigma = sigma
              )
         )
}
