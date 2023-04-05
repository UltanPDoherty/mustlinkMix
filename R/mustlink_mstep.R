#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link
#' constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param postprob_event Expanded observation posterior probability matrix.
#' @param postprob_block block posterior probability matrix.
#' @param block_num Number of blocks.
#' @param clust_num Number of clusters pre-specified.
#' @param event_num Number of observations in the dataset.
#' @param var_num Number of variables in the dataset.
#' @param model Model to be used. Either "vm" for Melnykov et al. or "ns" for
#' Shental et al.
#'
#' @return List containing prop, mu, sigma.
#' @export
mustlink_mstep <- function(data, postprob_event, postprob_block,
                           block_num = nrow(postprob_block),
                           clust_num = ncol(postprob_block),
                           event_num = nrow(data),
                           var_num = ncol(data),
                           model = c("vm", "ns")) {

  model <- rlang::arg_match(model)

  # block mixing proportions
  prop <- switch(model,
                 vm = colSums(postprob_event) / event_num,
                 ns = colSums(postprob_block) / block_num)

  postprob_event_sums <- colSums(postprob_event)
  postprob_event2 <- sweep(x = postprob_event, MARGIN = 2,
                           STATS = postprob_event_sums, FUN = "/")

  # Mean vector
  mu <- t(postprob_event2) %*% data

  # Covariance matrix
  data_mu <- array(dim = c(event_num, var_num, clust_num))
  sigma   <- array(dim = c(var_num, var_num, clust_num))
  for (k in 1:clust_num){
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
