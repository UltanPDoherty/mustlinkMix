#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param obs_pp Expanded observation posterior probability matrix.
#' @param chunk_pp Chunklet posterior probability matrix.
#' @param chunk_num Number of chunklets.
#' @param clust_num Number of clusters pre-specified.
#' @param obs_num Number of observations in the dataset.
#' @param var_num Number of variables in the dataset.
#'
#' @return List containing prop, mu, sigma.
#' @export
#'
#' @examples
#' chunks_of_25 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' chunk1  <- make_chunk(iris[, 1:4], chunk_labs = chunks_of_25)
#' params1 <- initialise_model(iris[, 1:4], clust_num = 3)
#' e_out1  <- mustlink_estep(as.matrix(iris[, 1:4]),
#'                           chunk = chunk1, params = params1,
#'                           obs_num = 150, var_num = 4, clust_num = 3)
#' obs_pp1 <- e_out1$obs_pp
#' chunk_pp1 <- e_out1$chunk_pp
#' mustlink_mstep(as.matrix(iris[, 1:4]),
#'                obs_pp = obs_pp1, chunk_pp = chunk_pp1)
mustlink_mstep <- function(data, obs_pp, chunk_pp,
                           chunk_num = nrow(chunk_pp),
                           clust_num = ncol(chunk_pp),
                           obs_num = nrow(data),
                           var_num = ncol(data)) {

  # Chunklet mixing proportions
  prop <- colSums(chunk_pp) / chunk_num

  obs_pp_sums <- colSums(obs_pp)
  obs_pp2 <- sweep(x = obs_pp, MARGIN = 2, STATS = obs_pp_sums, FUN = "/")

  # Mean vector
  mu <- t(obs_pp2) %*% data

  # Covariance matrix
  data_mu <- array(dim = c(obs_num, var_num, clust_num))
  sigma   <- array(dim = c(var_num, var_num, clust_num))
  for(k in 1:clust_num) {
    data_mu[, , k] <- sqrt(obs_pp2[, k]) * scale(data, center = mu[k, ], scale = FALSE)
    sigma[, , k] <- crossprod(data_mu[, , k])
  }

  return(list(prop = prop,
              mu = mu,
              sigma = sigma
              )
         )
}
