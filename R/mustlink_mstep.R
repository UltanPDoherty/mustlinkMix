#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param chunk Object from make_chunk.
#' @param obs_pp Expanded observation posterior probability matrix.
#' @param chunk_pp Chunklet posterior probability matrix.
#' @param obs_num Number of observations in the dataset.
#' @param var_num Number of varibles in the dataset.
#' @param clust_num Number of clusters pre-specified.
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
#' mustlink_mstep(as.matrix(iris[, 1:4]), chunk = chunk1,
#'                obs_pp = obs_pp1, chunk_pp = chunk_pp1,
#'                obs_num = 150, var_num = 4, clust_num = 3)
mustlink_mstep <- function(data, chunk, obs_pp, chunk_pp,
                           obs_num = nrow(data),
                           var_num = ncol(data),
                           clust_num = ncol(chunk_pp)) {

  # Chunklet mixing proportions
  prop <- colSums(chunk_pp) / chunk$num

  obs_pp_sums <- colSums(obs_pp)

  # Mean vector
  mu      <- t(obs_pp) %*% data / obs_pp_sums

  # Covariance matrix
  sigma0 <- data_mu <- list()
  sigma  <- array(NA, dim = c(var_num, var_num, clust_num))
  for(k in 1:clust_num) {
    sigma0[[k]]  <- array(NA, dim = c(var_num, var_num, obs_num))
    data_mu[[k]] <- matrix(NA, nrow = obs_num, ncol = var_num)
    for(i in 1:obs_num) {
      data_mu[[k]][i, ]  <- data[i, ] - mu[k, ]
      sigma0[[k]][, , i] <- obs_pp[i, k] * data_mu[[k]][i, ] %*% t(data_mu[[k]][i, ])
    }
    sigma[, , k] <- rowSums(sigma0[[k]], dims = 2) / obs_pp_sums[k]
  }

  return(list(prop = prop,
              mu = mu,
              sigma = sigma
              )
         )
}
