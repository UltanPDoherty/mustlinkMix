#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param chunk Object from make_chunk.
#' @param chunk_pp Posterior probability matrix.
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
#' chunk_pp1 <- e_out1$chunk_pp
#' mustlink_mstep(as.matrix(iris[, 1:4]), chunk = chunk1, chunk_pp <- chunk_pp1,
#'              obs_num = 150, var_num = 4, clust_num = 3)
mustlink_mstep <- function(data, chunk, chunk_pp,
                           obs_num = nrow(data),
                           var_num = ncol(data),
                           clust_num = ncol(chunk_pp)) {
 obs_pp <- matrix(NA, nrow = obs_num, ncol = clust_num)
  for(l in 1:chunk$num) {
    obs_pp[chunk$labs == l, ] <- matrix(1, chunk$size[l], 1) %*% t(chunk_pp[l, ])
  }

  obs_pp_sums <- colSums(obs_pp)
  mu      <- t(obs_pp) %*% data / obs_pp_sums

  sigma0 <- data_mu <- list()
  sigma  <- array(NA, dim = c(var_num, var_num, clust_num))

  for(k in 1:clust_num) {
    sigma0[[k]]  <- array(NA, dim = c(var_num, var_num, obs_num))
    data_mu[[k]] <- data - matrix(1, obs_num, 1) %*% t(mu[k, ])
    for(i in 1:obs_num) {
      sigma0[[k]][, , i] <- obs_pp[i, k] * data_mu[[k]][i, ] %*% t(data_mu[[k]][i, ])
    }
    sigma[, , k] <- rowSums(sigma0[[k]], dims = 2) / obs_pp_sums[k]
  }

  prop <- colSums(chunk_pp) / chunk$num

  return(list(prop = prop,
              mu = mu,
              sigma = sigma
              )
         )
}
