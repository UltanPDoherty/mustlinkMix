#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement M-step of EM algorithm for GMM with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param chunk_num Number of chunklets.
#' @param obs_pp Expanded observation posterior probability matrix.
#' @param chunk_pp Chunklet posterior probability matrix.
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
#' mustlink_mstep_mclust(as.matrix(iris[, 1:4]), chunk_num = chunk1$num,
#'                obs_pp = obs_pp1, chunk_pp = chunk_pp1)
mustlink_mstep_mclust <- function(data, chunk_num, obs_pp, chunk_pp) {

  # Chunklet mixing proportions
  prop <- colSums(chunk_pp) / chunk_num

  mc_params <- mclust::mstepVVV(data = data, z = obs_pp)$parameters

  return(list(prop = prop,
              mu = t(mc_params$mean),
              sigma = mc_params$variance$sigma
  )
  )
}
