#' Implement M-step of EM algorithm with positive constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param chunk Object from make_chunk.
#' @param pp Posterior probability matrix.
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
#' e_out1  <- mustlink_estep(as.matrix(iris[, 1:4]), chunk = chunk1, params = params1,
#'                         obs_num = 150, var_num = 4, clust_num = 3)
#' pp1 <- e_out1$pp
#' mustlink_mstep(as.matrix(iris[, 1:4]), chunk = chunk1, pp <- pp1,
#'              obs_num = 150, var_num = 4, clust_num = 3)
mustlink_mstep <- function(data, chunk, pp,
                         obs_num = nrow(data),
                         var_num = ncol(data),
                         clust_num = ncol(pp)) {

  sigma <- array(0, dim = c(var_num, var_num, clust_num))

  weight <- sweep(x = pp, MARGIN = 1, STATS = chunk$size, FUN = "*")
  weight_sums <- colSums(weight)
  prop   <- weight_sums / obs_num
  mu     <- t(weight) %*% chunk$mean / weight_sums

  #chunk_logic <- lapply(1:chunk$num, function(l){chunk$labs == l})
  chunk_list <- lapply(1:chunk$num, function(l){data[chunk$labs == l, , drop = FALSE]})


  for (k in 1:clust_num) {
    sigma_k <- list()
    chunk_mu_k <- list()
    for (l in 1:chunk$num) {
      chunk_mu_k[[l]]   <- chunk_list[[l]] - matrix(1, chunk$size[l], 1) %*% mu[k, ]
      sigma_k[[l]] <- pp[l, k] * t(chunk_mu_k[[l]]) %*% chunk_mu_k[[l]]
    }
    sigma[, , k] <- Reduce("+", sigma_k) / weight_sums[k]
  }

  return(list(prop = prop,
                 mu = mu,
                 sigma = sigma
              )
         )

}
