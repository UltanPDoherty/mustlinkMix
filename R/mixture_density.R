#' @export
mixture_density <- function(data, params) {
  clust_num <- length(params$prop)
  obs_num <- nrow(data)

  dens_mat <- matrix(nrow = obs_num, ncol = clust_num)
  for (k in seq_len(clust_num)) {
    dens_mat[, k] <- mvtnorm::dmvnorm(
      data,
      mean = params$mu[k, ],
      sigma = params$sigma[, , k]
    )
  }

  dens_vec <- as.numeric(dens_mat %*% params$prop)

  return(dens_vec)
}
