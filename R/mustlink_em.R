#' Fit a GMM that satisfies positive constraints via the EM algorithm.
#'
#' @param data Dataset as a data.frame or matrix.
#' @param clust_num Number of components of GMM.
#' @param chunk_labs Vector giving the chunklet label of each observation.
#' @param maxit Maximum number of EM iterations.
#' @param eps Convergence criterion for relative difference in log-likelihood.
#' @param start Initialisation option.
#'
#' @return A list consisting of a vector of cluster labels,
#'         a matrix of chunklet to cluster assignment probabilities,
#'         a list of model parameters,
#'         a vector of log-likelihood values,
#'         and a vector of times per iteration.
#' @export
#'
#' @examples
#' chunks_of_25 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' mustlink_em(iris[, 1:4], clust_num = 3, chunk_labs = chunks_of_25)
mustlink_em <- function(data, clust_num, chunk_labs,
                        maxit = 100, eps = 1e-10, start = "vanilla") {
  data <- as.matrix(data)
  chunk_time <- system.time({
    chunk <- make_chunk(data, chunk_labs)
  })
  init_time <- system.time({
    params  <- initialise_model(data, clust_num, start)
  })

  cat(paste0("chunk_time: ", chunk_time[3], ",\t",
             "init_time: ", init_time[3], "\n"
             )
      )

  ll <- e_time <- m_time <- mid_time <- c()
  ll_crit <- eps + 1
  it <- 0

  obs_num <- nrow(data)
  var_num <- ncol(params$mu)

  # EM algorithm
  while (ll_crit > eps && it < maxit) {
    it <- it + 1

    # E-step
    e_time[it] <- system.time({
      e_out <- mustlink_estep(data, chunk, params,
                            obs_num, var_num, clust_num)
    })[3]

    mid_time[it] <- system.time({
      ll <- append(ll, e_out$ll)
      pp <- e_out$pp

      if (it == 1) {
        ll_crit <- Inf
      } else if (ll[it - 1] == -Inf) {
        ll_crit <- Inf
      } else {
        ll_crit <- (ll[it] - ll[it - 1]) / abs(ll[it - 1])
      }
    })[3]


    # M-step
    m_time[it] <- system.time({
      params <- mustlink_mstep(data, chunk, pp,
                             obs_num, var_num, clust_num)
    })[3]

    cat(paste0("...EM-", it, ",\t",
               "e_time: ",   round(e_time[it], digits = 2), ",\t",
               "ll: ",       round(ll[it],     digits = 5), ",\t",
               "mid_time: ", round(mid_time[it], digits = 2), ",\t",
               "m_time: ",   round(m_time[it], digits = 2), ",\t",
               "\n"
               )
        )
  }
  # End of EM

  # Cluster labels
  clust_labs <- rep(NA, nrow(data))
  for (l in 1:chunk$num) {
    clust_labs[chunk$labs == l] <- which.max(pp[l, ])
  }

  res <- list(clust_labs = clust_labs,
              pp = pp,
              params = params,
              ll = ll,
              time = cbind(e_time, mid_time, m_time, e_time + mid_time + m_time)
              )

  return(res)
}
