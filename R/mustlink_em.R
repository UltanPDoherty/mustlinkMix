#' @title Must-Link / Positive Constraint EM GMM.
#'
#' @description
#' Use the EM algorithm to fit a Gaussian Mixture Model that satisfies
#' must-link / positive constraints.
#'
#'
#' @param data Dataset as a data.frame or matrix.
#' @param clust_num Number of components of GMM.
#' @param chunk_labs Vector giving the chunklet label of each observation.
#' @param maxit Maximum number of EM iterations.
#' @param eps Convergence criterion for relative difference in log-likelihood.
#' @param start Initialisation option.
#' @param init_seed Seed.
#'
#' @return A list consisting of a vector of cluster labels,
#'         a matrix of chunklet to cluster assignment probabilities,
#'         a list of model parameters,
#'         a vector of log-likelihood values,
#'         and a vector of times per iteration.
#' @export
#'
#' @examples
#' iris_frame <- flowCore::flowFrame(exprs = as.matrix(iris[, 1:4]))
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_table_labs <- table_to_label(iris_frame, type_marker = iris_tab)$labs
#' iris_chunk_labs <- chunklet_cores(iris[, 1:4], table_labs = iris_table_labs)$chunks
#' mustlink_em(iris[, 1:4], clust_num = 3, chunk_labs = iris_chunk_labs)
mustlink_em <- function(data, clust_num, chunk_labs,
                        maxit = 100, eps = 1e-10, start = "k-Means",
                        init_seed = NULL) {
  data <- as.matrix(data)
  chunk_time <- system.time({
    if (!is.matrix(data)) {
      if (inherits(data, "flowFrame")) {
        data <- flowCore::exprs(data)
      } else {
        data <- as.matrix(data)
      }
    }
    chunk <- make_chunk(data, chunk_labs)
  })
  init_time <- system.time({
    params  <- initialise_model(data, clust_num, start, init_seed)
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

  burnin <- 10

  # EM algorithm
    repeat {
      it <- it + 1

    # E-step
    e_time[it] <- system.time({
      e_out <- mustlink_estep(data, chunk, params,
                            obs_num, var_num, clust_num)
    })[3]

    # mid_time[it] <- system.time({
      ll <- append(ll, e_out$ll)

      # ll_crit is the relative increase in the log-likelihood.
      # if tree accounts for the log-likelihoods being Inf or -Inf.
      if (it >= burnin) {
        ll_diff <- ll[it] - ll[it - 1]
        if (ll[it - 1] == Inf) { ## (Inf, R), (Inf, +Inf), (Inf, -Inf)
          ll_crit <- NA
        } else if (ll[it - 1] == -Inf & ll[it] == -Inf) { ## (-Inf, -Inf)
          ll_crit <- NA
        } else if (ll_diff == Inf) { ## (-Inf, R), (-Inf, +Inf), (R, +Inf)
          ll_crit <- Inf
        }  else { ## (R, R), (R, +Inf)
          ll_crit <- ll_diff / abs(ll[it - 1])
        }
      }
    # })[3]


    # M-step
    m_time[it] <- system.time({
      params <- mustlink_mstep(data, chunk, obs_pp = e_out$obs_pp,
                               chunk_pp = e_out$chunk_pp,
                               obs_num, var_num, clust_num)
      # params <- mustlink_mstep_mclust(data, chunk_num = chunk$num,
      #                                 obs_pp = e_out$obs_pp,
      #                                 chunk_pp = e_out$chunk_pp)
    })[3]
      if (it == maxit) {
        warning(paste0("EM algorithm did not converge before ", maxit, " iterations."))
        if (!no_print) {
          cat(paste0("...EM stopped at ", Sys.time(), "\n"))
        }
        break
      } else if (!is.na(ll_crit) & ll_crit < eps) {
        if (!no_print) {
          cat(paste0("...EM converged at ", Sys.time(), "\n"))
        }
        break
      }

    cat(paste0("...EM-", it, ",\t",
               "e_time: ",   round(e_time[it], digits = 2), ",\t",
               "ll: ",       round(ll[it],     digits = 5), ",\t",
               # "mid_time: ", round(mid_time[it], digits = 2), ",\t",
               "m_time: ",   round(m_time[it], digits = 2), ",\t",
               "\n"
               )
        )
  }
  # End of EM

  label_time <- system.time({
    chunk_to_clust <- apply(X = e_out$chunk_pp, MARGIN = 1, FUN = which.max)
    clust_labs     <- chunk_to_clust[chunk$labs]

  res <- list(clust_labs = clust_labs,
              pp = e_out$chunk_pp,
              params = params,
              ll = ll,
              time = cbind(e_time, mid_time, m_time, e_time + mid_time + m_time)
              )

  return(res)
}
