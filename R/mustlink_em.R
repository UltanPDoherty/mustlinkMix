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
#' @param print_freq Controls how frequently the log-likelihood and time are
#' printed in the EM loop.
#' @param burnin Controls how many loops are completed before testing for
#' likelihood convergence.
#' @param no_print If TRUE, no output is printed.
#'
#' @return A list consisting of a vector of cluster labels,
#'         a matrix of chunklet to cluster assignment probabilities,
#'         a list of model parameters,
#'         a vector of log-likelihood values,
#'         and a vector of times.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_table_labs <- table_to_label(iris[, 1:4], type_marker = iris_tab)$labs
#' iris_chunk_labs <- chunklet_cores(iris[, 1:4], table_labs = iris_table_labs)$chunks
#' mustlink_em(iris[, 1:4], clust_num = 3, chunk_labs = iris_chunk_labs)
mustlink_em <- function(data, clust_num, chunk_labs,
                        maxit = 100, eps = 1e-10, start = "k-Means",
                        init_seed = NULL, print_freq = 10,
                        burnin = 10, no_print = FALSE) {

  setup_time <- system.time({
    if (!is.matrix(data)) {
      if (inherits(data, "flowFrame")) {
        data <- flowCore::exprs(data)
      } else {
        data <- as.matrix(data)
      }
    }
    chunk <- make_chunk(data, chunk_labs)
    params  <- initialise_model(data, clust_num, start, init_seed)

    ll <- c()
    ll_crit <- NA
    it <- 0

    obs_num <- nrow(data)
    var_num <- ncol(params$mu)
  })[3]

  # EM algorithm
  if (!no_print) {
    cat(paste0("...EM started at ", Sys.time(), "\n"))
  }
  em_time <- system.time({
    repeat {
      it <- it + 1

      e_out <- mustlink_estep(data, chunk, params,
                              obs_num, var_num, clust_num)

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

      if(!no_print & it %% print_freq == 1){
        cat(paste0("...No. of E-Steps: ", it,
                   ",\t log-likelihood: ", round(ll[it],     digits = 5),
                   ",\t Sys.time: ", Sys.time(), "\n"
                   )
            )
      }

      # EM has converged if the relative difference between consecutive values
      # of the log-likelihood, i.e. ll_crit, is not NA and is less than eps.
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

      params <- mustlink_mstep(data, chunk, obs_pp = e_out$obs_pp,
                                 chunk_pp = e_out$chunk_pp,
                                 obs_num, var_num, clust_num)
      }
  })[3]
  if (!no_print) {
    cat(paste0("...No. of E-Steps: ", it,
               ",\t log-likelihood: ", round(ll[it],     digits = 5),
               ",\t EM time: ", round(em_time, digits = 2) , "\n"))
  }

  label_time <- system.time({
    chunk_to_clust <- apply(X = e_out$chunk_pp, MARGIN = 1, FUN = which.max)
    clust_labs     <- chunk_to_clust[chunk$labs]
  })[3]

  times <- c(setup_time, it, em_time, label_time)
  names(times) <- c("setup", "iterations", "em", "label")

  res <- list(clust_labs = clust_labs,
              pp = e_out$chunk_pp,
              params = params,
              ll = ll,
              times = times
              )

  return(res)
}
