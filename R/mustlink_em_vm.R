#' @title Must-Link / Positive Constraint EM GMM.
#'
#' @param data Dataset in matrix form.
#' @param chunk_labels Output from chunklet_cores.
#' @param params Model parameters, for example, output from initialise_model.
#' @param clust_num Number of clusters.
#' @param burnin Number of iterations before likelihood convergence criterion is
#'               checked.
#' @param no_print If TRUE, no printing takes place.
#' @param print_freq Number of iterations between print statements during EM.
#' @param maxit Maximum number of EM iterations.
#' @param eps Likelihood convergence criterion threshold.
#' @param model Model to be used. Either "vm" for Melnykov et al. or "ns" for
#' Shental et al.
#'
#' @return List of chunklet posterior probability matrix, model parameters, and
#' vector of log-likelihood values for each iteration.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_zone_labels <- label_zones(iris[, 1:4], type_marker = iris_tab)$labels
#' iris_chunk_labels <- label_chunklets(iris[, 1:4],
#'                                    zone_labels = iris_zone_labels,
#'                                    zone_percent = 90)$chunk
#' iris_init <- initialise_model(iris[, 1:4], clust_num = 3,
#'                               start = "k-Means", init_seed = 123)
#' mustlink_em_vm(as.matrix(iris[, 1:4]), clust_num = 3,
#'              chunk_labels = iris_chunk_labels, params = iris_init)

mustlink_em_vm <- function(data, chunk_labels, params, clust_num,
                        burnin = 10, maxit = 1e4, eps = 1e-10,
                        no_print = FALSE, print_freq = 1,
                        model = "vm") {
  it <- 0
  loglike <- c()

  event_num <- nrow(data)
  var_num <- ncol(params$mu)

  chunk <- list(labels = chunk_labels,
                num = length(unique(chunk_labels)),
                size = as.numeric(table(chunk_labels)))

  repeat {
    it <- it + 1

    e_out <- mustlink_estep_vm(data, chunk = chunk, params = params,
                                 obs_num = obs_num, var_num = var_num,
                                 clust_num = clust_num)

    loglike <- append(loglike, e_out$loglike)

    if (!no_print && it %% print_freq == 1) {
      cat(paste0("...No. of E-Steps: ", it,
                 ",\t log-likelihood: ", round(loglike[it],     digits = 5),
                 ",\t Sys.time: ", Sys.time(), "\n"))
    }

    # loglike_crit is the relative increase in the log-likelihood.
    loglike_crit <- compute_loglike_crit(it = it, burnin = burnin,
                                         loglike = loglike)

    # EM has converged if the relative difference between consecutive values
    # of the log-likelihood, i.e. loglike_crit, is not NA and is less than eps.
    if (it == maxit) {
      warning(paste0("EM algorithm did not converge before ",
                     maxit, " iterations."))
      if (!no_print) {
        cat(paste0("...EM stopped at ", Sys.time(), "\n"))
      }
      break
    } else if (!is.na(loglike_crit) && loglike_crit < eps) {
      if (!no_print) {
        cat(paste0("...EM converged at ", Sys.time(), "\n"))
      }
      break
    }

    params <- mustlink_mstep_vm(data,
                                obs_pp = e_out$obs_pp,
                                chunk_pp = e_out$chunk_pp,
                                chunk_num = chunk$num, clust_num = clust_num,
                                obs_num = obs_num, var_num = var_num)
  }

  return(list(chunk_pp = e_out$chunk_pp,
              params = params,
              loglike = loglike))
}
