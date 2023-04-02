#' @title Must-Link / Positive Constraint EM GMM.
#'
#' @param data Dataset in matrix form.
#' @param block_labels Each event in a particular linked set has the same
#'                     number and every non-linked event has its own number.
#' @param params Model parameters, for example, output from initialise_model.
#' @param clust_num Number of clusters.
#' @param zone_num Number of zones.
#' @param burnin Number of iterations before likelihood convergence criterion is
#'               checked.
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
#' iris_block_labels <- label_chunklets(iris[, 1:4],
#'                                    zone_labels = iris_zone_labels,
#'                                    zone_percent = 90)$chunk
#' iris_init <- initialise_model(iris[, 1:4], clust_num = 3,
#'                               start = "k-Means", init_seed = 123)
#' mustlink_em_vm(as.matrix(iris[, 1:4]), clust_num = 3,
#'              block_labels = iris_block_labels, params = iris_init)

mustlink_em <- function(data, block_labels, params, clust_num, zone_num,
                        burnin = 2, maxit = 1e4, eps = 1e-10,
                        print_freq = 1,
                        model = c("vm", "ns")) {

  model <- rlang::arg_match(model)

  it <- 0
  loglike <- c()

  event_num <- nrow(data)
  var_num <- ncol(params$mu)

  block <- list(labels = block_labels,
                num = length(unique(block_labels)),
                size = as.numeric(table(block_labels)),
                zone_num = zone_num)

  repeat {
    it <- it + 1

    e_out <- mustlink_estep(data, block = block, params = params,
                            event_num = event_num, var_num = var_num,
                            clust_num = clust_num,
                            model = model)

    loglike <- append(loglike, e_out$loglike)

    if ((it %% print_freq) == 1) {
      cat(paste0(format(Sys.time(), "%H:%M:%S"),
                 "\t E-Step Number: ", it,
                 ",\t Log-likelihood: ", round(loglike[it], digits = 5), "\n"))
    }

    # loglike_crit is the relative increase in the log-likelihood.
    loglike_crit <- compute_loglike_crit(it = it, burnin = burnin,
                                         loglike = loglike)

    # EM has converged if the relative difference between consecutive values
    # of the log-likelihood, i.e. loglike_crit, is not NA and is less than eps.
    if (it == maxit) {
      warning(paste0("EM algorithm did not converge before ",
                     maxit, " iterations."))
      cat(paste0("...EM stopped at ", Sys.time(), "\n"))
      break
    } else if (!is.na(loglike_crit) && loglike_crit < eps) {
      cat(paste0("...EM converged at ", Sys.time(), "\n"))
      break
    }

    params <-  mustlink_mstep(data,
                              postprob_event = e_out$postprob_event,
                              postprob_block = e_out$postprob_block,
                              block_num = block$num, clust_num = clust_num,
                              event_num = event_num, var_num = var_num,
                              model = model)
  }

  return(list(block_pp = e_out$postprob_block,
              params = params,
              loglike = loglike))
}




compute_loglike_crit <- function(it, burnin, loglike) {

  # if tree accounts for the log-likelihoods being Inf or -Inf.
  if (it >= burnin) {
    loglike_diff <- loglike[it] - loglike[it - 1]
    if (loglike[it - 1] == Inf) { ## (Inf, R), (Inf, +Inf), (Inf, -Inf)
      loglike_crit <- NA
    } else if (loglike[it - 1] == -Inf && loglike[it] == -Inf) { ## (-Inf, -Inf)
      loglike_crit <- NA
    } else if (loglike_diff == Inf) { ## (-Inf, R), (-Inf, +Inf), (R, +Inf)
      loglike_crit <- Inf
    }  else { ## (R, R), (R, +Inf)
      loglike_crit <- loglike_diff / abs(loglike[it - 1])
    }
  } else {
    loglike_crit <- NA
  }

  return(loglike_crit)
}
