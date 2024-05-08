#' @title Must-Link / Positive Constraint EM GMM.
#'
#' @description
#' EM algorithm
#'
#' @inheritParams mustlink
#' @param constraints_unique vector of constrained set labels,
#'                           (unconstrained events have unique labels).
#' @param params Model parameters, e.g., output from `initial_parameters`.
#' @param zone_num Number of zones.
#'
#' @return List:
#' * postprob_block:
#' * params:
#' * loglike: vector of log-likelihood values for each iteration
#' @export
mustlink_em <- function(
    data,
    constraints_unique,
    params,
    clust_num,
    zone_num,
    burnin = 2,
    maxit = 1e4,
    eps = 1e-10,
    print_freq = 1,
    drop_cluster = FALSE) {
  it <- 0
  loglike <- c()

  event_num <- nrow(data)
  var_num <- ncol(params$mu)

  constraints_info <- list(
    labels = constraints_unique,
    num = length(unique(constraints_unique)),
    size = as.numeric(table(constraints_unique)),
    zone_num = zone_num
  )

  repeat {
    it <- it + 1

    e_out <- mustlink_estep(
      data,
      constraints_info = constraints_info,
      params = params,
      event_num = event_num,
      var_num = var_num,
      clust_num = clust_num
    )

    loglike <- append(loglike, e_out$loglike)

    if ((it %% print_freq) == 0) {
      cat(paste0(
        format(Sys.time(), "%H:%M:%S"),
        "\t E-Step Number: ", it,
        ",\t Log-likelihood: ", round(loglike[it], digits = 5), "\n"
      ))
    }

    # loglike_crit is the relative increase in the log-likelihood.
    loglike_crit <- compute_loglike_crit(
      it = it,
      burnin = burnin,
      loglike = loglike
    )

    # EM has converged if the relative difference between consecutive values
    # of the log-likelihood, i.e. loglike_crit, is not NA and is less than eps.
    if (it == maxit) {
      warning(paste0(
        "EM algorithm did not converge before ",
        maxit, " iterations."
      ))
      cat(paste0("...EM stopped at ", Sys.time(), "\n"))
      break
    } else if (!is.na(loglike_crit) && loglike_crit < eps) {
      cat(paste0("...EM converged at ", Sys.time(), "\n"))
      break
    }

    params <- mustlink_mstep(
      data,
      postprob_events = e_out$postprob_events,
      postprob_sets = e_out$postprob_sets,
      sets_num = constraints_info$num,
      clust_num = clust_num,
      event_num = event_num,
      var_num = var_num,
      drop_cluster = drop_cluster
    )

    params_na <- all(!is.na(params$mu)) & all(!is.na(params$sigma))

    stopifnot("No NA parameter values" = params_na)

    clust_num <- length(params$prop)
  }

  return(list(
    postprob_sets = e_out$postprob_sets,
    params = params,
    loglike = loglike
  ))
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
    } else { ## (R, R), (R, +Inf)
      loglike_crit <- loglike_diff / abs(loglike[it - 1])
    }
  } else {
    loglike_crit <- NA
  }

  return(loglike_crit)
}
