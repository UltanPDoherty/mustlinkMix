#' @title Must-Link / Positive Constraint EM GMM.
#'
#' @inheritParams mustlink
#' @param block_labels Each event in a particular linked set has the same
#'                     number and every non-linked event has its own number.
#' @param params Model parameters, for example, output from initialise_model.
#' @param zone_num Number of zones.
#'
#' @return List:
#' * postprob_block:
#' * params:
#' * loglike: vector of log-likelihood values for each iteration
#' @export
mustlink_em <- function(
    data,
    block_labels,
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

  block <- list(
    labels = block_labels,
    num = length(unique(block_labels)),
    size = as.numeric(table(block_labels)),
    zone_num = zone_num
  )

  repeat {
    it <- it + 1

    e_out <- mustlink_estep(
      data,
      block = block,
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
      postprob_event = e_out$postprob_event,
      postprob_block = e_out$postprob_block,
      block_num = block$num,
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
    postprob_block = e_out$postprob_block,
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
