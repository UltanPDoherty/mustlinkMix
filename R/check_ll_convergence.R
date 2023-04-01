#' Compute critical value for log-likelihood convergence test.
#'
#' @param it Iteration number.
#' @param burnin Number of iterations before likelihood convergence criterion is
#'               checked.
#' @param loglike Vector of log-likelihood values after each iteration.
#'
#' @return Critical value for log-likelihood convergence test.
#' @export
#'
#' @examples
#' check_loglike_convergence(2, 0, c(-100, -50))
check_loglike_convergence <- function(it, burnin, loglike) {
  loglike_crit <- NA
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
  }

  return(loglike_crit)
}
