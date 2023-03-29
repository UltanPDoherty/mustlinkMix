#' Compute critical value for log-likelihood convergence test.
#'
#' @param it Iteration number.
#' @param burnin Number of iterations before likelihood convergence criterion is
#'               checked.
#' @param ll Vector of log-likelihood values after each iteration.
#'
#' @return Critical value for log-likelihood convergence test.
#' @export
#'
#' @examples
#' check_ll_convergence(2, 0, c(-100, -50))
check_ll_convergence <- function(it, burnin, ll) {
  ll_crit <- NA
  # if tree accounts for the log-likelihoods being Inf or -Inf.
  if (it >= burnin) {
    ll_diff <- ll[it] - ll[it - 1]
    if (ll[it - 1] == Inf) { ## (Inf, R), (Inf, +Inf), (Inf, -Inf)
      ll_crit <- NA
    } else if (ll[it - 1] == -Inf && ll[it] == -Inf) { ## (-Inf, -Inf)
      ll_crit <- NA
    } else if (ll_diff == Inf) { ## (-Inf, R), (-Inf, +Inf), (R, +Inf)
      ll_crit <- Inf
    }  else { ## (R, R), (R, +Inf)
      ll_crit <- ll_diff / abs(ll[it - 1])
    }
  }

  return(ll_crit)
}
