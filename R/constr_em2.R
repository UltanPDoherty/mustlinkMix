#' Use the EM algorithm to fit a GMM that satisfies positive constraints.
#'
#' @param data Dataset as a data.frame or matrix.
#' @param clust_num Number of components of GMM.
#' @param chunk_labs Vector giving the chunklet label of each observation.
#' @param maxit Maximum number of EM iterations.
#' @param eps Convergence criterion for relative difference in log-likelihood.
#'
#' @return A list.
#' @export
#'
#' @examples
#' constr_em2(iris[, 1:4], 3, 1:150)
constr_em2 <- function(data, clust_num, chunk_labs, maxit = 100, eps = 1e-10, start = "vanilla") {
  # Dataset & its dimensions
  data <- as.matrix(data)
  n    <- nrow(data)
  p    <- ncol(data)

  chunk <- make_chunk(data, chunk_labs)

  # Log-likelihood and posterior probability matrix
  ll <- c(-Inf)
  weighted_pdf <- matrix(NA, chunk$num, clust_num)

  # Initialise model parameters
  init  <- initialise_model(clust_num, p)
  prop  <- init$prop
  mu    <- matrix(rep(colMeans(data), clust_num),
                  nrow = clust_num, byrow = FALSE)
  #mu    <- matrix(rep(colMeans(data), clust_num),
  #                nrow = clust_num, byrow = FALSE) + diag(stats::cov(data)) * init$mu
  sigma <- init$sigma

  ll_crit <- 1
  it <- 0
  # EM algorithm
  while (ll_crit > eps & it < maxit) {
    it <- it + 1

    # E-step
    for (l in 1:chunk$num) {
      for (k in 1:clust_num) {
        pdf_lk   <- norm_pdf(data[chunk$labs == l, ], mu[k, ], sigma[k, , ])
        like_lk  <- prod(pdf_lk)
        weighted_pdf[l, k] <- prop[k] * like_lk
      }
    }

    ll  <- append(ll, log(prod(rowSums(weighted_pdf))))
    if (ll[length(ll) - 1] != -Inf){
      ll_crit <- (ll[length(ll)] - ll[length(ll) - 1])/abs(ll[length(ll) - 1])
    } else {
      ll_crit <- Inf
    }

    pp <- weighted_pdf / rowSums(weighted_pdf) %*% matrix(1, 1, clust_num)

    # M-step
    prop <- chunk$size %*% pp / n
    for (k in 1:clust_num) {
      # prop[k] <- sum(pp[, k]) / chunk$num
      weight  <- pp[, k] * chunk$size
      mu[k, ] <- colSums(weight * chunk$mean) / sum(weight)
      sigma_k <- matrix(0, p, p)
      for (l in 1:chunk$num) {
        chunk_l    <- data[chunk$labs == l, ]
        centred_lk <- chunk_l - matrix(1, chunk$size[l], 1) %*% mu[k, ]
        sigma_k    <- sigma_k + pp[l, k] * t(centred_lk) %*% centred_lk
      }
      sigma[k, , ] <- sigma_k / sum(pp[, k] * chunk$size)
    }

    cat(paste0("...", it))
  }
  # End of EM
  cat("\n")

  # Cluster labels
  clust_labs <- rep(NA, n)
  for (l in 1:chunk$num) {
    clust_labs[chunk$labs == l] <- which.max(pp[l, ])
  }

  res <- list(clust_labs = clust_labs,
              pp = pp,
              prop = prop,
              mu = mu,
              sigma = sigma,
              ll = ll)

  return(res)
}
