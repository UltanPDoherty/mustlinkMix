#' Use the EM algorithm to fit a GMM that satisfies positive constraints.
#'
#' @param data Dataset as a data.frame or matrix.
#' @param clust_num Number of components of GMM.
#' @param chunk_labs Vector giving the chunklet label of each observation.
#' @param maxit Maximum number of EM iterations
#'
#' @return A list.
#' @export
#'
#' @examples
#' constr_em2(iris[, 1:4], 3, 1:150)
constr_em2 <- function(data, clust_num, chunk_labs, maxit = 30) {
  # Dataset & its dimensions
  data <- as.matrix(data)
  n    <- nrow(data)
  p    <- ncol(data)

  chunk <- make_chunk(data, chunk_labs)

  # Log-likelihood and posterior probability matrix
  ll <- c(-Inf)
  pp <- matrix(NA, chunk$num, clust_num)

  # Initialise model parameters
  init  <- initialise_model(clust_num, p)
  prop  <- init$prop
  mu    <- init$mu
  sigma <- init$sigma

  # EM algorithm
  for (it in 1:maxit) {
    cat(paste0("...", it))

    # E-step
    for (l in 1:chunk$num) {
      for (k in 1:clust_num) {
        pdf_lk   <- norm_pdf(data[chunk$labs == l, ], mu[k, ], sigma[k, , ])
        like_lk  <- prod(pdf_lk)
        pp[l, k] <- prop[k] * like_lk
      }
    }
    pp <- pp / rowSums(pp) %*% matrix(1, 1, clust_num)

    # M-step
    for (k in 1:clust_num) {
      prop[k] <- sum(pp[, k]) / chunk$num
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
