#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement E-step of EM algorithm for GMM with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param params List containing prop, mu, sigma.
#' @param chunk Object from make_chunk.
#' @param obs_num Number of observations in the dataset.
#' @param var_num Number of varibles in the dataset.
#' @param clust_num Number of clusters pre-specified.
#'
#' @return A list containing a log-likelihood value
#'         and a posterior probability matrix.
#' @export
#'
#' @examples
#' chunks_of_25 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' chunk1 <- make_chunk(iris[, 1:4], chunk_labs = chunks_of_25)
#' params1 <- initialise_model(iris[, 1:4], clust_num = 3)
#' mustlink_estep(as.matrix(iris[, 1:4]), chunk = chunk1, params = params1,
#'              obs_num = 150, var_num = 4, clust_num = 3)
mustlink_estep <- function(data, chunk, params,
                           obs_num = nrow(data), var_num = ncol(params$mu),
                           clust_num = nrow(params$mu)) {

  log_pdf <- vapply(1:clust_num, FUN.VALUE = double(obs_num),
                    FUN = function(k) {
                      mvtnorm::dmvnorm(data, mean = params$mu[k, ], sigma = params$sigma[, , k],
                                       log = TRUE)
                    }
                    )
   log_pdf_chunk <- matrix(NA, nrow = chunk$num, ncol = clust_num)
  for (l in 1:chunk$num) {
       log_pdf_chunk[l, ] <- colSums(log_pdf[chunk$labs == l, , drop = FALSE])
  }

   chunk_pp0 <- exp(sized_log_const + un_log_pdf_chunk)

  prop_mat    <- matrix(1, chunk$num, 1) %*% params$prop
  log_maxes   <- matrixStats::rowMaxs(log_pdf_chunk)
  log_max_mat <- log_maxes %*% matrix(1, 1, clust_num)

  log_pdf_chunk2 <- log_pdf_chunk - log_max_mat

  chunk_pp1 <- exp(log_pdf_chunk2 + log(prop_mat))

  chunk_pp  <- chunk_pp1 / (rowSums(chunk_pp1) %*% matrix(1, 1, clust_num))

  ll0 <- rowSums(exp(log_pdf_chunk + log(prop_mat) - log_max_mat))
  ll1 <- log(ll0) + log_maxes
  ll  <- sum(ll1)

  return(list(ll = ll,
              chunk_pp = chunk_pp))
}
