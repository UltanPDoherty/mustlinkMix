#' Implement E-step of EM algorithm with positive constraints.
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
  sigma_inv  <- apply(X = params$sigma, MARGIN = 3, simplify = FALSE,
                      FUN = solve)
  norm_const <- vapply(X = sigma_inv, FUN.VALUE = double(1),
                      FUN = function(x) {
                        (log(det(x)) - var_num * log(2 * pi)) / 2
                        }
                      )

  obs_pdf <- vapply(1:clust_num,
                 FUN = function(k) {
                   norm_const[k] + norm_pdf(data,
                                            params$mu[k, ],
                                            sigma_inv[[k]])
                   },
                 FUN.VALUE = double(obs_num)
                 )

  chunk_pdf <- matrix(NA, chunk$num, clust_num)
  for (l in 1:chunk$num) {
    chunk_pdf[l, ] <- colSums(obs_pdf[chunk$labs == l, , drop = FALSE])
  }

  chunk_pdf     <- sweep(x = chunk_pdf, MARGIN = 2,
                         STATS = log(params$prop), FUN = "+")
  chunk_pdf_max <- max(chunk_pdf)
  chunk_pdf     <- chunk_pdf - chunk_pdf_max

  ll <- chunk$num * chunk_pdf_max + sum(log(rowSums(exp(chunk_pdf))))
  pp <- exp(chunk_pdf) / rowSums(exp(chunk_pdf)) %*% matrix(1, 1, clust_num)

  return(list(ll = ll,
              pp = pp))
}
