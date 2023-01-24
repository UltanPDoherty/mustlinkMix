#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement E-step of EM algorithm for GMM
#' with positive / must-link constraints.
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
                      mvtnorm::dmvnorm(data, log = TRUE,
                                       mean = params$mu[k, ],
                                       sigma = params$sigma[, , k])
                      }
                    )

  log_maxes <- ll_vec <- chunk_unnorm_sums <- vector(mode = "numeric", length = chunk$num)

  lpdf_chunk <- chunk_unnorm <- chunk_pp <- matrix(NA,
                                                   nrow = chunk$num,
                                                   ncol = clust_num)
  obs_pp <- matrix(NA, nrow = obs_num, ncol = clust_num)

  singles_chunk <- chunk$size == 1
  singles_obs   <- chunk$labs %in% (1:chunk$num)[singles_chunk]

  lpdf_chunk[singles_chunk, ]  <- log_pdf[singles_obs, ]
  lpdf_chunk[!singles_chunk, ] <- rowsum(log_pdf[!singles_obs, ],
                                         group = chunk$labs[!singles_obs])

  # log_maxes <- matrixStats::rowMaxs(x = lpdf_chunk)
  #
  # chunk_unnorm <- exp(eachrow(x = lpdf_chunk - log_maxes,
  #                             y = log(params$prop),
  #                             oper = "+")
  #                     )
  # chunk_unnorm_sums <- rowSums(chunk_unnorm)
  #
  # chunk_pp <- chunk_unnorm / chunk_unnorm_sums
  #
  # ll_vec <- log(chunk_unnorm_sums) + log_maxes

  for (l in 1:chunk$num) {
    log_maxes[l]      <- max(lpdf_chunk[l, ])
    chunk_unnorm[l, ] <- exp(lpdf_chunk[l, ] - log_maxes[l] + log(params$prop))

    chunk_unnorm_sums[l] <- sum(chunk_unnorm[l, ])

    chunk_pp[l, ] <- chunk_unnorm[l, ] / chunk_unnorm_sums[l]

    ll_vec[l] <- log(chunk_unnorm_sums[l]) + log_maxes[l]
  }

  ll  <- sum(ll_vec)

  obs_pp <- chunk_pp[chunk$labs, ]


  return(list(ll = ll,
              chunk_pp = chunk_pp,
              obs_pp = obs_pp))
}
