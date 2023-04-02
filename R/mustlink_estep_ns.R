#' @title E-Step for Must-link Constrained GMM.
#'
#' @description
#' Implement E-step of EM algorithm for GMM
#' with positive / must-link constraints.
#'
#' @param data Dataset being clustered in matrix form.
#' @param params List containing prop, mu, sigma.
#' @param chunk Object from make_chunk.
#' @param event_num Number of observations in the dataset.
#' @param var_num Number of varibles in the dataset.
#' @param clust_num Number of clusters pre-specified.
#'
#' @return A list containing a log-likelihood value
#'         and a posterior probability matrix.
#' @export
#'
#' @examples
#' chunk_labels1 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' chunk1 <- list(labels = chunk_labels1,
#'                num = length(unique(chunk_labels1)),
#'                size = as.numeric(table(chunk_labels1)))
#' params1 <- initialise_model(iris[, 1:4], clust_num = 3)
#' mustlink_estep_ns(as.matrix(iris[, 1:4]), chunk = chunk1, params = params1,
#'                   event_num = 150, var_num = 4, clust_num = 3)
mustlink_estep_ns <- function(data, chunk, params,
                           event_num = nrow(data), var_num = ncol(params$mu),
                           clust_num = nrow(params$mu)) {

  # log_pdf is an event_num x clust_num matrix
  # it is the log pdf for each component evaluated at every point
  log_pdf <- vapply(1:clust_num, FUN.VALUE = double(event_num),
                    FUN = function(k) {
                      mvtnorm::dmvnorm(data, log = TRUE,
                                       mean = params$mu[k, ],
                                       sigma = params$sigma[, , k])
                      }
                    )

  # singles_chunk is a logical of length chunk$num
  # it identifies which chunklets are singletons,
  # i.e. unconstrained observations
  singles_chunk <- chunk$size == 1
  # singles_obs is a logical of length event_num
  # it identifies which observations correspond to the singleton chunklets
  singles_obs   <- chunk$labels %in% (1:chunk$num)[singles_chunk]

  nonsingle_chunknum <- (1:chunk$num)[!singles_chunk]
  single_chunknum   <- chunk$labels[singles_obs]

  # lpdf_chunk is the sum of log_pdf values within each chunklet
  lpdf_chunk <- matrix(NA, nrow = chunk$num, ncol = clust_num)
  lpdf_chunk[single_chunknum, ] <- log_pdf[singles_obs, ]
  lpdf_chunk[nonsingle_chunknum, ] <- rowsum(log_pdf[!singles_obs, ],
                                               group = chunk$labels[!singles_obs])

  # This for loop computes the chunklet posterior probability matrix, chunk_pp,
  chunk_unnorm <- chunk_pp <- matrix(NA, nrow = chunk$num, ncol = clust_num)
  log_maxes <- loglike_vec <- chunk_unnorm_sums <- vector(mode = "numeric",
                                                     length = chunk$num)

  for (l in 1:chunk$num) {

    # Add the log mixing proportions and then un-log this sum with exp.
    # Subtract lpdf_chunk row maxes to prevent exp mapping large values to Inf.
    log_maxes[l]      <- max(lpdf_chunk[l, ])
    chunk_unnorm[l, ] <- exp(lpdf_chunk[l, ] - log_maxes[l] + log(params$prop))

    # Normalise rows of chunk_unnorm to obtain chunk_pp.
    chunk_unnorm_sums[l] <- sum(chunk_unnorm[l, ])
    chunk_pp[l, ] <- chunk_unnorm[l, ] / chunk_unnorm_sums[l]

    loglike_vec[l] <- log(chunk_unnorm_sums[l]) + log_maxes[l]
  }

  loglike  <- sum(loglike_vec)

  obs_pp <- chunk_pp[chunk$labels, ]

  return(list(loglike = loglike,
              chunk_pp = chunk_pp,
              obs_pp = obs_pp))
}
