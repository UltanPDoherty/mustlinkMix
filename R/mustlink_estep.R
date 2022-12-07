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
  # cat("NaN: ", sum(is.nan(log_pdf)), "\n")
  # cat("Inf: ", sum(log_pdf == Inf) , "\n")
  # cat("-Inf: ", sum(log_pdf == -Inf) , "\n")

  log_pdf_chunk <- matrix(NA, nrow = chunk$num, ncol = clust_num)
  for (l in 1:chunk$num) {
       log_pdf_chunk[l, ] <- colSums(log_pdf[chunk$labs == l, , drop = FALSE])
  }
  # cat("NaN: ", sum(is.nan(log_pdf_chunk)), "\n")
  # cat("Inf: ", sum(log_pdf_chunk == Inf) , "\n")
  # cat("-Inf: ", sum(log_pdf_chunk == -Inf) , "\n")
  # cat(head(log_pdf_chunk), "\n")

  # prop_mat <- (matrix(1, obs_num, 1) %*% params$prop)
  # prop_pdf <- alt_pdf * prop_mat
  # prop_pdf_chunk <- matrix(NA, nrow = chunk$num, ncol = clust_num)
  # for (l in 1:chunk$num) {
  #   prop_pdf_chunk[l, ] <- matrixStats::colProds(prop_pdf[chunk$labs == l, , drop = FALSE])
  # }
  # alt_chunk_pp <- prop_pdf_chunk / (rowSums(prop_pdf_chunk) %*% matrix(1, 1, clust_num))
  #
  # ll <- sum(log(rowSums(prop_pdf_chunk)))
  # chunk_pp <- alt_chunk_pp


  # # create a list (length clust_num) of inverted covariance matrices
  # sigma_inv  <- apply(X = params$sigma, MARGIN = 3, simplify = FALSE,
  #                     FUN = solve)
  #
  # cat("sigma_inv: ", sigma_inv[[1]], "\n")
  # cat("sigma_inv: ", sigma_inv[[2]], "\n")
  # cat("sigma_inv: ", sigma_inv[[3]], "\n")
  #
  # # create a vector (length clust_num) of log normalising constants
  # log_const <- vapply(X = sigma_inv, FUN.VALUE = double(1),
  #                     FUN = function(x) {
  #                       (log(det(x)) - var_num * log(2 * pi)) / 2
  #                       }
  #                     )
  # cat("NaN: ", sum(is.nan(log_const)))


  # create a matrix (dimension obs_num x clust_num) containing the value of each
  # cluster's unnormalised log pdf evaluated at every observation.
  # un_log_pdf <- vapply(1:clust_num,
  #                   FUN = function(k) {unnorm_log_pdf(data,
  #                                                     params$mu[k, ],
  #                                                     sigma_inv[[k]])
  #                     },
  #                   FUN.VALUE = double(obs_num)
  #                   )
  # # log_pdf (dimension obs_num x clust_num) contains normalised log pdf values.
  # log_pdf <- un_log_pdf + (matrix(1, obs_num, 1) %*% log_const)
  # # pdf (dimension obs_num x clust_num) contains normalised pdf values.
  # pdf <- exp(log_pdf)
  # mm_pdf <- pdf %*% params$prop
  # ll <- sum(log(mm_pdf))

  # create a matrix (dimension chunk$num x clust_num) containing the sum of the
  # unnormalised log pdf values within each chunklet for every cluster.
  # un_log_pdf_chunk <- matrix(NA, chunk$num, clust_num)
  # for (l in 1:chunk$num) {
  #   un_log_pdf_chunk[l, ] <- colSums(un_log_pdf[chunk$labs == l, , drop = FALSE])
  # }

  # create a matrix (dimension chunk$num x clust_num) containing the chunklet
  # sizes multiplied by the log normalising constant of the pdf for each cluster.
  # sized_log_const <- (chunk$size %*% t(log_const))

  # # create a matrix (dimension chunk$num x clust_num) containing the mixing
  # # proportions to the power of the chunklet sizes.
  # prop_to_size <- vapply(X = 1:clust_num, FUN.VALUE = double(chunk$num),
  #                        FUN = function(k){
  #                          params$prop[k]^chunk$size
  #                          })

  # pp0, pp1, and pp have dimension chunk$num x clust_num.
  # chunk_pp0 <- exp(sized_log_const + un_log_pdf_chunk)

  prop_mat  <- matrix(1, chunk$num, 1) %*% params$prop
  # log_minmax <- min(matrixStats::rowMaxs(log_pdf_chunk))
  # cat("log_minmax: ", log_minmax, "\n")
  # log_max   <- max(log_pdf_chunk)
  # cat("log_max: ", log_max, "\n")
  # log_min   <- min(log_pdf_chunk)
  # cat("log_min: ", log_min, "\n")
  log_maxes <- matrixStats::rowMaxs(log_pdf_chunk)
  log_max_mat <- log_maxes %*% matrix(1, 1, clust_num)
  # exp(-746) == -Inf and exp(710) == Inf
  # on log scale, I have a range of 1454
  # log_max for hFD_demo is 2169, log_minmax is -103
  # i cannot subtract more than 641

  # const <- -745 - log_minmax
  log_pdf_chunk2 <- log_pdf_chunk - log_max_mat
  # log_pdf_chunk3 <- matrix(NA, nrow = chunk$num, ncol = clust_num)
  # log_pdf_chunk3[log_pdf_chunk2 >= 709] <- 709
  # log_pdf_chunk3[log_pdf_chunk2 < 709] <- log_pdf_chunk2[log_pdf_chunk2 < 709]

  chunk_pp1 <- exp(log_pdf_chunk2 + log(prop_mat))
  # cat("NaN: ", sum(is.nan(chunk_pp1)), "\n")
  # cat("Inf: ", sum(chunk_pp1 == Inf) , "\n")
  # cat("-Inf: ", sum(chunk_pp1 == -Inf) , "\n")

  # cat(sum(is.nan(chunk_pp1)))
  chunk_pp  <- chunk_pp1 / (rowSums(chunk_pp1) %*% matrix(1, 1, clust_num))
  # cat(sum(is.nan(chunk_pp)))

  ll0 <- rowSums(exp(log_pdf_chunk + log(prop_mat) - log_max_mat))
  # cat(max(ll0), "\n")
  ll1 <- log(ll0) + log_maxes
  #cat(ll1)
  ll <- sum(ll1)

  # obs_pp <- matrix(NA, nrow = obs_num, ncol = clust_num)
  # for(l in 1:chunk$num) {
  #   obs_pp[chunk$labs == l, ] <- matrix(1, chunk$size[l], 1) %*% t(chunk_pp[l, ])
  # }

  return(list(ll = ll,
              chunk_pp = chunk_pp))
}
