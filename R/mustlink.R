#' All-inclusive function for must-link GMM.
#'
#' Only requires flowFrame dataset, cell type - marker table, and number of clusters.
#'
#' @param data Dataset in matrix or data.frame format.
#' @param type_marker Matrix with entries +/-/0, rows for populations, columns for variables.
#' @param clust_num Number of clusters / components.
#' @param prob Probability density quantile threshold for chunklet core construction.
#' @param maxit Maximum number of EM iterations.
#' @param eps Convergence criterion for relative difference in log-likelihood.
#' @param start Initialisation option.
#' @param init_seed Seed.
#' @param print_freq Controls how frequently the log-likelihood and time are
#' printed in the EM loop.
#' @param burnin Controls how many loops are completed before testing for
#' likelihood convergence.
#' @param no_print If TRUE, no output is printed.
#'
#' @return A list consisting of a vector of cluster labels,
#'         a matrix of chunklet to cluster assignment probabilities,
#'         a list of model parameters,
#'         a vector of log-likelihood values,
#'         and a vector of times.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' mustlink(iris[, 1:4], iris_tab, 3)

mustlink <- function(data, type_marker, clust_num, prob,
                     maxit = 100, eps = 1e-10, start = "k-Means",
                     init_seed = NULL, print_freq = 10,
                     burnin = 10, no_print = FALSE) {

  if (is.null(type_marker)) {
    chunk_labs <- 1:nrow(data)
  } else {
    table_labs <- table_to_label(data = data,
                                 type_marker = type_marker)$labs

    chunk_labs <- chunklet_cores(data = data,
                                 table_labs = table_labs)$chunk
  }

  mustlink_em(data = data,
              clust_num = clust_num,
              chunk_labs = chunk_labs,
              maxit = maxit, eps = eps, start = start,
              init_seed = init_seed, print_freq = print_freq,
              burnin = burnin, no_print = no_print)
}
