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
#' mustlink(iris[, 1:4], type_marker = iris_tab, clust_num = 3, prob = 0.9)

mustlink <- function(data, type_marker, clust_num, prob = 0.9,
                     maxit = 100, eps = 1e-10, start = "k-Means",
                     init_seed = NULL, print_freq = 10,
                     burnin = 10, no_print = FALSE) {

  setup_time <- system.time({
    if (is.null(type_marker)) {
      chunk_labs <- 1:nrow(data)
    } else {
      zone_labs <- label_zones(data = data,
                                   type_marker = type_marker)$labs

      chunklets <- label_chunklets(data = data,
                                   zone_labs = zone_labs,
                                   prob = prob)
      chunk_labs <- chunklets$chunk
      core_labs  <- chunklets$core
    }

    if (!is.matrix(data)) {
      if (inherits(data, "flowFrame")) {
        data <- flowCore::exprs(data)
      } else {
        data <- as.matrix(data)
      }
    }

    #chunk <- make_chunk(data, chunk_labs)

    params  <- initialise_model(data, clust_num = clust_num,
                                start = start, init_seed = init_seed)

    obs_num <- nrow(data)
    var_num <- ncol(params$mu)
  })

  em_time <- system.time({
    em <- mustlink_em(data = data, chunk_labs = chunk_labs,
                      params = params, clust_num = clust_num,
                      maxit = maxit, eps = eps, burnin = burnin,
                      print_freq = print_freq, no_print = no_print)
  })

  label_time <- system.time({
    chunk_to_clust <- apply(X = em$chunk_pp, MARGIN = 1, FUN = which.max)
    clust_labs     <- chunk_to_clust[chunk_labs]
  })

  times <- rbind(setup_time, em_time, label_time)
  rownames(times) <- c("setup", "em", "label")

  res <- list(clust_labs = clust_labs,
              chunk_labs = chunk_labs,
              em = em,
              times = times
              )

}
