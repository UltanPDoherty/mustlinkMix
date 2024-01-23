#' All-inclusive function for must-link GMM.
#'
#' Only requires flowFrame dataset, cell type - marker table, number of
#' clusters, and size of chunklets.
#'
#' @param data Dataset in matrix or data.frame format.
#' @param clust_num Number of clusters / components.
#' @param zone_matrix Logical matrix with a column per zone and a row per event.
#' @param zone_percent Percentage of events in zone to be included in each
#'                     chunklet, either one value for all chunklets or one value
#'                     per chunklet.
#' @param maxit Maximum number of EM iterations.
#' @param eps Convergence criterion for relative difference in log-likelihood.
#' @param init_method Initialisation option.
#' @param init_labels Initial labels.
#' @param init_seed Seed.
#' @param print_freq Controls how frequently the log-likelihood and time are
#' printed in the EM loop.
#' @param burnin Controls how many loops are completed before testing for
#' likelihood convergence.
#' @param model Model to be used. Either "vm" for Melnykov et al. or "ns" for
#' Shental et al.
#'
#' @return A list consisting of a vector of cluster labels,
#'         a matrix of chunklet to cluster assignment probabilities,
#'         a list of model parameters,
#'         a vector of log-likelihood values,
#'         and a vector of times.
#' @export
mustlink <- function(data, clust_num,
                     zone_matrix = NULL, zone_percent = 100,
                     maxit = 1e4, eps = 1e-10, init_seed = NULL,
                     init_method = c("k-Means++", "k-Means",
                                     "Must-Link k-Means++",
                                     "Must-Link k-Means"),
                     init_labels = NULL,
                     print_freq = 10, burnin = 2,
                     model = c("vm", "ns")) {

  model <- rlang::arg_match(model)
  init_method <- rlang::arg_match(init_method)

  setup_time <- system.time({
    if (is.null(zone_matrix)) {
      block_labels <- linked_set_labels <- seq_len(nrow(data))
      zone_num <- 0
    } else {
      zone_num <- ncol(zone_matrix)

      constraints <- label_constraints(data = data, zone_matrix = zone_matrix,
                                       zone_percent = zone_percent)
      block_labels <- constraints$block
      linked_set_labels <- constraints$linked_set
    }

    if (!is.matrix(data)) {
      if (inherits(data, "flowFrame")) {
        data <- flowCore::exprs(data)
      } else {
        data <- as.matrix(data)
      }
    }

    if (is.null(init_labels)) {
      init_labels <- initial_partition(data, clust_num = clust_num,
                                       linked_set_labels = linked_set_labels,
                                       init_seed = init_seed,
                                       init_method = init_method)
    }
    init_params <- initial_parameters(data, init_labels = init_labels)
  })

  em_time <- system.time({
    em <- mustlink_em(data = data, block_labels = block_labels,
                      params = init_params,
                      clust_num = clust_num, zone_num = zone_num,
                      maxit = maxit, eps = eps, burnin = burnin,
                      print_freq = print_freq, model = model)
  })

  label_time <- system.time({
    block_to_clust <- apply(X = em$postprob_block, MARGIN = 1, FUN = which.max)
    clust_labels <- block_to_clust[block_labels]
  })

  times <- rbind(setup_time, em_time, label_time)
  rownames(times) <- c("setup", "em", "label")

  return(list(clust_labels = clust_labels,
              init_labels = init_labels,
              linked_set_labels = linked_set_labels,
              block_labels = block_labels,
              em = em,
              times = times))
}
