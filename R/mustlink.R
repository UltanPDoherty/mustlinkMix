#' @title Fit a constrained Gaussian mixture model.
#'
#' @description
#' Given a logical `zone_matrix` indicating which events are in each region for
#' which a constrained set should be constructed, choose `zone_percent`% of
#' those events to include in each constrained set.
#'
#' @param data Dataset in `matrix` or `data.frame` format.
#' @param clust_num Number of clusters / components.
#' @param zone_matrix Logical matrix with a column per zone and a row per event.
#' @param zone_percent Percentage of events in zone to be included in each
#'                     constrained set, either one value for all constrained
#'                     sets or one value per constrained set.
#' @param maxit Maximum number of EM iterations.
#' @param eps Convergence criterion for relative difference in log-likelihood.
#' @param init_method Initialisation option.
#' @param init_labels Initial labels.
#' @param init_seed Seed.
#' @param print_freq Controls how frequently the log-likelihood and time are
#'                   printed in the EM loop.
#' @param burnin Controls how many loops are completed before testing for
#'               likelihood convergence.
#' @param drop_cluster Whether empty clusters should be dropped.
#' @param prob_minimum Remove any event whose component assignment probability
#'                     is less than this threshold and classify as an outlier.
#'
#' @return A list
#' * clust_labels: vector of cluster labels
#' * init_labels: vector of initial labels
#' * constraints_common: vector of constrained set labels,
#'                       (unconstrained events have common label 0).
#' * constraints_unique: vector of constrained set labels,
#'                       (unconstrained events have unique labels).
#' * em: list, output from `mustlink_em`.
#' * times: runtimes for setup, em, and labelling.
#' @export
mustlink <- function(
    data,
    clust_num,
    zone_matrix = NULL,
    zone_percent = 100,
    maxit = 1e4,
    eps = 1e-10,
    init_seed = NULL,
    init_method = c("mlkmpp", "mlkm", "kmpp", "km", "old_mlkm", "old_mlkmpp"),
    init_labels = NULL,
    print_freq = 10,
    burnin = 2,
    drop_cluster = FALSE,
    prob_minimum = 0) {
  init_method <- rlang::arg_match(init_method)

  setup_time <- system.time({
    if (is.null(zone_matrix)) {
      constraints_unique <- constraints_common <- seq_len(nrow(data))
      zone_num <- 0
    } else {
      zone_num <- ncol(zone_matrix)

      constraints <- label_constraints(
        data = data,
        zone_matrix = zone_matrix,
        zone_percent = zone_percent
      )
      constraints_unique <- constraints$unique
      constraints_common <- constraints$common
    }

    if (!is.matrix(data)) {
      if (inherits(data, "flowFrame")) {
        data <- flowCore::exprs(data)
      } else {
        data <- as.matrix(data)
      }
    }

    if (is.null(init_labels)) {
      init_labels <- initial_partition(
        data,
        clust_num = clust_num,
        constraints_common = constraints_common,
        init_seed = init_seed,
        init_method = init_method
      )
    }
    init_params <- initial_parameters(data, init_labels = init_labels)
  })

  em_time <- system.time({
    em <- mustlink_em(
      data = data,
      constraints_unique = constraints_unique,
      params = init_params,
      clust_num = clust_num,
      zone_num = zone_num,
      maxit = maxit,
      eps = eps,
      burnin = burnin,
      print_freq = print_freq,
      drop_cluster = drop_cluster
    )
  })

  label_time <- system.time({
    sets_to_clust <- apply(X = em$postprob_sets, MARGIN = 1, FUN = which.max)
    clust_labels <- sets_to_clust[constraints_unique]
  })

  outliers <- apply(em$postprob_sets, 1, max) < prob_minimum
  clust_labels[outliers] <- 0

  times <- rbind(setup_time, em_time, label_time)
  rownames(times) <- c("setup", "em", "label")

  return(list(
    clust_labels = clust_labels,
    init_labels = init_labels,
    constraints_common = constraints_common,
    constraints_unique = constraints_unique,
    em = em,
    times = times
  ))
}
