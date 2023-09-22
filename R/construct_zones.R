#' @title Construct positive / negative region labels.
#'
#' @description
#' Convert a cell type vs marker table into a logical vector of region
#' membership for each observation.
#'
#' @param data Dataset in matrix or data.frame format.
#' @param type_marker Matrix, rows for populations, columns for variables.
#' @param custom_splits Vector of values for variables' bimodal thresholds.
#'
#' @return List containing
#' \itemize{
#' \item labels: a matrix of logical values with a row per observation and a
#'               column per population.
#' \item splits: a data frame containing upper, lower, and default bimodal
#'               threshold values for each marker from flowDensity::deGate.
#' }
#' @export
construct_zones <- function(data, type_marker, custom_splits = NULL) {

  event_num <- nrow(data)
  zone_num <- nrow(type_marker)

  if (zone_num > 1) {
    check_zone_overlap(type_marker = type_marker)
  }

  if (is.null(custom_splits)) {
    splits <- compute_splits(data, type_marker = type_marker)
  } else {
    splits <- custom_splits
  }

  zone_matrix <- construct_zone_matrix(data, type_marker, splits)

  zone_sizes <- colSums(zone_matrix)
  if (any(zone_sizes == 0)) {
    stop(paste0("The following zones do not contain any events: ",
                paste0(colnames(zone_matrix)[zone_sizes == 0],
                       collapse = ", "), "\n"))
  }

  return(list(matrix = zone_matrix,
              splits = splits))
}

construct_zone_matrix <- function(data, type_marker, splits) {

  event_num <- nrow(data)
  zone_num <- nrow(type_marker)
  zone_names <- rownames(type_marker)

  # create a +1/-1 matrix the same size as the data flowFrame
  is_event_positive <- t(apply(X = data, MARGIN = 1,
                               FUN = function(event) event > splits))
  event_plusminus <- (2 * is_event_positive) - 1

  zone_matrix <- matrix(NA, nrow = event_num, ncol = zone_num,
                        dimnames = list(NULL, zone_names))
  nonzero_markers <- list()
  for (j in 1:zone_num){
    nonzero_markers[[j]] <- type_marker[j, ] != 0
    zone_matrix[, j] <- apply(
      event_plusminus,
      MARGIN = 1,
      function(x) {
        all(x[nonzero_markers[[j]]] == type_marker[j, nonzero_markers[[j]]])
      }
    )
  }

  return(zone_matrix)
}

check_zone_overlap <- function(type_marker) {
  zone_num <- nrow(type_marker)
  zone_names <- rownames(type_marker)

  if (is.null(zone_names)) {
    zone_names <- 1:zone_num
  }

  pair_nums  <- utils::combn(1:zone_num, 2)
  pair_names <- utils::combn(zone_names, 2,
                             FUN = function(x) paste(x, collapse = " & "))

  rowdiffs <- abs(type_marker[pair_nums[1, ], ] - type_marker[pair_nums[2, ], ])

  if (zone_num == 2) {
    maxdiffs <- max(rowdiffs)
  } else {
    maxdiffs <- apply(X = rowdiffs, MARGIN = 1, FUN = max)
  }
  if (any(maxdiffs == 0)) {
    stop(paste0("The following population zones are identical:\n",
                paste("\t", pair_names[which(maxdiffs == 0)],
                      collapse = ",\n")))
  }
  if (any(maxdiffs == 1)) {
    message(paste0("The following population zones overlap:\n",
                   paste("\t", pair_names[which(maxdiffs == 1)],
                         collapse = ",\n")))
  }
}

#' @title Use flowDensity::deGate to split the markers.
#'
#' @param data Dataset in matrix or data.frame format.
#' @param type_marker Matrix, rows for populations, columns for variables.
#'
#' @return splits: a vector of bimodal threshold values for each marker
#' from flowDensity::deGate.
#'
#' @export
compute_splits <- function(data, type_marker) {
  var_num <- ncol(data)

  nonzero_markers <- lapply(seq_len(var_num),
                            FUN = function(j) type_marker[, j] != 0)
  to_be_split <- vapply(nonzero_markers, FUN.VALUE = logical(1), FUN = any)

  # create a vector of +/- splits to check every row of the data against
    }
  }

  return(splits)
}
