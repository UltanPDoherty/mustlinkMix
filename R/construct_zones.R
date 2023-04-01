#' @title Construct positive / negative region labels.
#'
#' @description
#' Convert a cell type vs marker table into a logical vector of region
#' membership for each observation.
#'
#' @param data Dataset in matrix or data.frame format.
#' @param type_marker Matrix, rows for populations, columns for variables.
#'
#' @return List containing labels, a matrix with rows for observations and columns
#'         for populations, and splits, a vector of the variable bimodal
#'         threshold values.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' out <- label_zones(iris[, 1:4], type_marker = iris_tab)
construct_zones <- function(data, type_marker) {

  event_num <- nrow(data)
  var_num <- ncol(data)
  zone_num <- nrow(type_marker)

  check_zone_overlap(type_marker = type_marker)

  to_be_split <- which(vapply(X = 1:var_num, FUN.VALUE = logical(1),
                              FUN = function(p) any(type_marker[, p] != 0)))

  # create a vector of +/- splits to check every row of the data against
  splits <- rep(NA, var_num)
  splits[to_be_split] <- vapply(to_be_split, FUN.VALUE = double(1),
                                FUN = function(p) flowDensity::deGate(data[, p]))

  # create a +1/-1 matrix the same size as the data flowFrame
  obs_tab_01 <- t(apply(X = data, MARGIN = 1,
                   FUN = function(vec) vec > splits))
  obs_tab_pm <- 2 * obs_tab_01 - 1

  zone_matrix <- matrix(NA, nrow = event_num, ncol = zone_num,
                      dimnames = list(NULL, zone_names))
  nonneutrals <- type_marker != 0
  for (j in 1:zone_num){
    zone_matrix[, j] <- vapply(X = 1:event_num,
                             FUN = function(i) {
                               all(obs_tab_pm[i, nonneutrals[j, ]]
                                   == type_marker[j, nonneutrals[j, ]])
                             }, FUN.VALUE = logical(1))
  }

  return(list(matrix = zone_matrix,
              splits = splits
              )
         )
}


check_zone_overlap <- function(type_marker) {
  zone_num <- nrow(type_marker)
  zone_names <- rownames(type_marker)

  if (is.null(zone_names)) {
    zone_names <- 1:zone_num
  }

  pairs <- utils::combn(1:zone_num, 2)
  pair_names <- matrix(zone_names[pairs], nrow = 2)

  rowdiffs <- abs(type_marker[pairs[1, ], ] - type_marker[pairs[2, ], ])
  rownames(rowdiffs) <- apply(X = pair_names, MARGIN = 2,
                              FUN = function(x) paste(x, collapse = "-"))

  maxdiffs <- apply(X = rowdiffs, MARGIN = 1, FUN = max)
  if (any(maxdiffs == 0)) {
    stop(paste0("The following pair(s) of population zones are identical: ",
                paste(names(maxdiffs)[which(maxdiffs == 0)], sep = ", "),
                ".\n"))
  }
  if (any(maxdiffs == 1)) {
    message(paste0("The following pair(s) of population zones overlap: ",
                   paste(names(maxdiffs)[which(maxdiffs == 1)], sep = ", "),
                   ".\n"))
  }
}
