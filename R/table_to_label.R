#' @title Construct positive / negative region labels.
#'
#' @description
#' Convert a cell type vs marker table into a logical vector of region
#' membership for each observation.
#'
#' @param data Dataset in matrix or data.frame format.
#' @param type_marker Matrix, rows for populations, columns for variables.
#'
#' @return List containing labs, a matrix with rows for observations and columns for populations,
#'         and splits, a vector of the variable bimodal threshold values.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' out <- label_zones(iris[, 1:4], type_marker = iris_tab)
label_zones <- function(data, type_marker) {

  obs_num  <- nrow(data)
  var_num  <- ncol(data)
  zone_num  <- nrow(type_marker)
  zone_names <- rownames(type_marker)

  # check if the regions overlap
  if (is.null(zone_names)) { zone_names <- 1:zone_num}
  pairs <- utils::combn(1:zone_num, 2)
  pair_names <- matrix(zone_names[pairs], nrow = 2)

  rowdiffs <- abs(type_marker[pairs[1, ], ] - type_marker[pairs[2, ], ])
  rownames(rowdiffs) <- apply(X = pair_names, MARGIN = 2,
                              FUN = function(x) {paste(x, collapse = "-")})

  maxdiffs <- apply(X = rowdiffs, MARGIN = 1, FUN = max)
  if (any(maxdiffs == 0)){
    stop(paste0("The following pair(s) of population regions are identical: ",
                  paste(names(maxdiffs)[which(maxdiffs == 0)], sep = ", "), ".\n"))
  }
  if (any(maxdiffs == 1)){
    message(paste0("The following pair(s) of population regions overlap: ",
                  paste(names(maxdiffs)[which(maxdiffs == 1)], sep = ", "), ".\n"))
  }


  # create a vector of +/- thresholds to check every row of the data against
  thresholds <- vapply(1:var_num, FUN.VALUE = double(1),
                       FUN = function(p) {
                         flowDensity::deGate(data[, p])
                       }
  )

  # create a +1/-1 matrix the same size as the data flowFrame
  obs_tab_01 <- t(apply(X = data, MARGIN = 1,
                   FUN = function(vec) {vec > thresholds}))
  obs_tab_pm <- 2*obs_tab_01 - 1

  zone_labs <- matrix(NA, nrow = obs_num, ncol = zone_num,
                      dimnames = list(NULL, zone_names))
  nonneutrals <- type_marker != 0
  for(j in 1:zone_num){
    zone_labs[, j] <- vapply(X = 1:obs_num,
                             FUN = function(i) {
                               all(obs_tab_pm[i, nonneutrals[j, ]] == type_marker[j, nonneutrals[j, ]])
                             }, FUN.VALUE = logical(1))
  }

  return(list(labs = zone_labs,
              splits = thresholds
              )
         )
}
