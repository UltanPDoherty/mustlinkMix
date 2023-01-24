#' @title Construct positive / negative region labels.
#'
#' @description
#' Convert a cell type vs marker table into a logical vector of region
#' membership for each observation.
#'
#' @param frame flowFrame object
#' @param type_marker Matrix, rows for populations, columns for variables.
#'
#' @return List containing labs, a matrix with rows for observations and columns for populations,
#'         and splits, a vector of the variable bimodal threshold values.
#' @export
#'
#' @examples
#' iris_frame <- flowCore::flowFrame(exprs = as.matrix(iris[, 1:4]))
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' out <- table_to_label(iris_frame, type_marker = iris_tab)
table_to_label <- function(frame, type_marker) {
  obs_num  <- nrow(frame)
  var_num  <- ncol(frame)
  pop_num  <- nrow(type_marker)

  # create a vector of +/- thresholds to check every row of the data against
  thresholds <- vapply(1:var_num, FUN.VALUE = double(1),
                       FUN = function(p) {
                         flowDensity::deGate(frame, channel = p)
                       }
  )

  # create a +1/-1 matrix the same size as the data flowFrame
  obs_tab_01 <- t(apply(X = flowCore::exprs(frame), MARGIN = 1,
                   FUN = function(vec) {vec > thresholds}))
  obs_tab_pm <- 2*obs_tab_01 - 1

  pops_labs <- matrix(NA, nrow = obs_num, ncol = pop_num,
                      dimnames = list(NULL, rownames(type_marker)))
  nonneutrals <- type_marker != 0
  for(j in 1:pop_num){
    pops_labs[, j] <- vapply(X = 1:obs_num,
                             FUN = function(i) {
                               all(obs_tab_pm[i, nonneutrals[j, ]] == type_marker[j, nonneutrals[j, ]])
                             }, FUN.VALUE = logical(1))
  }

  return(list(labs = pops_labs,
              splits = thresholds))
}
