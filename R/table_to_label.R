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
#' out <- table_to_labs(iris_frame, type_marker = iris_tab)
table_to_label <- function(frame, type_marker) {
  obs_num  <- nrow(frame)
  var_num  <- ncol(frame)
  pop_num  <- nrow(type_marker)

  type_marker[type_marker == 0] <- -3
  type_marker <- (type_marker + 1) / 2

  thresholds <- vapply(1:var_num, FUN.VALUE = double(1),
                       FUN = function(p) {
                         flowDensity::deGate(frame, channel = p)
                       }
  )

  obs_bool <- matrix(NA, nrow = obs_num, ncol = var_num)
  obs_orthant <- obs_binary <- c()
  for(i in 1:obs_num) {
    obs_bool[i, ]  <- as.integer(flowCore::exprs(frame)[i, ] > thresholds)
  }

  pops_labs <- list()
  for(j in 1:pop_num){
    nonneutrals  <- type_marker[j, ] != -1
    nn_num <- sum(nonneutrals)
    pops_labs[[j]] <- vapply(X = 1:obs_num,
                             FUN = function(i) {
                               sum(obs_bool[i, nonneutrals] == type_marker[j, nonneutrals]) == nn_num
                             }, FUN.VALUE = logical(1))
  }

  pops_labs2 <- Reduce(f = cbind, x = pops_labs)
  colnames(pops_labs2) <- rownames(type_marker)

  return(list(labs = pops_labs2,
              splits = thresholds))
}
