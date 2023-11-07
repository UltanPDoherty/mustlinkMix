trim_extreme <- function(data, alpha = 0.01) {
  extreme_mat <- array(dim = dim(data))
  maxes <- mins <- margins <- c()
  for (i in seq_len(ncol(data))) {
    maxes[i] <- max(data[, i])
    mins[i]  <- min(data[, i])
    margins[i] <- alpha * (maxes[i] - mins[i])
    extreme_mat[, i] <- data[, i] > maxes[i] - margins[i]
    extreme_mat[, i] <- extreme_mat[, i] | data[, i] < mins[i] + margins[i]
  }

  !apply(extreme_mat, 1, any)
}
