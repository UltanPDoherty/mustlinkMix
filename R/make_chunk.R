#' @title Construct chunklet for positive / must-link constraints.
#'
#' @description
#' Compute chunklet number, means, and sizes from labels.
#'
#' @param data Dataset to be clustered.
#' @param chunk_labs Vector of chunklet labels.
#'
#' @return A list consisting of chunklet labels, number of chunklets,
#'         a matrix of chunklet means, and a vector of chunklet sizes.
#' @export
#'
#' @examples
#' chunks_of_25 <- c(rep(1, 25), 2:26, rep(27, 25), 28:52, rep(53, 25), 54:78)
#' make_chunk(iris[, 1:4], chunk_labs = chunks_of_25)
make_chunk <- function(data, chunk_labs) {
  data       <- as.matrix(data)
  chunk_labs <- as.numeric(as.factor(chunk_labs))
  chunk_num  <- length(unique(chunk_labs))
  chunk_mean <- matrix(NA, chunk_num, ncol(data))
  chunk_size <- rep(NA, chunk_num)
  for (l in 1:chunk_num) {
    chunk_size[l] <- sum(chunk_labs == l)
    if (chunk_size[l] > 1) {
      chunk_mean[l, ] <- colMeans(data[chunk_labs == l, ])
    } else {
      chunk_mean[l, ] <- data[chunk_labs == l, ]
    }
  }

  chunk <- list(labs = chunk_labs,
                num  = chunk_num,
                mean = chunk_mean,
                size = chunk_size)
  return(chunk)
}
