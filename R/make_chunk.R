#' Compute chunklet number, means, and sizes from labels.
#'
#' @param data Dataset to be clustered.
#' @param chunk_labs Vector of chunklet labels.
#'
#' @return A list.
#' @export
#'
#' @examples
#' make_chunks(iris[, 1:4], c(rep(1, 25), 26:50, rep(51, 25), 76:100, rep(101, 25), 126:150))
make_chunk <- function(data, chunk_labs) {
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
