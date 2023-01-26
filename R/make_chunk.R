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

  data <- as.matrix(data)

  chunk_num  <- length(unique(chunk_labs))
  var_num    <- ncol(data)

  chunk_size <- as.numeric(table(chunk_labs))
  single_chunk <- chunk_size == 1
  single_obs   <- chunk_size[chunk_labs] == 1

  chunk_mean <- matrix(NA, nrow = chunk_num, ncol = var_num)
  chunk_mean[single_chunk, ] <- data[single_obs, ]
  chunk_sums <- rowsum(data[!single_obs, ], group = chunk_labs[!single_obs])
  chunk_mean[!single_chunk, ] <- chunk_sums / chunk_size[!single_chunk]

  chunk <- list(labs = chunk_labs,
                num  = chunk_num,
                mean = chunk_mean,
                size = chunk_size)
  return(chunk)
}
