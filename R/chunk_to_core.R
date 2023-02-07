#' Convert chunk_labs to core_labs, i.e. set all singleton chunklet labels to 0.
#'
#' @param chunk_labs Integer vector of chunklet labels.
#'
#' @return Integer vector of core_labels, same length as chunk_labs.
#' @export
chunk_to_core <- function(chunk_labs) {

  obs_num <- length(chunk_labs)

  chunk_tab <- table(chunk_labs)
  core_num <- sum( chunk_tab > 1)

  core_labs <- vector("integer", length = obs_num)
  core_labs[chunk_tab[chunk_labs] > 1] <- chunk_labs[chunk_tab[chunk_labs] > 1]
  core_labs[chunk_tab[chunk_labs] == 1] <- 0

  return(core_labs)
}
