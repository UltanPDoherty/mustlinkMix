#' Convert chunk_labels to core_labels, i.e. set all singleton chunklet labels to 0.
#'
#' @param chunk_labels Integer vector of chunklet labels.
#'
#' @return Integer vector of core_labels, same length as chunk_labels.
#' @export
chunk_to_core <- function(chunk_labels) {

  obs_num <- length(chunk_labels)

  chunk_tab <- table(chunk_labels)

  core_labels <- vector("integer", length = obs_num)
  core_labels[chunk_tab[chunk_labels] > 1] <- chunk_labels[chunk_tab[chunk_labels] > 1]
  core_labels[chunk_tab[chunk_labels] == 1] <- 0

  return(core_labels)
}
