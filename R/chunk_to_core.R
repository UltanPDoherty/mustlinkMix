#' @title Convert chunk_labels to core_labels.
#'
#' @description
#' Instead of giving each unconstrained event its own chunklet label, label all
#' unconstrained event as belonging to chunklet 0.
#'
#' @param chunk_labels Integer vector of chunklet labels.
#'
#' @return Integer vector of core_labels, same length as chunk_labels.
#' @export

chunk_to_core <- function(chunk_labels) {

  # Table counting number of events with each chunklet label.
  chunk_tab <- table(chunk_labels)

  # Logical vector indicating unconstrained events.
  unconstrained <- chunk_tab[chunk_labels] == 1

  # Relabel all unconstrained event as 0.
  core_labels <- chunk_labels
  core_labels[unconstrained] <- 0

  return(core_labels)
}
