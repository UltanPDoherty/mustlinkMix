#' Use read.table and readRDS to read files created by write.mustlink.
#'
#' @param file_prefix Desired file_prefix.
#' @param clust_labels Logical: should cluster labels be read from a file?
#' @param init_labels Logical: should initial labels be read from a file?
#' @param chunk_labels Logical: should chunklet labels be read from a file?
#' @param times Logical: should mustlink runtimes be read from a file?
#' @param em_loglike Logical: should log-likelihood values be read from a file?
#' @param em_chunk_pp Logical: should posterior probability matrix be read from
#'                    a file?
#' @param em_params Logical: should model parameters be read from a file?
#'
#' @return Same format as mustlink function.
#' @export

read_mustlink <- function(file_prefix,
                          clust_labels = TRUE, init_labels = TRUE,
                          chunk_labels = TRUE, times = TRUE, em_loglike = TRUE,
                          em_chunk_pp = TRUE, em_params = TRUE) {

  out <- list()

  if (clust_labels) {
    out$clust_labels <- utils::read.table(file = paste0(file_prefix,
                                                      "_clust_labels.txt"))$V1
  }

  if (init_labels) {
    out$init_labels <- utils::read.table(file = paste0(file_prefix,
                                                        "_init_labels.txt"))$V1
  }

  if (chunk_labels) {
    out$chunk_labels <- utils::read.table(file = paste0(file_prefix,
                                                      "_chunk_labels.txt"))$V1
  }

  if (times) {
    out$times <- utils::read.table(file = paste0(file_prefix, "_times.txt"),
                                   header = TRUE, row.names = 1)
  }

  if (any(c(em_loglike, em_chunk_pp, em_params))) {
    out$em <- list()
  }

  if (em_loglike) {
    out$em$loglike <- utils::read.table(file = paste0(file_prefix,
                                                 "_loglike.txt"))$V1
  }

  if (em_chunk_pp) {
    out$em$chunk_pp <- utils::read.table(file = paste0(file_prefix,
                                                       "_chunk_pp.txt"))
  }

  if (em_params) {
    out$em$params <- readRDS(file = paste0(file_prefix, "_params.rds"))
  }

  return(out)
}
