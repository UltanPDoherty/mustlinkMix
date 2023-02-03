#' Use read.table and readRDS to read files created by write.mustlink.
#'
#' @param file_prefix Desired file_prefix.
#' @param clust_labs Logical: should cluster labels be read from a file?
#' @param chunk_labs Logical: should chunklet labels be read from a file?
#' @param times Logical: should mustlink runtimes be read from a file?
#' @param em_ll Logical: should log-likelihood values be read from a file?
#' @param em_chunk_pp Logical: should posterior probability matrix be read from a file?
#' @param em_params Logical: should model parameters be read from a file?
#'
#' @return Same format as mustlink function.
#' @export
#'
#'

read_mustlink <- function(file_prefix,
                        clust_labs = TRUE, chunk_labs = TRUE, times = TRUE,
                        em_ll = TRUE, em_chunk_pp = TRUE, em_params = TRUE) {

  out <- list()

  if (clust_labs) {
    out$clust_labs <- utils::read.table(file = paste0(file_prefix, "_clust_labs.txt"))$V1
  }

  if (chunk_labs) {
    out$chunk_labs <- utils::read.table(file = paste0(file_prefix, "_chunk_labs.txt"))$V1
  }

  if (times) {
    out$times <- utils::read.table(file = paste0(file_prefix, "_times.txt"),
                                   header = TRUE, row.names = 1)
  }

  if (any(c(em_ll, em_chunk_pp, em_params))) {
    out$em <- list()
  }

  if (em_ll) {
    out$em$ll <- utils::read.table(file =paste0(file_prefix, "_loglike.txt"))$V1
  }

  if (em_chunk_pp) {
    out$em$chunk_pp <- utils::read.table(file = paste0(file_prefix, "_chunk_pp.txt"))
  }

  if (em_params) {
    out$em$params <- readRDS(file = paste0(file_prefix, "_params.rds"))
  }

  return(out)
}
