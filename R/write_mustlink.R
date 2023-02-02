#' Use write.table and saveRDS to save the output from a mustlink call.
#'
#' @param mustlink_out Output from mustlink function.
#' @param file_prefix Desired file_prefix.
#' @param clust_labs Logical: should cluster labels be written to a file?
#' @param chunk_labs Logical: should chunklet labels be written to a file?
#' @param times Logical: should runtimes be written to a file?
#' @param em_ll Logical: should log-likelihood values be written to a file?
#' @param em_chunk_pp Logical: should posterior probability matrix be written to a file?
#' @param em_params Logical: should model parameters be written to a file?
#'
#' @return Function does not have a return value.
#' @export
#'
write_mustlink <- function(mustlink_out, file_prefix,
                           clust_labs = TRUE, chunk_labs = TRUE, times = TRUE,
                           em_ll = TRUE, em_chunk_pp = TRUE, em_params = TRUE) {

  if (clust_labs) {
    utils::write.table(mustlink_out$clust_labs,
                       file = paste0(file_prefix, "_clust_labs.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (chunk_labs) {
    utils::write.table(mustlink_out$chunk_labs,
                       file = paste0(file_prefix, "_chunk_labs.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (times) {
    utils::write.table(mustlink_out$times,
                       file = paste0(file_prefix, "_times.txt"),
                       row.names = TRUE)
  }

  if (em_ll) {
    utils::write.table(mustlink_out$em$ll,
                       file = paste0(file_prefix, "_loglike.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (em_chunk_pp) {
    utils::write.table(mustlink_out$em$chunk_pp,
                       file = paste0(file_prefix, "_chunk_pp.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (em_params) {
    saveRDS(mustlink_out$em$params,
            file = paste0(file_prefix, "_params.rds"))
  }

}
