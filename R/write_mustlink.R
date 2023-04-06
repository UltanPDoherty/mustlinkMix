#' Use write.table and saveRDS to save the output from a mustlink call.
#'
#' @param mustlink_out Output from mustlink function.
#' @param file_prefix Desired file_prefix.
#' @param clust_labels Logical: should cluster labels be written to a file?
#' @param init_labels Logical: should initial labels be written to a file?
#' @param block_labels Logical: should block labels be written to a file?
#' @param linked_set_labels Logical: should linked set labels be written to a
#'                          file?
#' @param times Logical: should runtimes be written to a file?
#' @param em_loglike Logical: should log-likelihood values be written to a file?
#' @param em_postprob_block Logical: should posterior probability matrix be
#'                          written to a file?
#' @param em_params Logical: should model parameters be written to a file?
#'
#' @return Function does not have a return value.
#' @export
#'
write_mustlink <- function(mustlink_out, file_prefix,
                           clust_labels = TRUE, init_labels = TRUE,
                           block_labels = TRUE, linked_set_labels = TRUE,
                           times = TRUE, em_loglike = TRUE,
                           em_postprob_block = TRUE, em_params = TRUE) {

  if (clust_labels) {
    utils::write.table(mustlink_out$clust_labels,
                       file = paste0(file_prefix, "_clust_labels.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (init_labels) {
    utils::write.table(mustlink_out$init_labels,
                       file = paste0(file_prefix, "_init_labels.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (block_labels) {
    utils::write.table(mustlink_out$block_labels,
                       file = paste0(file_prefix, "_block_labels.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (linked_set_labels) {
    utils::write.table(mustlink_out$linked_set_labels,
                       file = paste0(file_prefix, "_link_labels.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (times) {
    utils::write.table(mustlink_out$times,
                       file = paste0(file_prefix, "_times.txt"),
                       row.names = TRUE)
  }

  if (em_loglike) {
    utils::write.table(mustlink_out$em$loglike,
                       file = paste0(file_prefix, "_loglike.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (em_postprob_block) {
    utils::write.table(mustlink_out$em$postprob_block,
                       file = paste0(file_prefix, "_postprob_block.txt"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (em_params) {
    saveRDS(mustlink_out$em$params,
            file = paste0(file_prefix, "_params.rds"))
  }

}
