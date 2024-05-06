#' Use read.table and readRDS to read files created by write.mustlink.
#'
#' @param file_prefix Desired file_prefix.
#' @param clust_labels Logical: should cluster labels be read from a file?
#' @param init_labels Logical: should initial labels be read from a file?
#' @param block_labels Logical: should block labels be read from a file?
#' @param linked_set_labels Logical: should linked set labels be read from a
#'                          file?
#' @param times Logical: should mustlink runtimes be read from a file?
#' @param em_loglike Logical: should log-likelihood values be read from a file?
#' @param em_postprob_block Logical: should posterior probability matrix be read
#'                          from a file?
#' @param em_params Logical: should model parameters be read from a file?
#'
#' @return Same format as mustlink function.
#' @export

read_mustlink <- function(file_prefix,
                          clust_labels = TRUE, init_labels = TRUE,
                          block_labels = TRUE, linked_set_labels = TRUE,
                          times = TRUE, em_loglike = TRUE,
                          em_postprob_block = TRUE, em_params = TRUE) {

  out <- list()

  if (clust_labels) {
    out$clust_labels <- utils::read.table(
      file = paste0(file_prefix, "_clust_labels.txt")
    )$V1
  }

  if (init_labels) {
    out$init_labels <- utils::read.table(
      file = paste0(file_prefix, "_init_labels.txt")
    )$V1
  }

  if (linked_set_labels) {
    out$linked_set_labels <- utils::read.table(
      file = paste0(file_prefix, "_link_labels.txt")
    )$V1
  }

  if (block_labels) {
    out$block_labels <- utils::read.table(
      file = paste0(file_prefix, "_block_labels.txt")
    )$V1
  }

  if (times) {
    out$times <- utils::read.table(
      file = paste0(file_prefix, "_times.txt"), header = TRUE, row.names = 1
    )
  }

  if (any(c(em_loglike, em_postprob_block, em_params))) {
    out$em <- list()
  }

  if (em_loglike) {
    out$em$loglike <- utils::read.table(
      file = paste0(file_prefix, "_loglike.txt")
    )$V1
  }

  if (em_postprob_block) {
    out$em$postprob_block <- utils::read.table(
      file = paste0(file_prefix, "_postprob_block.txt")
    )
  }

  if (em_params) {
    out$em$params <- readRDS(file = paste0(file_prefix, "_params.rds"))
  }

  return(out)
}
