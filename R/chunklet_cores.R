#' Convert +/- Table Labels into Chunklets
#'
#' @param data Dataset in matrix form.
#' @param table_labs Output from table_to_labs.
#'
#' @return A list of two label vectors, labs with all non-core points labelled 0,
#'         chunks with all non-core points given their own chunklet label.
#' @export
#'
#' @examples
#' iris_frame <- flowCore::flowFrame(exprs = as.matrix(iris[, 1:4]))
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_table_labs <- table_to_labs(iris_frame, type_marker = iris_tab)
#' iris_chunk_labs <- chunklet_cores(iris[, 1:4], table_labs = iris_table_labs)
chunklet_cores <- function(data, table_labs) {
  chunk_num <- ncol(table_labs)
  obs_num   <- nrow(data)

  regions <- list()
  means <- list()
  sigmas <- list()
  densities <- list()
  cores <- list()
  for(l in 1:chunk_num) {
    regions[[l]]   <- data[table_labs[, l],]
    means[[l]]     <- colMeans(regions[[l]])
    sigmas[[l]]    <- stats::cov(regions[[l]])
    densities[[l]] <- mvtnorm::dmvnorm(data,
                                     mean = means[[l]],
                                     sigma = sigmas[[l]])
    cores[[l]] <- densities[[l]] > stats::quantile(densities[[l]], 0.95)
  }

  core_mat <- Reduce(f = cbind, x = cores)
  core_mat_sums <- rowSums(core_mat)
  non_cores <- core_mat_sums == 0
  if (max(core_mat_sums) > 1) {warning("cores overlap")}

  core_labs <- c()
  core_labs[non_cores] <- 0
  core_labs[!non_cores] <- apply(X = core_mat[!non_cores, ], MARGIN = 1, FUN = which.max)

  for(l in 1:chunk_num) {
    core_labs[core_labs == l] <- l * table_labs[core_labs == l, l]
  }

  core_chunks <- core_labs
  core_chunks[core_chunks == 0] <- (chunk_num + 1):(chunk_num + sum(core_chunks == 0))

  return(list(labs = core_labs,
              chunks = core_chunks))
}
