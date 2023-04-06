#' Convert +/- Table Labels into Chunklets
#'
#' @param data Dataset in matrix or data.frame form.
#' @param zone_matrix Output from construct_zones.
#' @param zone_percent Percentage of events in zone to be included in each
#'                     chunklet, either one value for all chunklets or one value
#'                     per chunklet.
#'
#' @return A list of two label vectors, labels with all non-core points labelled
#'         0, chunks with all non-core points given their own chunklet label.
#' @export
label_constraints <- function(data, zone_matrix, zone_percent) {

  zone_num <- ncol(zone_matrix)

  if (length(zone_percent) == 1) {
    zone_percent <- rep(zone_percent, times = zone_num)
  } else if (length(zone_percent) != zone_num) {
    stop("zone_percent must either be a vector with one entry per chunklet
           or a single value to be used for all chunklets.")
  }

  if (any(zone_percent < 0 | zone_percent > 100)) {
    stop("zone_percent must not be less than 0 or greater than 100")
  } else {
    prob <- (100 - zone_percent) / 100
  }

  zones <- densities <- list()
  linked_set_matrix <- zone_matrix
  quantiles <- vector("numeric", length = zone_num)

  # only points in zone l can be assigned to chunklet l
  # the highest Gaussian density points in region l are selected
  for (l in 1:zone_num) {
    zones[[l]] <- data[zone_matrix[, l], ]

    densities[[l]] <- mclust::dmvnorm(zones[[l]],
                                       mean = colMeans(zones[[l]]),
                                       sigma = stats::cov(zones[[l]]))
    quantiles[l] <- stats::quantile(densities[[l]], prob[l])

    linked_set_matrix[zone_matrix[, l], l]  <- densities[[l]] >= quantiles[l]
  }

  # zones may overlap but cores are prevented from doing so
  check_linked_set_overlap(linked_set_matrix)

  linked_set_labels <- apply(X = linked_set_matrix, MARGIN = 1,
                          FUN = function(x) {
                            label <- which(x == max(x))
                            ifelse(length(label) == 1, label, 0)
                          })

  block_labels <- linked_set_labels
  block_labels[linked_set_labels == 0] <- zone_num +
                                            seq_len(sum(linked_set_labels == 0))

  return(list(linked_set = linked_set_labels,
              block = block_labels))
}

check_linked_set_overlap <- function(linked_set_matrix) {
  unique_rows <- unique(linked_set_matrix)
  check_overlap <- rowSums(unique_rows) > 1

  if (any(check_overlap)) {

    unique_overlaps <- unique_rows[check_overlap, , drop = FALSE]
    overlap_names <- apply(unique_overlaps, MARGIN = 1,
                           FUN = function(x) {
                             paste(colnames(linked_set_matrix)[x],
                                   collapse = " & ")
                           })

    message(paste0("Initial constrained sets overlapped for these populations:\n",
                   paste("\t", overlap_names,
                         collapse = ",\n")))
    message("Points in intersections were excluded from all final constrained sets.")
  }
}
