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
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_zone_matrix <- label_zones(iris[, 1:4], type_marker = iris_tab)$labels
#' iris_chunk_labels <- label_chunklets(iris[, 1:4],
#'                                      zone_matrix = iris_zone_matrix,
#'                                      zone_percent = 90)
label_constraints <- function(data, zone_matrix, zone_percent) {

  zone_num <- ncol(zone_matrix)
  event_num <- nrow(data)

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

    densities[[l]] <- mvtnorm::dmvnorm(zones[[l]],
                                       mean = colMeans(zones[[l]]),
                                       sigma = stats::cov(zones[[l]]))
    quantiles[l] <- stats::quantile(densities[[l]], prob[l])

    linked_set_matrix[zone_matrix[, l], l]  <- densities[[l]] >= quantiles[l]
  }

  # zones may overlap but cores are prevented from doing so
  check_linked_set_overlap(linked_set_matrix)

  linked_set_labels <- chunk_labels <- vector("integer", length = event_num)

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

    unique_overlaps <- unique_rows[check_overlap, ]
    overlap_names <- apply(unique_overlaps, MARGIN = 1,
                           FUN = function(x) {
                             paste(colnames(linked_set_matrix)[x],
                                   collapse = "-")
                           })

    message(paste0("Initial chunklets overlapped for these populations: ",
                   paste(overlap_names, sep = ", "),
                   ".\n",
                   "Points in intersection were excluded from all chunklets.",
                   "\n \n"))
  }
}
