#' Convert +/- Table Labels into Chunklets
#'
#' @param data Dataset in matrix or data.frame form.
#' @param zone_labels Output from label_zones.
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
#' iris_zone_labels <- label_zones(iris[, 1:4], type_marker = iris_tab)$labels
#' iris_chunk_labels <- label_chunklets(iris[, 1:4],
#'                                      zone_labels = iris_zone_labels,
#'                                      zone_percent = 90)
label_chunklets <- function(data, zone_labels, zone_percent) {

  chunk_num <- ncol(zone_labels)
  event_num <- nrow(data)

  if (length(zone_percent) == 1) {
    zone_percent <- rep(zone_percent, times = chunk_num)
  } else if (length(zone_percent) != chunk_num) {
    stop("zone_percent must either be a vector with one entry per chunklet
           or a single value to be used for all chunklets.")
  }

  if (any(zone_percent < 0 | zone_percent > 100)) {
    stop("zone_percent must not be less than 0 or greater than 100")
  } else {
    prob <- (100 - zone_percent) / 100
  }

  zones <- densities <- list()
  cores <- matrix(NA, nrow = event_num, ncol = chunk_num)
  quantiles <- vector("numeric", length = chunk_num)

  # only points in zone l can be assigned to chunklet l
  # the highest Gaussian density points in region l are selected
  for (l in 1:chunk_num) {
    zones[[l]]   <- data[zone_labels[, l], ]

    densities[[l]] <- mvtnorm::dmvnorm(zones[[l]],
                                       mean = colMeans(zones[[l]]),
                                       sigma = stats::cov(zones[[l]]))
    quantiles[l] <- stats::quantile(densities[[l]], prob[l])

    cores[zone_labels[, l], l]  <- densities[[l]] >= quantiles[l]
    cores[!zone_labels[, l], l] <- FALSE
  }

  # zones may overlap but cores are prevented from doing so

  cores_count   <- rowSums(cores)
  cores_overlap <- cores_count > 1
  cores_single  <- cores_count == 1

  if (any(cores_overlap)) {
    find_overlap <- apply(X = cores[cores_overlap, ],
                          MARGIN = 1, simplify = FALSE,
                          FUN = function(x) colnames(zone_labels)[x])
    overlap_names <- unlist(unique(lapply(X = find_overlap,
                                          FUN = function(x) {
                                            paste(x, collapse = "-")
                                            })))

    message(paste0("Initial chunklets overlapped for these populations: ",
                   paste(overlap_names, sep = ", "),
                   ".\n",
                   "Points in intersection were excluded from all chunklets.",
                   "\n \n"))
    }

  core_labels <- chunk_labels <- vector("integer", length = event_num)

  # remove overlap
  core_labels[!cores_single]  <- 0L

  core_labels[cores_single] <- apply(X = cores[cores_single, ],
                                   MARGIN = 1, FUN = which.max)

  chunk_labels[cores_single]  <- core_labels[cores_single]
  chunk_labels[!cores_single] <- chunk_num + 1:sum(!cores_single)

  return(list(core = core_labels,
              chunk = chunk_labels))
}
