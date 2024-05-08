#' Produce constraint label vectors based on `zone_matrix`.
#'
#' @inheritParams mustlink
#'
#' @return List:
#' * common: vector of constrained set labels,
#'           (unconstrained events have common label 0).
#' * unique: vector of constrained set labels,
#'           (unconstrained events have unique labels).
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
  constraints_matrix <- zone_matrix
  quantiles <- vector("numeric", length = zone_num)

  # only points in zone l can be assigned to chunklet l
  # the highest Gaussian density points in region l are selected
  for (l in 1:zone_num) {
    zones[[l]] <- data[zone_matrix[, l], ]

    densities[[l]] <- mvtnorm::dmvnorm(zones[[l]],
      mean = colMeans(zones[[l]]),
      sigma = stats::cov(zones[[l]])
    )
    quantiles[l] <- stats::quantile(densities[[l]], prob[l])

    constraints_matrix[zone_matrix[, l], l] <- densities[[l]] >= quantiles[l]
  }

  # zones may overlap but cores are prevented from doing so
  if (zone_num == 1) {
    constraints_common <- as.integer(constraints_matrix[, 1])
  } else if (zone_num > 1) {
    check_constraints_overlap(constraints_matrix)
    constraints_common <- apply(
      X = constraints_matrix, MARGIN = 1,
      FUN = function(x) {
        label <- which(x == max(x))
        ifelse(length(label) == 1, label, 0)
      }
    )
  }

  constraints_unique <- constraints_common
  constraints_unique[constraints_common == 0] <-
    zone_num + seq_len(sum(constraints_common == 0))

  return(list(
    common = constraints_common,
    unique = constraints_unique
  ))
}

check_constraints_overlap <- function(constraints_matrix) {
  unique_rows <- unique(constraints_matrix)
  check_overlap <- rowSums(unique_rows) > 1

  if (any(check_overlap)) {
    unique_overlaps <- unique_rows[check_overlap, , drop = FALSE]
    overlap_names <- apply(unique_overlaps,
      MARGIN = 1,
      FUN = function(x) {
        paste(colnames(constraints_matrix)[x],
          collapse = " & "
        )
      }
    )

    message(paste0(
      "Initial constrained sets overlapped for these populations:",
      paste("\n\t", overlap_names, collapse = ",\n")
    ))
    message(paste0(
      "Points in intersections were excluded from all final ",
      "constrained sets."
    ))
  }
}
