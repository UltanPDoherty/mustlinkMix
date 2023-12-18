label_prezones <- function(data, splits) {
  split_matrix <- array(dim = dim(data))
  split_num <- sum(!is.na(splits))

  combos <- expand.grid(rep(list(c(T, F)), split_num))
  colnames(combos) <- paste0("split", seq_len(split_num))

  for (j in seq_along(splits)) {
    if (!is.na(splits[j])) {
      split_matrix[, j] <- data[, j] > splits[j]
    }
  }

  combo_match <- matrix(nrow = nrow(split_matrix),
                        ncol = nrow(combos))
  pz_labels <- c()
  for (i in seq_len(nrow(split_matrix))) {
    combo_match[i, ] <- apply(combos, 1,
                              function(x) all(split_matrix[i, !is.na(splits)] == x))
    pz_labels[i] <- which(combo_match[i, ])
  }
  pz_labels
}

weight_prezones <- function(pz_labels) {
  pz_num <- length(unique(pz_labels))
  pz_sizes <- as.numeric(table(pz_labels))

  pz_weights <- 1 / (pz_num * pz_sizes[pz_labels])

  while (max(pz_weights) > 1 / 1000) {
    pz_weights[pz_weights > 1 / 1000] <- 1 / 1000
    pz_weights <- pz_weights / sum(pz_weights)
  }

  return(pz_weights)
}

# find_peaks <- function(dens, width_percent = 0.01, height_percent = 0) {
#
#   # w is the number of density points in width_percent of the range
#   w <- round(length(dens$x) * width_percent)
#
#   is_peak <- peak_left <- peak_right <- c()
#   for (i in seq(w + 1, length(dens$y) - w)) {
#     peak_left[i - w] <- all(dens$y[i] > dens$y[(i - w - 1):(i - 1)])
#     peak_right[i - w] <- all(dens$y[i] > dens$y[(i + 1):(i+w)])
#     is_peak[i - w] <- peak_left[i - w] & peak_right[i - w]
#   }
#
#   min_peak <- height_percent * max(dens$y)
#   is_peak <- c(rep(F, w), is_peak, rep(F, w)) & dens$y > min_peak
#
#   return(is_peak)
# }

find_peaks <- function(dens, width_percent = 0.01, min_height = 0.01) {

  # w is the number of density points in width_percent of the range
  w <- round(length(dens$x) * width_percent)

  is_peak <- peak_left <- peak_right <- c()
  for (i in seq(w + 1, length(dens$y) - w)) {
    peak_left[i - w] <- all(dens$y[i] > dens$y[(i - w - 1):(i - 1)])
    peak_right[i - w] <- all(dens$y[i] > dens$y[(i + 1):(i+w)])
    is_peak[i - w] <- peak_left[i - w] & peak_right[i - w]
  }

  is_peak <- c(rep(F, w), is_peak, rep(F, w)) & dens$y > min_height

  return(is_peak)
}

find_valley <- function(dens, is_peak) {
  peak1_ind <- which.max(dens$y)
  peak1 <- data.frame(x = dens$x[peak1_ind],
                      y = dens$y[peak1_ind])
  peak2_ind <- which(is_peak & dens$y != max(dens$y))
  peak2 <- data.frame(x = dens$x[peak2_ind],
                      y = dens$y[peak2_ind])

  scores <- valleys <- valley_ind <- c()
  interpeak <- matrix(nrow = length(dens$x), ncol = length(peak2$x))
  for (i in seq_along(peak2$x)) {
    if (peak2$x[i] < peak1$x) {
      interpeak[, i] <- dens$x > peak2$x[i] & dens$x < peak1$x
      valley_ind[i] <- peak2_ind[i] + which.min(dens$y[interpeak[, i]])
      valleys[i] <- dens$y[valley_ind[i]]
    } else {
      interpeak[, i] <- dens$x > peak1$x & dens$x < peak2$x[i]
      valley_ind[i] <- peak1_ind + which.min(dens$y[interpeak[, i]])
      valleys[i] <- dens$y[valley_ind[i]]
    }
    scores[i] <- peak2$y[i] - valleys[i]
  }

  return(dens$x[valley_ind[which.max(scores)]])
}

advanced_split <- function(data, type_marker, splits,
                          width_percent = 0.01) {
  var_num <- ncol(data)

  nonzero_markers <- lapply(seq_len(var_num),
                            FUN = function(j) type_marker[, j] != 0)
  to_be_split <- vapply(nonzero_markers, FUN.VALUE = logical(1), FUN = any)

  if (!is.null(splits)) {
    pz_labels <- label_prezones(data, splits)
    pz_weights <- weight_prezones(pz_labels)
  } else {
    pz_weights <- NULL
    splits <- rep(NA, var_num)
  }

  pz_density <- pz_peaks <- list()
  pz_splits <- splits

  for (j in seq_along(splits)) {
    if (is.na(pz_splits[j]) & to_be_split[j]) {
      cat("advanced_split for marker ", j, "\n")
      pz_density[[j]] <- density(data[, j], weights = pz_weights,
                                 warnWbw = FALSE)
      pz_peaks[[j]] <- find_peaks(pz_density[[j]], width_percent)
      if (sum(pz_peaks[[j]]) == 1) {
        message("Only one peak found.")
        pz_splits[j] <- NA
      } else {
        pz_splits[j] <- find_valley(pz_density[[j]], pz_peaks[[j]])
      }
    }
  }

  return(pz_splits)
}
