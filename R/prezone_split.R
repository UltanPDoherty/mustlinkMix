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

  for (i in seq_along(pz_weights)) {
    if (pz_weights[i] > 1 / 1000) {
      pz_weights[i] <- 1 / 1000
    }
  }

  pz_weights / sum(pz_weights)
}

find_valley <- function(dens, width_percent = 0.01, peak_quantile = 0.75) {
  is_peak <- peak_left <- peak_right <- c()
  w <- round(length(dens$x) * width_percent)
  for (i in seq(w + 1, length(dens$y) - w)) {
    peak_left[i - w] <- all(dens$y[i] > dens$y[(i - w - 1):(i - 1)])
    peak_right[i - w] <- all(dens$y[i] > dens$y[(i + 1):(i+w)])
    is_peak[i - w] <- peak_left[i - w] & peak_right[i - w]
  }

  not_low <- dens$y > quantile(dens$y, peak_quantile)
  is_big_peak <- c(rep(F, w), is_peak, rep(F, w)) & not_low

  if (sum(is_big_peak) == 1) {
    message("Only one peak found.")
    return(NA)
  }

  peak1 <- data.frame(x = dens$x[which.max(dens$y)],
                      y = max(dens$y))
  peak2 <- data.frame(x = dens$x[is_big_peak & dens$y != max(dens$y)],
                      y = dens$y[is_big_peak & dens$y != max(dens$y)])

  scores <- valleys <- c()
  for (i in seq_along(peak2$x)) {
    if (peak2$x[i] < peak1$x) {
      valleys[i] <- min(dens$y[dens$x > peak2$x[i] & dens$x < peak1$x])
    } else {
      valleys[i] <- min(dens$y[dens$x > peak1$x & dens$x < peak2$x[i]])
    }
    scores[i] <- peak2$y[i] - valleys[i]
  }

  return(dens$x[dens$y == valleys[which.max(scores)]])
}

prezone_split <- function(data, splits, type_marker,
                          width_percent = 0.01, peak_quantile = 0.75) {

  pz_labels <- label_prezones(data, splits)
  pz_weights <- weight_prezones(pz_labels)

  nonzero_markers <- lapply(seq_len(ncol(type_marker)),
                            FUN = function(j) type_marker[, j] != 0)
  to_be_split <- vapply(nonzero_markers, FUN.VALUE = logical(1), FUN = any)

  pz_density <- list()
  pz_splits <- splits

  for (j in seq_along(splits)) {
    if (is.na(pz_splits[j]) & to_be_split[j]) {
      cat("prezone_split for marker ", j, "\n")
      pz_density[[j]] <- density(data[, j], weights = pz_weights,
                                 warnWbw = FALSE)
      pz_splits[j] <- find_valley(pz_density[[j]],
                                  width_percent, peak_quantile)
    }
  }

  return(pz_splits)
}
