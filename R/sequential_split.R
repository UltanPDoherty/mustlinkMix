#' Convert +/- Table Labels into Chunklets
#'
#' @param x Dataset in matrix or data.frame form.
#' @param typemarker Cell-type marker table.
#' @param min_height Minimum height for a peak to be recognised by find_peaks.
#' @param min_score Minimum score for a split to be returned by find_valley.
#' @param plot Logical value.
#'
#' @return splits, typemarker, subsetter
#' @export
sequential_split <- function(x, typemarker, min_height, plot = TRUE){
  progress <- splits <- scores <- array(dim = dim(typemarker))
  G <- nrow(typemarker)
  P <- ncol(typemarker)
  N <- nrow(x)
  subsetter <- matrix(TRUE, nrow = N, ncol = G)
  paused <- array(FALSE, dim = dim(typemarker))

  for (g in 1:G){
    for (p in 1:P){
      if (typemarker[g, p] != 0) {
        progress[g, p] = FALSE
      }
    }
  }

  row_plot_num <- floor(sqrt(G))
  col_plot_num <- ceiling(G / row_plot_num)
  round_count <- 0

  while (any(!progress, na.rm = TRUE)) {
    round_count <- round_count + 1
    graphics::par(mfrow = c(row_plot_num, col_plot_num))
    for (g in 1:G){
      if (any(paused[g, !is.na(progress[g, ]) & !progress[g, ]])){
        next
      }
      proposals <- matrix(nrow = 2, ncol = P)
      for (p in 1:P){
        if (!is.na(progress[g, p]) & !progress[g, p]){
          proposals[, p] <- find_valley(
            stats::density(x[subsetter[, g], p]),
            score = TRUE,
            min_score = min_score,
            min_height = min_height)
        }
      }

      if (all(is.na(proposals[1, ]))){
        paused[g, !is.na(progress[g, ]) & !progress[g, ]] <- TRUE
        progress[g, !is.na(progress[g, ])] <- TRUE
        next
      } else {
        p_choice <- which.max(proposals[2, ])
        splits[g, p_choice] <- proposals[1, p_choice]
        scores[g, p_choice] <- proposals[2, p_choice]
        progress[g, p_choice] <- TRUE

        plot(stats::density(x[subsetter[, g], p_choice]),
             main = paste0("Round ", round_count,
                           ": g = ", g, ", p = ", p_choice,
                           ", score = ", round(scores[g, p_choice], 3)))
        graphics::abline(v = splits[g, p_choice])

        if (typemarker[g, p_choice] == +1) {
          subsetter[, g] <- subsetter[, g] & x[, p_choice] > splits[g, p_choice]
        } else {
          subsetter[, g] <- subsetter[, g] & x[, p_choice] < splits[g, p_choice]
        }
      }
    }
  }

  if (G > 1) {
    equal_subsets <- matrix(nrow = G, ncol = G)
    is_a_duplicate <- rep(FALSE, G)
    for (g in 1:(G-1)) {
      for (h in (g + 1):G) {
        equal_subsets[g, h] <- all(subsetter[, g] == subsetter[, h])
        if (equal_subsets[g, h]) {
          print(paste0("Failed to distinguish between populations ",
                       g, " & ", h, "."))
          is_a_duplicate[h] <- TRUE
        }
      }
    }
    to_be_deleted <- which(is_a_duplicate)

    subsetter <- subsetter[, -to_be_deleted]
    progress <- progress[-to_be_deleted, ]
    splits <- splits[-to_be_deleted, ]
    scores <- scores[-to_be_deleted, ]
    paused <- paused[-to_be_deleted, ]
    typemarker <- typemarker[-to_be_deleted, ]
    G <- G - sum(is_a_duplicate)
  }

  for (g in 1:G) {
    for (p in 1:P) {
      if (is.na(splits[g, p])) {
        typemarker[g, p] <- 0
      }
    }
  }

  return(list(splits = splits,
              typemarker = typemarker,
              subsetter = subsetter))
}
