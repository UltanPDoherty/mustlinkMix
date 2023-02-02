#' Compute f1 given two vectors of labels
#'
#' @param clust_labs Vector of cluster labels to evaluate.
#' @param true_labs Vector of true / reference labels to compare to.
#'
#' @return List of mean f1, mean precision, mean recall, and f1 matrix.
#' @export
#'
#' @examples
compute_f1 <- function(clust_labs, true_labs){

  true_num  <- length(unique(true_labs))
  clust_num <- length(unique(clust_labs))

  tab <- table(true_labs, clust_labs)

  re_mat <- sweep(tab, STAT = rowSums(tab), FUN = "/", MARGIN = 1)
  pr_mat <- sweep(tab, STAT = colSums(tab), FUN = "/", MARGIN = 2)

  f1_mat <- 2 * re_mat * pr_mat / (re_mat + pr_mat)
  f1_mat[f1_mat == "NaN"] <- 0

  if (true_num <= clust_num) {
    true_to_clust <- as.numeric(clue::solve_LSAP(f1_mat[1:true_num, ], maximum = TRUE))
    other_clust   <- setdiff(1:clust_num, true_to_clust)
    true_to_clust <- c(true_to_clust, other_clust)
  } else {
    clust_to_true <- as.numeric(clue::solve_LSAP(t(f1_mat[1:true_num, ]), maximum = TRUE))
    other_true    <- setdiff(1:true_num, clust_to_true)
    true_to_clust <- (1:true_num)[order(c(clust_to_true, other_true))]
    true_to_clust[true_to_clust > clust_num] <- NA
  }

  f1_mat2 <- f1_mat[, true_to_clust]
  f1_mat2[, is.na(true_to_clust)] <- 0

  f1_vec <- pr_vec <- re_vec <- rep(NA, true_num)
  for (i in 1:true_num) {
    if (!is.na(true_to_clust)[i]) {
      f1_vec[i] <- f1_mat[i, true_to_clust[i]]
      pr_vec[i] <- pr_mat[i, true_to_clust[i]]
      re_vec[i] <- re_mat[i, true_to_clust[i]]
    } else {
      f1_vec[i] <- 0
      pr_vec[i] <- 0
      re_vec[i] <- 0
    }
  }

  matched_vals <- c(sum(f1_vec), sum(pr_vec), sum(re_vec)) / sum(!is.na(true_to_clust))
  names(matched_vals) <- c("f1", "pr", "re")

  return(list(mean_f1 = mean(f1_vec),
              mean_pr = mean(pr_vec),
              mean_re = mean(re_vec),
              f1_mat = f1_mat2
              )
         )
}
