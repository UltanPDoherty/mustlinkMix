#' Compute f1 given two vectors of labels
#'
#' @param clust_labs Vector of cluster labels to evaluate.
#' @param true_labs Vector of true / reference labels to compare to.
#' @param exclude_from_true Vector of true labels to be excluded.
#'                           Ignored if NULL.
#' @param prec_vec Logical: Should vectors of precision and recall values,
#'                 with assignment based on f1, be returned?
#'
#' @return List of mean f1, mean precision, mean recall, and f1 matrix.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_out <- mustlink(iris[, 1:4], type_marker = iris_tab,
#'                      clust_num = 3, prob = 0.9)
#' compute_f1(clust_labs = iris_out$clust_labs,
#'            true_labs = iris$Species,
#'            excluded_from_true = "setosa")

compute_f1 <- function(clust_labs, true_labs,
                       exclude_from_true = NULL,
                       prec_rec = FALSE){

  if (!is.null(exclude_from_true)) {
    excluded   <- true_labs %in% exclude_from_true
    true_labs  <- true_labs[!excluded]
    clust_labs <- clust_labs[!excluded]
    if (is.factor(true_labs)) {
      true_labs <- droplevels(true_labs)
    }
    if (is.factor(clust_labs)) {
      clust_labs <- droplevels(clust_labs)
    }
  }

  true_num  <- length(unique(true_labs))
  clust_num <- length(unique(clust_labs))

  tab <- table(true_labs, clust_labs)

  re_mat <- sweep(tab, STATS = rowSums(tab), FUN = "/", MARGIN = 1)
  pr_mat <- sweep(tab, STATS = colSums(tab), FUN = "/", MARGIN = 2)

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

  out <- list(f1_mat = f1_mat2,
              f1_vec = f1_vec)

  if (prec_rec) {
    out$pr_vec <- pr_vec
    out$re_vec <- re_vec
  }

  return(out)
}
