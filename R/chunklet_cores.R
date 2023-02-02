#' Convert +/- Table Labels into Chunklets
#'
#' @param data Dataset in matrix or data.frame form.
#' @param zone_labs Output from table_to_labs.
#' @param prob Probability cut-off for cores.
#'
#' @return A list of two label vectors, labs with all non-core points labelled 0,
#'         chunks with all non-core points given their own chunklet label.
#' @export
#'
#' @examples
#' iris_tab   <- rbind(se = c(-1, +1, -1, -1),
#'                     ve = c(00, -1, +1, +1),
#'                     vi = c(+1, 00, +1, +1))
#' iris_zone_labs <- label_zones(iris[, 1:4], type_marker = iris_tab)$labs
#' iris_chunk_labs <- label_chunklets(iris[, 1:4], zone_labs = iris_zone_labs)
label_chunklets <- function(data, zone_labs, prob = 0.9) {

  chunk_num <- ncol(zone_labs)
  obs_num   <- nrow(data)
  var_num   <- ncol(data)

  regions <- densities <- list()
  means   <- matrix(NA, nrow = chunk_num, ncol = var_num)
  sigmas  <- array(NA, dim = c(var_num, var_num, chunk_num))
  cores   <- matrix(NA, nrow = obs_num, ncol = chunk_num)
  quants  <- vector("numeric", length = chunk_num)

  # only points in region l can be assigned to chunklet l
  # the highest Gaussian density points in region l are selected
  for(l in 1:chunk_num) {
    regions[[l]]   <- data[zone_labs[, l], ]

    means[l, ]     <- colMeans(regions[[l]])
    sigmas[, , l]  <- stats::cov(regions[[l]])

    densities[[l]] <- mvtnorm::dmvnorm(regions[[l]],
                                       mean = means[l, ],
                                       sigma = sigmas[, , l])
    quants[l] <- stats::quantile(densities[[l]], prob)

    cores[zone_labs[, l], l]  <- densities[[l]] > quants[l]
    cores[!zone_labs[, l], l] <- FALSE
  }

  # regions may not be disjoint, so the cores could overlap, prevent this

  cores_count   <- rowSums(cores)
  cores_overlap <- cores_count > 1
  cores_single  <- cores_count == 1

  if (any(cores_overlap)) {
    find_overlap <- apply(X = cores[cores_overlap, ], MARGIN = 1,
                          FUN = function(x) {colnames(zone_labs)[x]})
    overlap_names <- unique(apply(X = find_overlap, MARGIN = 2,
                                  FUN = function(x) {paste(x, collapse = "-")}))

    message(paste0("Initial cores overlapped for the following pair(s): ",
                   paste(overlap_names, sep = ", "), ".\n",
                   "Points in intersection were excluded from all cores."))
    }

  core_labs <- chunk_labs <- vector("integer", length = obs_num)

  core_labs[!cores_single]  <- 0L # overlap removed

  core_labs[cores_single] <- apply(X = cores[cores_single, ],
                                   MARGIN = 1, FUN = which.max)

  # for(l in 1:chunk_num) {
  #   core_labs[core_labs == l] <- l * zone_labs[core_labs == l, l]
  # }

  chunk_labs[cores_single]  <- core_labs[cores_single]
  chunk_labs[!cores_single] <- chunk_num + 1:sum(!cores_single)

  return(list(core = core_labs,
              chunk = chunk_labs))
}
