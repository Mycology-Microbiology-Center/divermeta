#' Average redundancy
#'
#' Computes the average redundancy \eqn{\overline{Re}}{mean Re} of every sample: the average,
#' across the units present in the sample, of the functional redundancy \eqn{Re}{Re} of each unit
#' computed as in [redundancy()] over its own subunits only. Only the distances between subunits
#' of the same unit are used, so units that are not comparable with each other (e.g. different
#' enzyme families) can be summarized together, and no unit dominates the index because of its
#' size or abundance.
#'
#' @inheritParams relative.multiplicity
#' @inheritParams redundancy
#' @param normalize Logical (default `FALSE`). If `TRUE`, the redundancy of every unit is divided
#'   by its Simpson's diversity \eqn{S}{S}, giving \eqn{1 - Q / S}{1 - Q / S}, before averaging
#'   (see [redundancy()]).
#'
#' @return Named numeric vector with the average redundancy of every sample. Samples with zero
#'   total abundance get `NA`.
#'
#' @details
#' Let \eqn{Re(f)}{Re(f)} be the functional redundancy of the composition \eqn{f}{f} of a unit in
#' the sample, computed from the relative abundances of its subunits within the unit, and
#' \eqn{N_f}{N_f} its total abundance. Then
#' \deqn{\overline{Re} = \frac{1}{\left|\{f : N_f > 0\}\right|} \sum_{f \,:\, N_f > 0} Re(f),}{mean Re = (1/|{f : N_f > 0}|) sum_{f : N_f > 0} Re(f),}
#' with \eqn{Re(f) = S(f) - Q(f)}{Re(f) = S(f) - Q(f)}, or
#' \eqn{Re(f) = 1 - Q(f) / S(f)}{Re(f) = 1 - Q(f) / S(f)} with `normalize = TRUE`.
#'
#' A unit with a single present subunit has no redundancy: it contributes 0 to the average, also
#' when normalized. Only units present in the sample are averaged; with `include_absent = TRUE`,
#' the average is instead taken over all units, and absent units contribute 0. As in
#' [redundancy()], distances greater than 1 are capped at 1, so the index is never negative.
#'
#' @seealso [redundancy()], [relative.multiplicity()]
#'
#' @export
#'
#' @examples
#' # Two units: A with subunits a1, a2, a3 and B with subunits b1, b2
#' clust <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B")
#' diss <- data.frame(
#'   ID1 = c("a1", "a1", "a2", "b1"),
#'   ID2 = c("a2", "a3", "a3", "b2"),
#'   Distance = c(0.3, 0.8, 0.6, 0.4)
#' )
#' abund <- matrix(
#'   c(5, 3, 0, 2, 0,
#'     0, 1, 4, 2, 6,
#'     2, 0, 0, 0, 1),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("S1", "S2", "S3"), names(clust))
#' )
#'
#' average.redundancy(abund, diss, clust)
#' average.redundancy(abund, diss, clust, normalize = TRUE)
#'
average.redundancy <- function(
  abund,
  diss,
  clust,
  normalize = FALSE,
  include_absent = FALSE,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_flag(normalize, "normalize")
  .check_flag(include_absent, "include_absent")

  units <- .rm_prepare(abund, diss, clust, NULL, chunk_size, check_large_distance_file, header)
  acc <- .average.redundancy_acc(units$ab, units$cl, normalize, include_absent)
  .consume_pairs(units$src, list(acc), units$ab$subunits, units$cl)[[1]]
}


# Internal helpers
# ----------------------------

# Accumulates the average redundancy of every sample over the pairs read. `ab`
# includes the subunits only listed in `clust` (see .extend_abund)
.average.redundancy_acc <- function(ab, cl, normalize, include_absent) {
  units_acc <- .ar_unit_redundancy_acc(ab$A, cl, normalize)

  finish <- function() {
    res <- units_acc$finish()
    .average_units(res$re, res$N > 0, ab, include_absent)
  }

  list(scope = "within", update = units_acc$update, finish = finish)
}


# Redundancy of every unit (columns) in every row of A, reading the pairs from
# `src`. With N = sum a, S2 = sum a^2 and
#   acc = sum over pairs (i, j) of the unit of a_i a_j min(d_ij, 1),
# Simpson's diversity is (N^2 - S2) / N^2 and Rao's Q is 2 acc / N^2, so
#   Re = (N^2 - S2 - 2 acc) / N^2, and normalized Re = 1 - 2 acc / (N^2 - S2)
# N^2 - S2 is exactly 0 when at most one subunit of the unit is present, where
# the redundancy is 0. Returns the redundancies and the unit totals N
.ar_unit_redundancy_acc <- function(A, cl, normalize) {
  n_units <- length(cl$unit_ids)
  acc <- matrix(0, nrow = nrow(A), ncol = n_units)

  update <- function(pairs) {
    within <- cl$unit_of[pairs$i] == cl$unit_of[pairs$j]
    i <- pairs$i[within]
    j <- pairs$j[within]
    d <- pmin(pairs$d[within], 1)
    acc <<- acc + .pair_sums(A, i, j, d, group = cl$unit_of[i], n_groups = n_units)
  }

  finish <- function() {
    G <- .unit_indicator(cl$unit_of, n_units)
    N <- as.matrix(A %*% G)
    S2 <- as.matrix(A^2 %*% G)
    simpson <- N^2 - S2
    re <- if (normalize) 1 - 2 * acc / simpson else (simpson - 2 * acc) / N^2
    re[simpson == 0] <- 0
    list(re = re, N = N)
  }

  list(scope = "within", update = update, finish = finish)
}
