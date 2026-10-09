#' Proportional abundance
#'
#' Computes the proportional abundance of every sample: the abundance of every element in the
#' sample, divided by a reference level of that element across all samples, added up over the
#' elements. The reference level of an element is a percentile of its nonzero abundances in all
#' samples, so every element is measured on its own scale and abundant elements do not dominate
#' the index.
#'
#' The function does not know whether the columns of `abund` are units or subunits: it treats
#' every column as one element.
#'
#' @param abund Numeric matrix, data.frame or sparse `Matrix` of abundances with samples in rows
#'   and elements (units or subunits) in columns. Column names are the element identifiers.
#' @param percentile Single number in (0, 100] (default `100`, the maximum). Percentile of the
#'   nonzero abundances of every element across the samples used as its reference level.
#' @param cap_at_one Logical (default `FALSE`). If `TRUE`, the ratio of every element is capped
#'   at 1. Only relevant when `percentile < 100`, where ratios can exceed 1.
#'
#' @return Named numeric vector with the proportional abundance of every sample. Samples with
#'   zero total abundance get 0.
#'
#' @details
#' Let \eqn{x_{ij}}{x_ij} be the abundance of element \eqn{j}{j} in sample \eqn{i}{i}, and
#' \eqn{P_j}{P_j} the `percentile` (type 7 quantile, see [stats::quantile()]) of its nonzero
#' abundances across samples. Then
#' \deqn{PA_i = \sum_{j \,:\, x_{ij} > 0} \frac{x_{ij}}{P_j},}{PA_i = sum_{j : x_ij > 0} x_ij / P_j,}
#' or, with `cap_at_one = TRUE`,
#' \deqn{PA_i = \sum_{j \,:\, x_{ij} > 0} \min\left(\frac{x_{ij}}{P_j}, 1\right).}{PA_i = sum_{j : x_ij > 0} min(x_ij / P_j, 1).}
#'
#' With `percentile = 100`, \eqn{P_j}{P_j} is the maximum abundance of the element, every ratio is
#' in (0, 1] and the index is at most the number of elements present. Elements absent from all
#' samples have no reference level and contribute 0. The index does not change if all abundances
#' are multiplied by the same constant.
#'
#' @seealso [multiplicity.inventory()], [relative.multiplicity()]
#'
#' @export
#'
#' @examples
#' abund <- matrix(
#'   c(10, 0, 4, 6,
#'     5,  8, 3, 0,
#'     0, 12, 7, 9),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("S1", "S2", "S3"), c("e1", "e2", "e3", "e4"))
#' )
#'
#' proportional.abundance(abund)
#' proportional.abundance(abund, percentile = 50)
#' proportional.abundance(abund, percentile = 50, cap_at_one = TRUE)
#'
proportional.abundance <- function(abund, percentile = 100, cap_at_one = FALSE) {
  .check_percentile(percentile)
  .check_flag(cap_at_one, "cap_at_one")
  ab <- .parse_abund(abund)
  .proportional.abundance_core(ab, percentile, cap_at_one)
}


# Proportional abundance of every row of a parsed abundance table (see
# .parse_abund), so divermeta() can pass the abundances of the units
.proportional.abundance_core <- function(ab, percentile, cap_at_one) {
  n_samples <- length(ab$samples)
  nz <- .nonzero(ab$A)
  sums <- numeric(n_samples)
  if (length(nz$x) > 0) {
    ref <- .column_percentiles(nz, ncol(ab$A), percentile)
    ratio <- nz$x / ref[nz$j]
    if (cap_at_one) {
      ratio <- pmin(ratio, 1)
    }
    by_row <- rowsum(ratio, nz$i)
    sums[as.integer(rownames(by_row))] <- by_row[, 1]
  }
  # Empty samples have no elements to add up: 0
  names(sums) <- ab$samples
  sums
}


# Type 7 quantile (see stats::quantile) of the nonzero entries `nz` (see
# .nonzero) of every column, all columns at once: with the n entries of a column
# sorted, h = 1 + (n - 1) p / 100 and the quantile interpolates between the
# entries floor(h) and ceiling(h). Columns without nonzero entries get NA
.column_percentiles <- function(nz, n_cols, percentile) {
  ord <- order(nz$j, nz$x)
  x <- nz$x[ord]
  n <- tabulate(nz$j, nbins = n_cols)
  offset <- cumsum(n) - n
  present <- n > 0

  index <- 1 + (n[present] - 1) * percentile / 100
  lo <- floor(index)
  h <- index - lo
  x_lo <- x[offset[present] + lo]
  x_hi <- x[offset[present] + ceiling(index)]
  q <- x_lo
  interp <- h > 0 & x_hi != x_lo
  q[interp] <- (1 - h[interp]) * x_lo[interp] + h[interp] * x_hi[interp]

  ref <- rep(NA_real_, n_cols)
  ref[present] <- q
  ref
}


# `percentile` must be a single finite number in (0, 100]
.check_percentile <- function(percentile) {
  if (!is.numeric(percentile) || length(percentile) != 1 || !is.finite(percentile) ||
    percentile <= 0 || percentile > 100) {
    stop("`percentile` must be a single number in (0, 100]")
  }
}
