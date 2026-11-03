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
#' @param cap_at `NULL` (default) or a single positive number. If `NULL`, the ratios are not
#'   capped. Otherwise, every element's ratio higher than `cap_at` is set to `cap_at`. With
#'   `percentile = 100` every ratio is at most 1, so only a `cap_at` below 1 changes the result.
#'
#' @return Named numeric vector with the proportional abundance of every sample. Samples with
#'   zero total abundance get 0.
#'
#' @details
#' Let \eqn{x_{ij}}{x_ij} be the abundance of element \eqn{j}{j} in sample \eqn{i}{i}, and
#' \eqn{P_j}{P_j} the `percentile` (type 7 quantile, see [stats::quantile()]) of its nonzero
#' abundances across samples. Then
#' \deqn{PA_i = \sum_{j \,:\, x_{ij} > 0} \frac{x_{ij}}{P_j},}{PA_i = sum_{j : x_ij > 0} x_ij / P_j,}
#' or, with a numeric `cap_at` \eqn{c}{c},
#' \deqn{PA_i = \sum_{j \,:\, x_{ij} > 0} \min\left(\frac{x_{ij}}{P_j}, c\right).}{PA_i = sum_{j : x_ij > 0} min(x_ij / P_j, c).}
#'
#' With `percentile = 100`, \eqn{P_j}{P_j} is the maximum abundance of the element, every ratio is
#' in (0, 1] and the index is at most the number of elements present. Elements absent from all
#' samples have no reference level and contribute 0. The index does not change if all abundances
#' are multiplied by the same constant.
#'
#' @seealso [proportional.abundance.ref_abu()] when the reference abundances are known,
#'   [multiplicity.inventory()], [relative.multiplicity()]
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
#' proportional.abundance(abund, percentile = 50, cap_at = 1)
#'
proportional.abundance <- function(abund, percentile = 100, cap_at = NULL) {
  .check_percentile(percentile)
  .check_cap_at(cap_at)
  ab <- .parse_abund(abund)
  .proportional.abundance_core(ab, percentile, cap_at)
}


#' Proportional abundance with known reference abundances
#'
#' Computes the proportional abundance of every sample, as [proportional.abundance()] does, when
#' the reference abundance of every element is already known (e.g. from an external database or
#' a larger collection of samples) instead of being a percentile of the given samples.
#'
#' @inheritParams proportional.abundance
#' @param ref_abu Numeric vector of strictly positive reference abundances, one per element
#'   (column of `abund`). If named, values are matched to the columns by name (extra names are
#'   ignored); otherwise they follow the order of the columns. With `cap_at`, reference abundances
#'   may also be 0: the ratio of such an element is `cap_at` in the samples where it is present,
#'   and 0 where it is absent.
#'
#' @return Named numeric vector with the proportional abundance of every sample. Samples with
#'   zero total abundance get 0.
#'
#' @details
#' With \eqn{R_j}{R_j} the reference abundance of element \eqn{j}{j},
#' \deqn{PA_i = \sum_{j \,:\, x_{ij} > 0} \frac{x_{ij}}{R_j},}{PA_i = sum_{j : x_ij > 0} x_ij / R_j,}
#' or, with a numeric `cap_at` \eqn{c}{c},
#' \deqn{PA_i = \sum_{j \,:\, x_{ij} > 0} \min\left(\frac{x_{ij}}{R_j}, c\right).}{PA_i = sum_{j : x_ij > 0} min(x_ij / R_j, c).}
#'
#' Ratios are not bounded by 1: abundances above their reference give ratios above 1, unless
#' capped. Elements absent from a sample contribute 0 to it.
#'
#' @seealso [proportional.abundance()]
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
#' # Reference abundances from an external database
#' proportional.abundance.ref_abu(abund, ref_abu = c(e1 = 8, e2 = 10, e3 = 5, e4 = 12))
#'
#' # Element e4 without reference abundance: its ratio is cap_at where present
#' proportional.abundance.ref_abu(abund, ref_abu = c(e1 = 8, e2 = 10, e3 = 5, e4 = 0), cap_at = 1)
#'
proportional.abundance.ref_abu <- function(abund, ref_abu, cap_at = NULL) {
  .check_cap_at(cap_at)
  ab <- .parse_abund(abund)
  ref <- .rm_per_unit_values(
    ref_abu, ab$subunits, ab$subunits, "ref_abu",
    allow_zero = !is.null(cap_at),
    zero_hint = " (zero values are allowed only with `cap_at`)",
    what = "element"
  )
  .pa_sum_ratios(ab, .nonzero(ab$A), ref, cap_at)
}


# Proportional abundance of every row of a parsed abundance table (see
# .parse_abund), so divermeta() can pass the abundances of the units
.proportional.abundance_core <- function(ab, percentile, cap_at) {
  nz <- .nonzero(ab$A)
  ref <- if (length(nz$x) > 0) .column_percentiles(nz, ncol(ab$A), percentile)
  .pa_sum_ratios(ab, nz, ref, cap_at)
}


# Adds up x / ref over the nonzero entries `nz` (see .nonzero) of every row,
# capped at cap_at. A zero reference (only allowed with cap_at) gives x / 0 =
# Inf, which becomes cap_at; absent elements are not in `nz` and contribute 0
.pa_sum_ratios <- function(ab, nz, ref, cap_at) {
  sums <- numeric(length(ab$samples))
  if (length(nz$x) > 0) {
    ratio <- nz$x / ref[nz$j]
    if (!is.null(cap_at)) {
      ratio <- pmin(ratio, cap_at)
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
