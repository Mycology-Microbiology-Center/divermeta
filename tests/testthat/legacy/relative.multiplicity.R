# LEGACY implementation (pre sample x subunit / three-column diss migration).
# Not part of the package: sourced by tests/testthat/helper-legacy.R and used
# as a reference in the tests. Do not modify.

#' Relative multiplicity
#'
#' Computes relative multiplicity \eqn{RM}{RM} for a single sample: the average, across the units
#' present in the sample, of each unit's distance-based functional diversity \eqn{\delta D_{\sigma}}{delta D_sigma} relative
#' to the diversity of a reference composition of that same unit. Because every unit is compared
#' only against itself, units that are not comparable with each other (e.g. different enzyme
#' families) can be summarized together, and no unit dominates the index because of its size or
#' abundance.
#'
#' @param ab List of numeric vectors, one per unit, with the abundances of the unit's subunits in
#'   the sample. A single numeric vector is treated as a sample with one unit.
#' @param diss List of numeric matrices or `dist` objects, one per unit, with the pairwise
#'   dissimilarities among the unit's subunits. `diss[[i]]` must be of size
#'   `length(ab[[i]])` and in the same order.
#' @param ab_ref List of numeric vectors, one per unit, with the abundances of the subunits in the
#'   unit's reference composition.
#' @param diss_ref List of numeric matrices or `dist` objects, one per unit, with the pairwise
#'   dissimilarities among the reference subunits. `diss_ref[[i]]` must be of size
#'   `length(ab_ref[[i]])` and in the same order.
#' @param sigma Numeric cutoff \eqn{\sigma}{sigma} (default `1`) at which two subunits are considered
#'   maximally different. Either a single value, applied to every unit, or a vector of strictly
#'   positive values with one value per unit.
#' @param cap_at_one Logical (default `FALSE`). If `TRUE`, each unit's ratio is capped at 1, for
#'   the cases where the unit is more diverse in the sample than in its reference.
#' @param include_absent Logical (default `FALSE`). If `FALSE`, only the units present in the
#'   sample are averaged. If `TRUE`, every unit counts in the average, and units absent from the
#'   sample contribute 0.
#'
#' @return Numeric scalar, the relative multiplicity \eqn{RM}{RM}. `NA` if no unit has abundance
#'   in the sample and `include_absent = FALSE`.
#'
#' @details
#' Let \eqn{f}{f} be the composition of a unit in the sample, \eqn{N_f}{N_f} its total abundance
#' and \eqn{\tilde{f}}{f~} its reference composition. Then
#' \deqn{RM = \frac{1}{\left|\{f : N_f > 0\}\right|} \sum_{f \,:\, N_f > 0} \frac{\delta D_{\sigma}(f)}{\delta D_{\sigma}(\tilde{f})},}{RM = (1/|{f : N_f > 0}|) sum_{f : N_f > 0} delta D_sigma(f) / delta D_sigma(f~),}
#' with \eqn{\delta D_{\sigma}}{delta D_sigma} computed by [diversity.functional()].
#'
#' Only units present in the sample are averaged: a unit with zero total abundance is ignored, so
#' declaring extra units the sample does not have leaves the value unchanged. With
#' `include_absent = TRUE`, the average is instead taken over all \eqn{k}{k} units,
#' \deqn{RM = \frac{1}{k} \sum_{f=1}^{k} \frac{\delta D_{\sigma}(f)}{\delta D_{\sigma}(\tilde{f})},}{RM = (1/k) sum_f delta D_sigma(f) / delta D_sigma(f~),}
#' where absent units contribute 0, so a sample lacking some units gets a lower value. Each
#' reference composition must have a positive total abundance.
#'
#' @seealso [relative.multiplicity.ref_div()] when the reference diversities are already known,
#'   [relative.multiplicity.by_blocks()] for many samples at once,
#'   [diversity.functional()], [multiplicity.distance()]
#'
#' @export
#'
#' @examples
#' # Maximally different subunits with equal abundances: the diversity of a unit is its
#' # number of subunits. The three units have 4, 8 and 16 subunits in their reference.
#' max_diss <- function(n) {
#'   m <- matrix(1, n, n)
#'   diag(m) <- 0
#'   m
#' }
#' sizes <- c(4, 8, 16)
#' ab_ref <- lapply(sizes, function(n) rep(1, n))
#' diss_ref <- lapply(sizes, max_diss)
#'
#' # Half of every unit is present: RM = 0.5
#' ab <- lapply(sizes, function(n) rep(c(1, 0), n / 2))
#' relative.multiplicity(ab, diss_ref, ab_ref, diss_ref)
#'
#' # Only the last unit is present, but complete. Absent units are ignored: RM = 1
#' ab <- list(rep(0, 4), rep(0, 8), rep(1, 16))
#' relative.multiplicity(ab, diss_ref, ab_ref, diss_ref)
#'
#' # Counting the absent units as zero: RM = 1/3
#' relative.multiplicity(ab, diss_ref, ab_ref, diss_ref, include_absent = TRUE)
#'
relative.multiplicity <- function(
  ab,
  diss,
  ab_ref,
  diss_ref,
  sigma = 1,
  cap_at_one = FALSE,
  include_absent = FALSE
) {
  units <- .rm_check_units(ab, diss, "")
  refs <- .rm_check_units(ab_ref, diss_ref, "_ref")
  k <- length(units$ab)

  if (length(refs$ab) != k) {
    stop(paste0(
      "Actual and reference inputs must have the same number of units. Actual: ",
      k, ". Reference: ", length(refs$ab)
    ))
  }

  sigma <- .rm_check_sigma(sigma, k)
  .rm_check_flag(cap_at_one, "cap_at_one")
  .rm_check_flag(include_absent, "include_absent")

  ref_totals <- vapply(refs$ab, sum, numeric(1))
  if (any(ref_totals == 0)) {
    stop(paste0(
      "Reference total abundance cannot be zero. Units: ",
      paste(which(ref_totals == 0), collapse = ", ")
    ))
  }

  ref_div <- vapply(seq_len(k), function(i) {
    .rm_unit_diversity(refs$ab[[i]], refs$diss[[i]], sigma[i])
  }, numeric(1))

  .rm_ref_div_core(units$ab, units$diss, ref_div, sigma, cap_at_one, include_absent)
}


#' Relative multiplicity from reference diversities
#'
#' Computes relative multiplicity \eqn{RM}{RM} for a single sample when the reference diversity
#' \eqn{\delta D_{\sigma}(\tilde{f})}{delta D_sigma(f~)} of every unit is already known. This is the
#' same index as [relative.multiplicity()], without recomputing the reference diversities.
#'
#' @inheritParams relative.multiplicity
#' @param ref_div Numeric vector of strictly positive reference diversities, one per unit and in
#'   the same order as `ab`. They should be computed with the same `sigma` as the one given here,
#'   e.g. with [diversity.functional()].
#'
#' @return Numeric scalar, the relative multiplicity \eqn{RM}{RM}. `NA` if no unit has abundance
#'   in the sample and `include_absent = FALSE`.
#'
#' @details
#' `sigma` is only used for the diversities of the sample. As in [relative.multiplicity()], a unit
#' with zero total abundance in the sample is ignored, unless `include_absent = TRUE`.
#'
#' @seealso [relative.multiplicity()], [relative.multiplicity.by_blocks.ref_div()],
#'   [diversity.functional()]
#'
#' @export
#'
#' @examples
#' # Two units whose reference diversities are known to be 4 and 10
#' ab <- list(c(5, 5), c(1, 2, 3))
#' diss <- list(
#'   matrix(c(0, 1, 1, 0), 2, 2),
#'   matrix(c(0, 0.5, 1, 0.5, 0, 0.2, 1, 0.2, 0), 3, 3)
#' )
#' relative.multiplicity.ref_div(ab, diss, ref_div = c(4, 10))
#'
relative.multiplicity.ref_div <- function(
  ab,
  diss,
  ref_div,
  sigma = 1,
  cap_at_one = FALSE,
  include_absent = FALSE
) {
  units <- .rm_check_units(ab, diss, "")
  k <- length(units$ab)

  ref_div <- .rm_check_ref_div(ref_div, k)
  sigma <- .rm_check_sigma(sigma, k)
  .rm_check_flag(cap_at_one, "cap_at_one")
  .rm_check_flag(include_absent, "include_absent")

  .rm_ref_div_core(units$ab, units$diss, ref_div, sigma, cap_at_one, include_absent)
}

# Internal helpers
# ----------------------------

# Diversity of a single unit. Zero if the unit is absent. Zero-abundance
# subunits are dropped before computing, so only present subunits are used.
.rm_unit_diversity <- function(ab, diss, sig) {
  keep <- ab > 0
  n_keep <- sum(keep)

  if (n_keep == 0) {
    return(0)
  }
  if (n_keep == 1) {
    return(1)
  }

  if (n_keep < length(ab)) {
    if (inherits(diss, "dist")) {
      diss <- as.matrix(diss)
    }
    diss <- diss[keep, keep, drop = FALSE]
    ab <- ab[keep]
  }

  diversity.functional(ab, diss, sig)
}


# Averages the per-unit ratios for one sample (inputs already validated). Only
# the units present in the sample are averaged (NA if none is), unless
# include_absent, where absent units contribute 0 and count in the average
.rm_ref_div_core <- function(ab, diss, ref_div, sigma, cap_at_one, include_absent) {
  present <- which(vapply(ab, sum, numeric(1)) > 0)
  n_avg <- if (include_absent) length(ab) else length(present)
  if (length(present) == 0) {
    return(if (include_absent) 0 else NA_real_)
  }

  actual <- vapply(present, function(i) {
    .rm_unit_diversity(ab[[i]], diss[[i]], sigma[i])
  }, numeric(1))

  ratios <- actual / ref_div[present]
  if (cap_at_one) {
    ratios <- pmin(ratios, 1)
  }

  sum(ratios) / n_avg
}

# Wraps a single unit into a list and validates abundances and matrices
.rm_check_units <- function(ab, diss, suffix) {
  ab_name <- paste0("`ab", suffix, "`")
  diss_name <- paste0("`diss", suffix, "`")

  if (!is.list(ab)) {
    ab <- list(ab)
  }
  if (!is.list(diss) || is.data.frame(diss)) {
    diss <- list(diss)
  }

  if (length(ab) == 0) {
    stop(paste0(ab_name, " must contain at least one unit"))
  }
  if (length(ab) != length(diss)) {
    stop(paste0(
      ab_name, " and ", diss_name, " must have the same number of units. ",
      ab_name, ": ", length(ab), ". ", diss_name, ": ", length(diss)
    ))
  }

  for (i in seq_along(ab)) {
    a <- ab[[i]]
    d <- diss[[i]]

    if (!is.numeric(a)) {
      stop(paste0("Abundance vector of unit ", i, " in ", ab_name, " must be numeric"))
    }
    if (any(is.na(a))) {
      stop(paste0("Abundance vector of unit ", i, " in ", ab_name, " contains NA values"))
    }
    if (any(a < 0)) {
      stop(paste0("Abundances of unit ", i, " in ", ab_name, " must be non-negative"))
    }

    if (is.data.frame(d)) {
      d <- as.matrix(d)
      diss[[i]] <- d
    }

    if (inherits(d, "dist")) {
      n <- attr(d, "Size")
    } else {
      dims <- dim(d)
      if (length(dims) != 2 || dims[1] != dims[2]) {
        stop(paste0(
          "Distance matrix of unit ", i, " in ", diss_name, " must be square. Matrix: ",
          paste(dims, collapse = "x"), "."
        ))
      }
      n <- dims[1]
    }

    if (!is.numeric(d)) {
      stop(paste0("Distance matrix of unit ", i, " in ", diss_name, " must be numeric"))
    }
    if (n != length(a)) {
      stop(paste0(
        "Abundance vector and matrix of unit ", i, " must have compatible sizes. Matrix: ",
        n, "x", n, ". Vector: ", length(a)
      ))
    }
    if (any(is.na(d))) {
      stop(paste0("Distance matrix of unit ", i, " in ", diss_name, " contains NA values"))
    }
    if (any(d < 0)) {
      stop(paste0("Distances of unit ", i, " in ", diss_name, " must be non-negative"))
    }
  }

  list(ab = ab, diss = diss)
}


# Returns one strictly positive value per unit, matched by name when possible
.rm_per_unit_values <- function(x, k, unit_ids, name) {
  if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) || any(x <= 0)) {
    stop(paste0("`", name, "` must contain strictly positive finite numeric values"))
  }

  if (!is.null(unit_ids) && !is.null(names(x))) {
    missing <- setdiff(unit_ids, names(x))
    if (length(missing) > 0) {
      stop(paste0(
        "`", name, "` is named but missing values for units: ",
        paste(missing, collapse = ", ")
      ))
    }
    return(unname(x[unit_ids]))
  }

  if (length(x) != k) {
    stop(paste0(
      "`", name, "` must have one value per unit. Units: ", k, ". Values: ", length(x)
    ))
  }

  unname(x)
}


.rm_check_sigma <- function(sigma, k, unit_ids = NULL) {
  if (is.numeric(sigma) && length(sigma) == 1) {
    if (!is.finite(sigma) || sigma <= 0) {
      stop("`sigma` must contain strictly positive finite numeric values")
    }
    return(rep(sigma, k))
  }
  .rm_per_unit_values(sigma, k, unit_ids, "sigma")
}


.rm_check_ref_div <- function(ref_div, k, unit_ids = NULL) {
  .rm_per_unit_values(ref_div, k, unit_ids, "ref_div")
}


.rm_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(paste0("`", name, "` must be a single TRUE or FALSE value"))
  }
}
