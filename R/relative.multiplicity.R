#' Relative multiplicity
#'
#' Computes relative multiplicity \eqn{RM}{RM} of every sample: the average, across the units
#' present in the sample, of each unit's distance-based functional diversity
#' \eqn{\delta D_{\sigma}}{delta D_sigma} relative to the diversity of a reference composition of
#' that same unit. Because every unit is compared only against itself, units that are not
#' comparable with each other (e.g. different enzyme families) can be summarized together, and no
#' unit dominates the index because of its size or abundance. The reference composition of every
#' unit is built from the whole study, as set by `assume_max_reference_distance` and
#' `assume_homogeneous_abundance`, and its diversity is computed only once.
#'
#' @inheritParams raoQuadratic
#' @param diss Data frame with three columns (taken by position): the identifiers of two subunits
#'   and the distance between them, or the path to a file with them (see [distance-files]).
#'   Distances between every pair of subunits of the same unit must be listed, in either
#'   orientation; rows for subunits of different units are ignored, also when their distance is
#'   invalid (NA, infinite or negative). Pairs with a subunit that has zero abundance in every sample may be
#'   missing, with a warning, except with `assume_homogeneous_abundance = TRUE` and
#'   `assume_max_reference_distance = FALSE`, where every subunit counts in the reference.
#' @param clust Vector or factor of cluster memberships: the unit of every subunit. If named, names
#'   are subunit identifiers; every column of `abund` must be listed, and listed subunits that are
#'   not in `abund` have zero abundance in every sample. Otherwise it must follow the columns of
#'   `abund`. A data frame is also accepted, with two columns (subunit identifiers, units) or one
#'   column of units (named by its row names, if set).
#' @param sigma Numeric cutoff \eqn{\sigma}{sigma} (default `1`) at which two subunits are
#'   considered maximally different. Either a single value or a vector of strictly positive values
#'   with one value per unit. If named, values are matched to units by name (numeric names are
#'   matched as numbers, so `"1e+05"` matches the unit `100000`); otherwise they follow the order
#'   of `unique(clust)`.
#' @param cap_at_one Logical (default `FALSE`). If `TRUE`, each unit's ratio is capped at 1, for
#'   the cases where the unit is more diverse in the sample than in its reference.
#' @param include_absent Logical (default `FALSE`). If `FALSE`, each sample averages only the
#'   units present in it. If `TRUE`, every unit counts in the average, and units absent from the
#'   sample contribute 0.
#' @param assume_max_reference_distance Logical (default `FALSE`). If `FALSE`, the reference
#'   distances of a unit are the distances among all its subunits. If `TRUE`, every pair of
#'   subunits in the reference is at distance \eqn{\sigma}{sigma}.
#' @param assume_homogeneous_abundance Logical (default `FALSE`). If `FALSE`, the reference
#'   abundances of a unit are those of its subunits pooled across all samples. If `TRUE`, all
#'   subunits of the unit have the same abundance in the reference.
#' @param check_large_distance_file Logical (default `FALSE`). If `TRUE` and `diss` is a file, the
#'   file is fully checked with [check_distance_file()] before computing the index: every pair of
#'   subunits of the same unit must be listed exactly once, with a valid distance. This reads the
#'   file once more and uses temporary disk space. If `FALSE`, only the number of rows is checked
#'   and a warning is given. Ignored when `diss` is in memory, which is always fully checked. See
#'   [distance-files].
#'
#' @return Named numeric vector with the relative multiplicity \eqn{RM}{RM} of every sample.
#'   Samples with zero total abundance get `NA`.
#'
#' @details
#' Let \eqn{f}{f} be the composition of a unit in the sample, \eqn{N_f}{N_f} its total abundance
#' and \eqn{\tilde{f}}{f~} its reference composition. Then
#' \deqn{RM = \frac{1}{\left|\{f : N_f > 0\}\right|} \sum_{f \,:\, N_f > 0} \frac{\delta D_{\sigma}(f)}{\delta D_{\sigma}(\tilde{f})},}{RM = (1/|{f : N_f > 0}|) sum_{f : N_f > 0} delta D_sigma(f) / delta D_sigma(f~),}
#' with \eqn{\delta D_{\sigma}}{delta D_sigma} computed as in [diversity.functional()].
#'
#' Only units present in the sample are averaged: a unit with zero total abundance is ignored. With
#' `include_absent = TRUE`, the average is instead taken over all \eqn{k}{k} units,
#' \deqn{RM = \frac{1}{k} \sum_{f=1}^{k} \frac{\delta D_{\sigma}(f)}{\delta D_{\sigma}(\tilde{f})},}{RM = (1/k) sum_f delta D_sigma(f) / delta D_sigma(f~),}
#' where absent units contribute 0; in that case, with the pooled reference abundances
#' (`assume_homogeneous_abundance = FALSE`), a unit with zero abundance across the whole study has
#' no reference and is dropped with a warning.
#'
#' With `assume_max_reference_distance = TRUE`, the reference diversity of a unit is the inverse
#' Simpson index of its reference abundances; if `assume_homogeneous_abundance` is also `TRUE`, it
#' is simply the number of subunits of the unit.
#'
#' Every unit's diversity depends on the distances only through the sum over its pairs of
#' subunits \eqn{\sum_{i<j} a_i a_j (\min(d_{ij}, \sigma) - \sigma)}{sum_{i<j} a_i a_j
#' (min(d_ij, sigma) - sigma)}, accumulated while reading the distance table. The samples and the
#' reference are accumulated in the same read.
#'
#' @seealso [relative.multiplicity.ref_div()] when the reference diversities are already known,
#'   [diversity.functional()], [multiplicity.distance()]
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
#' relative.multiplicity(abund, diss, clust)
#' relative.multiplicity(
#'   abund, diss, clust,
#'   assume_max_reference_distance = TRUE,
#'   assume_homogeneous_abundance = TRUE
#' )
#'
relative.multiplicity <- function(
  abund,
  diss,
  clust,
  sigma = 1,
  cap_at_one = FALSE,
  include_absent = FALSE,
  assume_max_reference_distance = FALSE,
  assume_homogeneous_abundance = FALSE,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_flag(cap_at_one, "cap_at_one")
  .check_flag(include_absent, "include_absent")
  .check_flag(assume_max_reference_distance, "assume_max_reference_distance")
  .check_flag(assume_homogeneous_abundance, "assume_homogeneous_abundance")

  # The homogeneous reference with the listed distances uses every subunit of
  # the unit, also those with zero abundance in every sample
  need_all <- assume_homogeneous_abundance && !assume_max_reference_distance
  units <- .rm_prepare(
    abund, diss, clust, sigma, chunk_size, check_large_distance_file, header, need_all
  )
  n_units <- length(units$cl$unit_ids)

  # Reference abundances of every subunit
  if (assume_homogeneous_abundance) {
    ref <- rep(1, length(units$cl$unit_of))
  } else {
    ref <- as.numeric(Matrix::colSums(units$A))
  }
  ref <- matrix(ref, nrow = 1)

  # Units that never occur have no reference. They are absent from every sample,
  # so dropping them only changes the result when absent units are included
  G <- .unit_indicator(units$cl$unit_of, n_units)
  keep <- as.numeric(as.matrix(ref %*% G)) > 0
  if (!any(keep)) {
    # Every sample is empty: the distances are not read
    units$src$close()
    return(.by_sample(rep(NA_real_, nrow(units$ab$A)), units$ab))
  }
  if (include_absent && !all(keep)) {
    warning(paste0(
      "Units with zero abundance across all samples have no reference and are dropped: ",
      paste(units$cl$unit_ids[!keep], collapse = ", ")
    ))
  }

  # Reference diversities, computed once per unit
  if (assume_max_reference_distance) {
    # All pairs at distance sigma: delta D_sigma is the inverse Simpson index
    N <- as.numeric(as.matrix(ref %*% G))
    S2 <- as.numeric(as.matrix(ref^2 %*% G))
    ref_div <- N^2 / S2
    div <- .rm_unit_diversities(units$A, units$src, units$cl, units$sigma)
  } else {
    # The reference is one more row, accumulated in the same read as the samples
    div <- .rm_unit_diversities(rbind(units$A, ref), units$src, units$cl, units$sigma)
    ref_div <- as.numeric(div[nrow(div), ])
    div <- div[-nrow(div), , drop = FALSE]
  }

  .rm_average_ratios(div[, keep, drop = FALSE], ref_div[keep], units$ab, cap_at_one, include_absent)
}


#' Relative multiplicity from reference diversities
#'
#' Computes relative multiplicity \eqn{RM}{RM} of every sample, as [relative.multiplicity()]
#' does, when the reference diversity \eqn{\delta D_{\sigma}(\tilde{f})}{delta D_sigma(f~)} of
#' every unit is already known (e.g. from an external database).
#'
#' @inheritParams relative.multiplicity
#' @param ref_div Numeric vector of strictly positive reference diversities, one per unit. If named,
#'   values are matched to units by name (extra names are ignored, and numeric names are matched as
#'   numbers, as for `sigma`); otherwise they follow the order of `unique(clust)`. They should be computed with the same `sigma` as the one given here, e.g.
#'   with [diversity.functional()].
#'
#' @return Named numeric vector with the relative multiplicity \eqn{RM}{RM} of every sample.
#'   Samples with zero total abundance get `NA`.
#'
#' @details
#' `sigma` is only used for the diversities of the samples. As in [relative.multiplicity()], only
#' the units present in a sample count in its average, unless `include_absent = TRUE`, where every
#' unit in `clust` counts, including units with zero abundance in every sample.
#'
#' @seealso [relative.multiplicity()], [diversity.functional()]
#'
#' @export
#'
#' @examples
#' clust <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B")
#' diss <- data.frame(
#'   ID1 = c("a1", "a1", "a2", "b1"),
#'   ID2 = c("a2", "a3", "a3", "b2"),
#'   Distance = c(0.3, 0.8, 0.6, 0.4)
#' )
#' abund <- matrix(
#'   c(5, 3, 0, 2, 0,
#'     0, 1, 4, 2, 6),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(c("S1", "S2"), names(clust))
#' )
#'
#' # Reference diversities from an external database
#' relative.multiplicity.ref_div(abund, diss, clust, ref_div = c(A = 2.5, B = 1.8))
#'
relative.multiplicity.ref_div <- function(
  abund,
  diss,
  clust,
  ref_div,
  sigma = 1,
  cap_at_one = FALSE,
  include_absent = FALSE,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_flag(cap_at_one, "cap_at_one")
  .check_flag(include_absent, "include_absent")

  units <- .rm_prepare(abund, diss, clust, sigma, chunk_size, check_large_distance_file, header)
  ref_div <- .rm_check_ref_div(ref_div, units$cl$unit_ids, units$input_order)

  div <- .rm_unit_diversities(units$A, units$src, units$cl, units$sigma)
  .rm_average_ratios(div, ref_div, units$ab, cap_at_one, include_absent)
}


# Internal helpers
# ----------------------------

# Validates the inputs and opens the distances. Subunits listed in a named
# `clust` but not in `abund` are added with zero abundance in every sample.
# `sigma` is NULL for indices without a cutoff. The pairs of subunits with zero
# abundance in every sample are only required with `need_all`
.rm_prepare <- function(abund, diss, clust, sigma, chunk_size, check, header = TRUE,
                        need_all = FALSE) {
  ab <- .parse_abund(abund)
  clust <- .as_clust_vector(clust)
  A <- ab$A
  subunits <- ab$subunits

  if (!is.null(names(clust))) {
    extra <- setdiff(.as_ids(names(clust)), subunits)
    if (length(extra) > 0) {
      zeros <- matrix(0, nrow = nrow(A), ncol = length(extra))
      if (inherits(A, "Matrix")) {
        zeros <- Matrix::Matrix(zeros, sparse = TRUE)
      }
      A <- cbind(A, zeros)
      subunits <- c(subunits, extra)
      colnames(A) <- subunits
    }
  }

  # Unnamed per-unit values follow unique(clust) as given, before .parse_clust
  # reorders a named `clust` into the columns of `abund`
  input_order <- unique(.as_ids(clust))

  cl <- .parse_clust(clust, subunits)
  if (!is.null(sigma)) {
    sigma <- .rm_check_sigma(sigma, cl$unit_ids, input_order)
  }
  need <- if (need_all) NULL else as.numeric(Matrix::colSums(A)) > 0
  src <- .diss_source(diss, subunits, chunk_size, check, "within", cl, header = header, need = need)

  list(ab = ab, A = A, cl = cl, sigma = sigma, src = src, input_order = input_order)
}


# Diversity delta D_sigma of every unit (columns) in every row of A, zero for
# absent units, reading the pairs from `src`. With
# Q = (sigma (N^2 - sum a^2) + 2 acc) / N^2, where
#   acc = sum over pairs (i, j) of the unit of a_i a_j (min(d_ij, sigma) - sigma),
# delta D_sigma = 1 / (1 - Q / sigma) = sigma N^2 / (sigma sum a^2 - 2 acc)
.rm_unit_diversities <- function(A, src, cl, sigma) {
  acc <- .rm_unit_diversities_acc(A, cl, sigma)
  .consume_pairs(src, list(acc), colnames(A), cl)[[1]]
}


.rm_unit_diversities_acc <- function(A, cl, sigma) {
  n_units <- length(cl$unit_ids)
  acc <- matrix(0, nrow = nrow(A), ncol = n_units)

  update <- function(pairs) {
    within <- cl$unit_of[pairs$i] == cl$unit_of[pairs$j]
    i <- pairs$i[within]
    j <- pairs$j[within]
    g <- cl$unit_of[i]
    w <- pmin(pairs$d[within], sigma[g]) - sigma[g]
    acc <<- acc + .pair_sums(A, i, j, w, group = g, n_groups = n_units)
  }

  finish <- function() {
    G <- .unit_indicator(cl$unit_of, n_units)
    N <- as.matrix(A %*% G)
    S2 <- as.matrix(A^2 %*% G)
    sig <- matrix(sigma, nrow = nrow(N), ncol = n_units, byrow = TRUE)
    div <- sig * N^2 / (sig * S2 - 2 * acc)
    div[N == 0] <- 0
    div
  }

  list(scope = "within", update = update, finish = finish)
}


# Averages the per-unit ratios of every sample. Only the units present in the
# sample are averaged, unless include_absent, where absent units contribute 0
# and count in the average
.rm_average_ratios <- function(div, ref_div, ab, cap_at_one, include_absent) {
  ratios <- sweep(div, 2, ref_div, "/")
  if (cap_at_one) {
    ratios <- pmin(ratios, 1)
  }

  .average_units(ratios, div > 0, ab, include_absent)
}


# Averages the per-unit values `vals` of every sample over the units `present`
# in it, unless include_absent, where absent units contribute 0 and count in
# the average. Shared by relative multiplicity and average redundancy
.average_units <- function(vals, present, ab, include_absent) {
  n_present <- rowSums(present)
  n_avg <- if (include_absent) ncol(vals) else n_present

  res <- rowSums(vals * present) / n_avg
  res[n_present == 0] <- if (include_absent) 0 else NA_real_
  .by_sample(res, ab)
}


# Returns one strictly positive value per unit, matched by name when possible.
# Unnamed values follow `input_order`, the order of unique(clust)
.rm_per_unit_values <- function(x, unit_ids, input_order, name) {
  k <- length(unit_ids)
  if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) || any(x <= 0)) {
    stop(paste0("`", name, "` must contain strictly positive finite numeric values"))
  }

  if (!is.null(names(x))) {
    pos <- .match_names(unit_ids, names(x))
    if (anyNA(pos)) {
      stop(paste0(
        "`", name, "` is named but missing values for units: ",
        paste(unit_ids[is.na(pos)], collapse = ", ")
      ))
    }
    return(unname(x[pos]))
  }

  if (length(x) != k) {
    stop(paste0(
      "`", name, "` must have one value per unit. Units: ", k, ". Values: ", length(x)
    ))
  }

  unname(x[match(unit_ids, input_order)])
}


.rm_check_sigma <- function(sigma, unit_ids, input_order) {
  if (is.numeric(sigma) && length(sigma) == 1) {
    if (!is.finite(sigma) || sigma <= 0) {
      stop("`sigma` must contain strictly positive finite numeric values")
    }
    return(rep(sigma, length(unit_ids)))
  }
  .rm_per_unit_values(sigma, unit_ids, input_order, "sigma")
}


.rm_check_ref_div <- function(ref_div, unit_ids, input_order) {
  .rm_per_unit_values(ref_div, unit_ids, input_order, "ref_div")
}