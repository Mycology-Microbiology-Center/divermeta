#' Distance-based functional diversity (Chiu & Chao 2014)
#'
#' Computes distance-based functional diversity \eqn{\delta D_{\sigma}}{delta D_sigma} following Chiu & Chao (2014)
#' for every sample. Pairwise distances are capped at the cutoff \eqn{\sigma}{sigma}.
#'
#' @inheritParams raoQuadratic
#' @param sig Numeric cutoff \eqn{\sigma}{sigma} at which two units are considered different (default `1`).
#'
#' @return Named numeric vector with the distance-based functional diversity
#'   \eqn{\delta D_{\sigma}}{delta D_sigma} of every sample. Samples with zero total abundance get `NA`.
#' @references
#' \itemize{
#' \item Chiu CH, Chao A (2014) Distance-based functional diversity measures and their decomposition: A framework based on Hill numbers. PLOS ONE 9(7). \doi{10.1371/journal.pone.0100014}. \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014}
#' }
#' @seealso [raoQuadratic()], [diversity.functional.traditional()]
#' @export
#'
#' @examples
#' abund <- matrix(c(2, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("S1", "S2"), c("a", "b", "c")))
#' diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))
#' diversity.functional(abund, diss, sig = 0.8)
#'
diversity.functional <- function(
  abund,
  diss,
  sig = 1,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_sig(sig)
  ab <- .parse_abund(abund)
  src <- .diss_source(
    diss, ab$subunits, chunk_size, check_large_distance_file, "all",
    header = header, need = .present_subunits(ab)
  )
  .diversity.functional_core(ab, src, sig)
}


.diversity.functional_core <- function(ab, src, sig) {
  .consume_pairs(src, list(.diversity.functional_acc(ab, sig)), ab$subunits)[[1]]
}


# Accumulates Rao's quadratic entropy with the distances capped at sig
.diversity.functional_acc <- function(ab, sig) {
  P <- .relative(ab)
  Q <- numeric(nrow(P))
  list(
    scope = "all",
    update = function(pairs) {
      pairs$d <- pmin(pairs$d, sig)
      Q <<- Q + .rao_q(P, pairs)
    },
    finish = function() .by_sample(1 / (1 - Q / sig), ab)
  )
}



#' Distance-based functional diversity (order q) (Chiu & Chao 2014)
#'
#' Computes distance-based functional diversity \eqn{^{q}FD}{FD^q} following Chiu & Chao (2014)
#' for every sample. This corresponds to D(Q) in their paper and \eqn{^{q}FD}{FD^q} in the
#' divermeta manuscript. For \eqn{q = 1}{q = 1}, the analytic limit is used.
#'
#' @inheritParams raoQuadratic
#' @param q Numeric order of the Hill number (non-negative, default `1`).
#'
#' @return Named numeric vector with the distance-based functional diversity \eqn{^{q}FD}{FD^q}
#'   of every sample. Samples with zero total abundance get `NA`.
#'
#' @details
#' Only the subunits present in a sample (nonzero abundance) contribute, also for \eqn{q = 0}{q = 0}.
#'
#' When Rao's quadratic entropy \eqn{Q}{Q} of a sample is 0 (a single present subunit, or all
#' distances between its present subunits 0), the formula is undefined (0/0) and the value is set
#' to 1, as [diversity.functional()] gives for the same sample. Values of \eqn{q}{q} within
#' `1e-8` of 1 use the analytic limit at \eqn{q = 1}{q = 1}.
#'
#' @references
#' \itemize{
#' \item Chiu CH, Chao A (2014) Distance-based functional diversity measures and their decomposition: A framework based on Hill numbers. PLOS ONE 9(7). \doi{10.1371/journal.pone.0100014}. \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014}
#' }
#' @seealso [diversity.functional()], [raoQuadratic()]
#' @export
#'
#' @examples
#' abund <- matrix(c(2, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("S1", "S2"), c("a", "b", "c")))
#' diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))
#' diversity.functional.traditional(abund, diss, q = 2)
#'
diversity.functional.traditional <- function(
  abund,
  diss,
  q = 1,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_q(q)
  ab <- .parse_abund(abund)
  src <- .diss_source(
    diss, ab$subunits, chunk_size, check_large_distance_file, "all",
    header = header, need = .present_subunits(ab)
  )
  .diversity.functional.traditional_core(ab, src, q)
}


.diversity.functional.traditional_core <- function(ab, src, q) {
  .consume_pairs(src, list(.diversity.functional.traditional_acc(ab, q)), ab$subunits)[[1]]
}


# Accumulates Rao's quadratic entropy and the sum of order q in the same read
.diversity.functional.traditional_acc <- function(ab, q) {
  P <- .relative(ab)
  Q <- numeric(nrow(P))
  acc <- numeric(nrow(P))
  q_one <- abs(q - 1) < .q_one_tol

  if (q_one) {
    # Analytic limit q -> 1: exp(-sum_ij d_ij p_i p_j ln(p_i) / Q) (square root
    # already included), with w = p ln(p)
    W <- P
    if (inherits(W, "Matrix")) {
      nz <- W@x > 0
      W@x[nz] <- W@x[nz] * log(W@x[nz])
    } else {
      nz <- W > 0
      W[nz] <- W[nz] * log(W[nz])
    }
  } else {
    # Only present subunits contribute (0^0 would be 1)
    Pq <- if (q == 0) (P > 0) * 1 else P^q
  }

  update <- function(pairs) {
    Q <<- Q + .rao_q(P, pairs)
    if (q_one) {
      # sum over ordered pairs of w_i d_ij p_j
      acc <<- acc + .pair_sums(W, pairs$i, pairs$j, pairs$d, Y = P)[, 1] +
        .pair_sums(W, pairs$j, pairs$i, pairs$d, Y = P)[, 1]
    } else {
      acc <<- acc + 2 * .pair_sums(Pq, pairs$i, pairs$j, pairs$d)[, 1]
    }
  }

  finish <- function() {
    if (q_one) {
      vals <- exp(-acc / Q)
    } else {
      vals <- sqrt((acc / Q)^(1 / (1 - q)))
    }
    # Q = 0 (one present subunit, or all distances 0): a single functional type.
    # Q is exactly 0 then, and a tolerance would break the invariance to the
    # scale of the distances
    vals[Q <= 0] <- 1
    .by_sample(vals, ab)
  }

  list(scope = "all", update = update, finish = finish)
}



#' Functional redundancy (Re) (Ricotta & Pavoine 2025)
#'
#' Computes functional redundancy \eqn{Re}{Re}, a measure of the degree to which
#' distinct elements are functionally similar given their abundances and
#' pairwise dissimilarities, for every sample. This implementation follows the Simpson–Rao
#' family and corresponds to the `q = 2` case.
#'
#' @inheritParams raoQuadratic
#' @param diss Data frame with three columns (taken by position): the identifiers of two subunits
#'   and the dissimilarity between them, scaled to the range \[0, 1\] (distances above 1 are
#'   capped at 1), or the path to a file with them (see [distance-files]). Every pair of subunits
#'   present in at least one sample must be listed, in either orientation; pairs with a subunit
#'   that has zero abundance in every sample may be missing, with a warning.
#'
#' @param normalize Logical (default `FALSE`). If `TRUE`, the redundancy is divided by Simpson's
#'   diversity \eqn{S}{S}, giving \eqn{(S - Q) / S = 1 - Q / S}{(S - Q) / S = 1 - Q / S}.
#'
#' @return Named numeric vector with the functional redundancy `Re` of every sample. Samples
#'   with zero total abundance get `NA`.
#'
#' @details
#' With relative abundances \eqn{p_i}{p_i}, Simpson's diversity is
#' \eqn{S = 1 - \sum_i p_i^2}{S = 1 - sum_i p_i^2} and Rao's quadratic entropy is
#' \eqn{Q = \sum_{i \neq j} p_i p_j d_{ij}}{Q = sum_{i != j} p_i p_j d_ij}. Then
#' \deqn{Re = S - Q,}{Re = S - Q,}
#' and, with `normalize = TRUE`,
#' \deqn{Re = \frac{S - Q}{S} = 1 - \frac{Q}{S}.}{Re = (S - Q) / S = 1 - Q / S.}
#' A sample with a single present subunit has no redundancy: its value is 0, also when normalized.
#'
#' Distances greater than 1 are capped at 1 (\eqn{d_{ij} = \min(d_{ij}, 1)}{d_ij = min(d_ij, 1)})
#' before computing \eqn{Q}{Q}, so that \eqn{Q \le S}{Q <= S} and the redundancy is never negative:
#' it lies between 0 and \eqn{S}{S}, or between 0 and 1 when normalized.
#'
#' @references
#' \itemize{
#' \item Ricotta C, Pavoine S (2025) What do functional diversity, redundancy, rarity, and originality actually measure? A theoretical guide for ecologists and conservationists. Ecological Complexity 61. \doi{10.1016/j.ecocom.2025.101116}. \url{https://www.sciencedirect.com/science/article/pii/S1476945X25000017}
#' \item Rao CR (1982) Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology 21. \doi{10.1016/0040-5809(82)90004-1}.
#' }
#' @seealso [raoQuadratic()], [average.redundancy()]
#' @export
#'
#' @examples
#' abund <- matrix(c(2, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("S1", "S2"), c("a", "b", "c")))
#' diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))
#' redundancy(abund, diss)
#' redundancy(abund, diss, normalize = TRUE)
#'
redundancy <- function(
  abund,
  diss,
  normalize = FALSE,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_flag(normalize, "normalize")
  ab <- .parse_abund(abund)
  src <- .diss_source(
    diss, ab$subunits, chunk_size, check_large_distance_file, "all",
    header = header, need = .present_subunits(ab)
  )
  .redundancy_core(ab, src, normalize)
}


.redundancy_core <- function(ab, src, normalize = FALSE) {
  .consume_pairs(src, list(.redundancy_acc(ab, normalize)), ab$subunits)[[1]]
}


# Accumulates Rao's quadratic entropy with the distances capped at 1
.redundancy_acc <- function(ab, normalize = FALSE) {
  P <- .relative(ab)
  Q <- numeric(nrow(P))
  list(
    scope = "all",
    update = function(pairs) {
      pairs$d <- pmin(pairs$d, 1)
      Q <<- Q + .rao_q(P, pairs)
    },
    finish = function() {
      # Simpson's diversity (q = 2): D = 1 - sum p_i^2
      d <- 1 - as.numeric(Matrix::rowSums(P^2))

      # Functional redundancy
      re <- if (normalize) (d - Q) / d else d - Q

      # A single present subunit has no redundancy. Counted, since d is only
      # approximately 0 in floating point (e.g. 49 * (1 / 49) != 1)
      re[as.numeric(Matrix::rowSums(ab$A > 0)) <= 1] <- 0
      .by_sample(re, ab)
    }
  )
}


#' Metagenomic Alpha-Diversity Index (MAD) (Finn 2024)
#'
#' Computes the Metagenomic Alpha-Diversity Index (MAD) of every sample, a metric that measures
#' the average dissimilarity of elements (e.g., protein-encoding genes) within clusters
#' relative to cluster representatives. Unlike multiplicity, MAD does not account for
#' element abundances and decreases as the number of elements per cluster increases.
#'
#' @inheritParams raoQuadratic
#' @param diss Data frame with three columns (taken by position): the identifiers of two subunits
#'   and the dissimilarity between them, scaled to the range \[0, 1\], or the path to a file with
#'   them (see [distance-files]). The distance between every representative and every other
#'   subunit of its cluster present in a sample must be listed. Only these pairs are kept while
#'   reading, so missing pairs and pairs listed again with a different distance are always
#'   detected, also in files (there is no need for `check_large_distance_file`).
#' @param clust Vector or factor of cluster memberships for each subunit. If named, names are
#'   subunit identifiers; otherwise it follows the columns of `abund`. A data frame is also accepted, with two columns (subunit identifiers, units) or one
#'   column of units (named by its row names, if set).
#' @param representatives Optional named vector mapping cluster names to the identifier of the
#'   representative subunit of each cluster. Numeric identifiers are matched as numbers, so a name
#'   `"1e+05"` matches the cluster `100000`. If `NULL` (default), the first subunit (in the column
#'   order of `abund`) of each cluster present in the sample is used as the representative.
#'
#' @return Named numeric vector with the Metagenomic Alpha-Diversity Index (MAD) of every sample.
#'   Higher values indicate greater average dissimilarity within clusters. Samples with zero total
#'   abundance get `NA`.
#'
#' @details
#' Abundances are only used to decide which subunits are present (nonzero abundance) in every
#' sample: the index of a sample is computed over its present subunits and clusters. A given
#' representative is used as the reference of its cluster even when it is absent from a sample.
#'
#' Note: Unlike multiplicity indices, MAD does not incorporate element abundances and
#' decreases as cluster size increases, which may not reflect biological complexity.
#'
#' Note on duplicated pairs: an in-memory `diss` is validated as a whole, so a pair listed again
#' with a different distance is an error even when MAD does not use that pair. When `diss` is a
#' file, only the pairs MAD uses (representative and subunit of its cluster) are checked for
#' conflicting duplicates, and the others are ignored. In both cases every distance must be a
#' valid (non-negative, not NA) number.
#'
#' @references
#' \itemize{
#' \item Finn DR (2024) A metagenomic alpha-diversity index for microbial functional
#'   biodiversity. FEMS Microbiology Ecology 100(3):fiae019.
#'   \doi{10.1093/femsec/fiae019}
#' }
#'
#' @seealso [multiplicity.distance()] for abundance-weighted distance-based multiplicity,
#'   [multiplicity.inventory()] for inventory multiplicity
#'
#' @export
#'
#' @examples
#' # Example: Compute MAD for gene clusters
#' clust <- c(g1 = 1, g2 = 1, g3 = 2, g4 = 2, g5 = 2, g6 = 3)
#' diss <- data.frame(
#'   ID1 = c("g1", "g3", "g3", "g4"),
#'   ID2 = c("g2", "g4", "g5", "g5"),
#'   Distance = c(0.2, 0.1, 0.15, 0.12)
#' )
#' abund <- matrix(
#'   c(1, 1, 1, 1, 1, 1,
#'     1, 0, 0, 1, 1, 1),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(c("S1", "S2"), names(clust))
#' )
#'
#' # Use default (first present element as representative)
#' metagenomic.alpha.index(abund, diss, clust)
#'
#' # Specify representatives explicitly
#' reps <- c("1" = "g1", "2" = "g3", "3" = "g6")
#' metagenomic.alpha.index(abund, diss, clust, representatives = reps)
#'
metagenomic.alpha.index <- function(
  abund,
  diss,
  clust,
  representatives = NULL,
  chunk_size = 1e6,
  header = TRUE
) {
  ab <- .parse_abund(abund)
  cl <- .parse_clust(clust, ab$subunits)
  rep_of_unit <- .mad_representatives(representatives, ab, cl)
  src <- .diss_source(diss, ab$subunits, chunk_size, scope = "none", cl = cl, header = header)
  .metagenomic.alpha.index_core(ab, src, cl, rep_of_unit)
}


.metagenomic.alpha.index_core <- function(ab, src, cl, rep_of_unit) {
  .consume_pairs(src, list(.metagenomic.alpha.index_acc(ab, cl, rep_of_unit)), ab$subunits, cl)[[1]]
}


# Representative subunit (column index) of every unit, or NULL to use the first
# present subunit of the unit in every sample
.mad_representatives <- function(representatives, ab, cl) {
  if (is.null(representatives)) {
    return(NULL)
  }
  n_units <- length(cl$unit_ids)
  if (is.null(names(representatives))) {
    stop("`representatives` must be a named vector (names are cluster identifiers)")
  }
  pos <- .match_names(cl$unit_ids, names(representatives))
  if (anyNA(pos)) {
    stop(paste0(
      "Missing representatives for clusters: ",
      paste(cl$unit_ids[is.na(pos)], collapse = ", ")
    ))
  }
  rep_ids <- .as_ids(representatives[pos])
  rep_of_unit <- .match_names(rep_ids, ab$subunits)
  if (anyNA(rep_of_unit)) {
    stop(paste0(
      "Representatives must be subunits of `abund`: ",
      paste(rep_ids[is.na(rep_of_unit)], collapse = ", ")
    ))
  }
  wrong <- cl$unit_of[rep_of_unit] != seq_len(n_units)
  if (any(wrong)) {
    stop(paste0(
      "Representatives must belong to their cluster: ",
      paste(rep_ids[wrong], collapse = ", ")
    ))
  }
  rep_of_unit
}


# Looks up, while the pairs are read, the distance between every present subunit
# and the representative of its cluster in its sample. Only these pairs are kept,
# so missing and conflicting pairs are found exactly, also in files
.metagenomic.alpha.index_acc <- function(ab, cl, rep_of_unit) {
  n_units <- length(cl$unit_ids)

  # Present subunits of every sample
  nz <- .nonzero(ab$A)
  s <- nz$i
  m <- nz$j
  u <- cl$unit_of[m]
  su <- (s - 1) * as.numeric(n_units) + u

  # Representative of every present subunit's cluster in its sample
  if (is.null(rep_of_unit)) {
    o <- order(su, m)
    first <- !duplicated(su[o])
    rep <- m[o][first][match(su, su[o][first])]
  } else {
    rep <- rep_of_unit[u]
  }

  # Needed pairs (representative, subunit), as keys of the pair i < j
  n <- as.numeric(length(ab$subunits))
  other <- rep != m
  query <- pmin(rep, m)[other] + (pmax(rep, m)[other] - 1) * n
  needed <- unique(query)
  found <- rep(NA_real_, length(needed))

  update <- function(pairs) {
    pos <- match(pairs$i + (pairs$j - 1) * n, needed)
    hit <- !is.na(pos)
    if (!any(hit)) {
      return(invisible(NULL))
    }
    pos <- pos[hit]
    d <- pairs$d[hit]
    # A pair listed again must have the same distance
    prev <- c(found[pos], d[match(pos, pos)])
    conflict <- !is.na(prev) & abs(prev - c(d, d)) > sqrt(.Machine$double.eps)
    if (any(conflict)) {
      k <- pos[(which(conflict)[1] - 1) %% length(pos) + 1]
      i <- (needed[k] - 1) %% n + 1
      j <- (needed[k] - 1) %/% n + 1
      stop(paste0(
        "`diss` lists conflicting distances for the same pair: ",
        ab$subunits[i], "-", ab$subunits[j]
      ))
    }
    found[pos] <<- d
    invisible(NULL)
  }

  finish <- function() {
    # Distance to the representative (zero for the representative itself)
    d <- numeric(length(m))
    d_other <- found[match(query, needed)]
    if (anyNA(d_other)) {
      k <- which(other)[which(is.na(d_other))[1]]
      stop(paste0(
        "Missing distances in `diss` between representatives and subunits of their cluster (e.g. ",
        ab$subunits[rep[k]], "-", ab$subunits[m[k]], ")"
      ))
    }
    d[other] <- d_other

    # Per sample: sum over present clusters of (1 + mean distance), over present subunits
    groups <- unique(su)
    g <- match(su, groups)
    mean_d <- as.vector(rowsum(d, g, reorder = FALSE)) / tabulate(g)
    sample_of_group <- s[match(groups, su)]
    n_samples <- length(ab$samples)
    total <- numeric(n_samples)
    if (length(groups) > 0) {
      by_sample <- rowsum(1 + mean_d, sample_of_group)
      total[as.integer(rownames(by_sample))] <- by_sample[, 1]
    }
    n_present <- tabulate(s, n_samples)

    .by_sample(total / n_present, ab)
  }

  list(scope = "none", update = update, finish = finish)
}
