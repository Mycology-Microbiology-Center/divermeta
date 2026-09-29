# Shared internals
# ----------------------------
# Every index receives the same three kinds of input:
#   - abund: abundance table with samples in rows and subunits in columns
#   - clust: vector assigning every subunit (column) to its unit
#   - diss:  three-column table (id_1, id_2, distance) with pairwise distances
# The helpers below validate and index them once, and accumulate sums over the
# pairs of the distance table chunk by chunk, so no distance matrix is built.

# Number of cells (samples x pairs) processed at a time by .pair_sums
.chunk_cells <- 1e7

# Orders q closer than this to 1 use the analytic limit at q = 1, since the
# general formula loses precision there
.q_one_tol <- 1e-8


# Identifiers as character, without scientific notation for whole numbers
.as_ids <- function(ids) {
  if (is.factor(ids)) {
    return(as.character(ids))
  }
  if (is.numeric(ids) && all(ids == round(ids), na.rm = TRUE)) {
    return(sprintf("%.0f", ids))
  }
  as.character(ids)
}


# Validates the abundance table. Returns the abundances (base matrix or sparse
# dgCMatrix), the sample and subunit identifiers and the total of every sample
.parse_abund <- function(abund) {
  is_sparse <- inherits(abund, "Matrix")
  if (is.data.frame(abund)) {
    abund <- as.matrix(abund)
  }
  if (!is_sparse && !is.matrix(abund)) {
    stop("`abund` must be a matrix, data.frame or sparse Matrix (samples in rows, subunits in columns)")
  }
  if (!is_sparse && !is.numeric(abund)) {
    stop("`abund` must be numeric")
  }
  if (is_sparse) {
    abund <- methods::as(methods::as(methods::as(abund, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    vals <- abund@x
  } else {
    vals <- abund
  }
  if (anyNA(vals)) {
    stop("`abund` contains NA values")
  }
  if (any(is.infinite(vals))) {
    stop("`abund` contains infinite values")
  }
  if (any(vals < 0)) {
    stop("Abundances must be non-negative")
  }

  subunits <- colnames(abund)
  if (is.null(subunits)) {
    stop("`abund` must have column names with the subunit identifiers")
  }
  if (anyDuplicated(subunits) > 0) {
    stop("Column names of `abund` (subunit identifiers) must be unique")
  }
  samples <- rownames(abund)
  if (is.null(samples)) {
    samples <- as.character(seq_len(nrow(abund)))
  }

  list(
    A = abund,
    samples = samples,
    subunits = subunits,
    N = as.numeric(Matrix::rowSums(abund))
  )
}


# Subunits present (nonzero abundance) in at least one sample. The others need
# no distances
.present_subunits <- function(ab) {
  as.numeric(Matrix::colSums(ab$A)) > 0
}


# Clustering as a vector, named by subunit when the names are known. Accepts a
# vector or factor, a one-column matrix or data frame (units, with the subunit
# identifiers as row names, if set) or a two-column data frame (subunit
# identifiers, units)
.as_clust_vector <- function(clust) {
  if (is.data.frame(clust)) {
    if (ncol(clust) == 2) {
      return(stats::setNames(clust[[2]], as.character(clust[[1]])))
    }
    if (ncol(clust) != 1) {
      stop(paste0(
        "A data frame `clust` must have one column (units, with the subunit identifiers as ",
        "row names) or two columns (subunit identifiers, units)"
      ))
    }
    # Automatic row names (1, 2, ...) are not subunit identifiers
    ids <- if (.row_names_info(clust) > 0) rownames(clust) else NULL
    return(stats::setNames(clust[[1]], ids))
  }
  if (!is.null(dim(clust))) {
    if (ncol(clust) != 1) {
      stop("A matrix `clust` must have a single column (units, with the subunit identifiers as row names)")
    }
    return(stats::setNames(as.vector(clust), rownames(clust)))
  }
  clust
}


# Validates the clustering of the subunits. A named `clust` is matched to the
# subunits by name; an unnamed one must follow the order of the subunits.
# Returns the unit (integer) of every subunit and the unit identifiers
.parse_clust <- function(clust, subunits) {
  if (is.null(clust) || length(clust) == 0) {
    stop("`clust` cannot be empty")
  }
  clust <- .as_clust_vector(clust)
  if (anyNA(clust)) {
    stop("`clust` contains NA values")
  }

  if (!is.null(names(clust))) {
    if (anyNA(names(clust)) || anyDuplicated(names(clust)) > 0) {
      stop("`clust` names (subunit identifiers) must be unique and not NA")
    }
    missing <- setdiff(subunits, names(clust))
    if (length(missing) > 0) {
      stop(paste0(
        "`clust` is named but missing some subunits present in `abund`: ",
        paste(utils::head(missing, 5), collapse = ", "),
        if (length(missing) > 5) ", ..." else ""
      ))
    }
    clust <- clust[subunits]
  } else if (length(clust) != length(subunits)) {
    stop(paste0(
      "`abund` and `clust` must have the same number of subunits. `abund`: ",
      length(subunits), ". `clust`: ", length(clust)
    ))
  }

  clust_chr <- .as_ids(clust)
  unit_ids <- unique(clust_chr)
  list(unit_of = match(clust_chr, unit_ids), unit_ids = unit_ids)
}


# Positions of `x` in the identifiers `ids`. Numeric `x` is matched numerically
# when every identifier is a distinct number, so 100000 matches "1e+05"
.match_ids <- function(x, ids) {
  if (is.numeric(x)) {
    ids_num <- suppressWarnings(as.numeric(ids))
    if (!anyNA(ids_num) && anyDuplicated(ids_num) == 0) {
      return(match(as.numeric(x), ids_num))
    }
  }
  match(.as_ids(x), ids)
}


# Positions of the identifiers `x` in the names `table`. Identifiers without an
# exact match are matched as numbers, so "1e+05" matches "100000", unless two
# names in `table` are the same number
.match_names <- function(x, table) {
  x <- .as_ids(x)
  table <- as.character(table)
  pos <- match(x, table)
  todo <- is.na(pos)
  if (any(todo)) {
    x_num <- suppressWarnings(as.numeric(x[todo]))
    t_num <- suppressWarnings(as.numeric(table))
    if (anyDuplicated(t_num, incomparables = NA) == 0) {
      pos[todo] <- match(x_num, t_num, incomparables = NA)
    }
  }
  pos
}


# Indexes rows of distances by `ids`: rows with unknown identifiers and self
# pairs are dropped, and so are the rows outside `scope` ("within" or "between"
# units of the clustering `cl`; "all" and "none" keep every row), so their
# distances are never checked. Returns the pairs as indices i < j into `ids`
# with their distances d, and the number of rows pairing two known subunits
# (`n_known`, before the scope). Shared by in-memory tables and the chunks
# read from files
.index_pairs <- function(id1, id2, d, ids, name = "diss", scope = "all", cl = NULL) {
  arg <- paste0("`", name, "`")
  if (!is.numeric(d)) {
    stop(paste0("The distance column (third) of ", arg, " must be numeric"))
  }

  i1 <- .match_ids(id1, ids)
  i2 <- .match_ids(id2, ids)
  keep <- !is.na(i1) & !is.na(i2) & i1 != i2
  n_known <- sum(keep)
  if (scope %in% c("within", "between")) {
    same_unit <- cl$unit_of[i1[keep]] == cl$unit_of[i2[keep]]
    keep[keep] <- if (scope == "within") same_unit else !same_unit
  }
  i1 <- i1[keep]
  i2 <- i2[keep]
  d <- d[keep]

  if (anyNA(d)) {
    stop(paste0(arg, " contains NA distances"))
  }
  if (any(is.infinite(d))) {
    stop(paste0(arg, " contains infinite distances"))
  }
  if (any(d < 0)) {
    stop("Distances must be non-negative")
  }

  list(i = pmin(i1, i2), j = pmax(i1, i2), d = as.numeric(d), n_known = n_known)
}


# Validates the three-column distance table and indexes it by `ids`. Rows with
# unknown identifiers, self pairs and rows outside `scope` are dropped (see
# .index_pairs), and every unordered pair is kept once (listing it again is
# allowed only with the same distance).
# Returns the pairs as indices i < j into `ids` with their distances d
.parse_diss <- function(diss, ids, name = "diss", scope = "all", cl = NULL) {
  arg <- paste0("`", name, "`")
  if (!is.data.frame(diss) && !is.matrix(diss)) {
    stop(paste0(
      arg, " must be a data frame with three columns (id_1, id_2, distance) ",
      "or the path to a file with them"
    ))
  }
  if (ncol(diss) != 3) {
    stop(paste0(arg, " must have three columns: id_1, id_2, distance"))
  }
  diss <- as.data.frame(diss, stringsAsFactors = FALSE)

  pairs <- .index_pairs(diss[[1]], diss[[2]], diss[[3]], ids, name, scope, cl)
  i <- pairs$i
  j <- pairs$j
  d <- pairs$d

  # A pair may be listed more than once (e.g. in both orientations) only with the same distance
  key <- i + (j - 1) * as.numeric(length(ids))
  dup <- duplicated(key)
  if (any(dup)) {
    first <- match(key, key)
    conflict <- which(abs(d - d[first]) > sqrt(.Machine$double.eps))
    if (length(conflict) > 0) {
      k <- conflict[1]
      stop(paste0(
        arg, " lists conflicting distances for the same pair: ",
        ids[i[k]], "-", ids[j[k]]
      ))
    }
    i <- i[!dup]
    j <- j[!dup]
    d <- d[!dup]
  }

  list(i = i, j = j, d = d)
}


# Subset of the pairs
.subset_pairs <- function(pairs, keep) {
  list(i = pairs$i[keep], j = pairs$j[keep], d = pairs$d[keep])
}


# Returns a pair of `members` (sorted indices) missing from the pairs, all of
# which are pairs of members, or NULL if every pair of members is listed
.missing_pair <- function(i, j, members) {
  m <- length(members)
  pos <- match(i, members)
  counts <- tabulate(pos, m)
  short <- which(counts < m - seq_len(m))
  if (length(short) == 0) {
    return(NULL)
  }
  k <- short[1]
  partner <- setdiff(members[(k + 1):m], j[pos == k])[1]
  c(members[k], partner)
}


# Stops unless every pair of the `n` identifiers is listed
.require_all_pairs <- function(pairs, ids, name = "diss", what = "subunits") {
  n <- length(ids)
  expected <- n * (n - 1) / 2
  if (length(pairs$i) < expected) {
    ex <- .missing_pair(pairs$i, pairs$j, seq_len(n))
    stop(paste0(
      "Missing distances in `", name, "` for ", expected - length(pairs$i),
      " pair(s) of ", what, " (e.g. ", ids[ex[1]], "-", ids[ex[2]], ")"
    ))
  }
  invisible(NULL)
}


# Keeps the pairs of subunits of the same unit, and stops unless every pair of
# subunits of every unit is listed
.within_unit_pairs <- function(pairs, cl, ids, name = "diss") {
  within <- cl$unit_of[pairs$i] == cl$unit_of[pairs$j]
  pairs <- .subset_pairs(pairs, within)

  n_units <- length(cl$unit_ids)
  sizes <- tabulate(cl$unit_of, n_units)
  counts <- tabulate(cl$unit_of[pairs$i], n_units)
  short <- which(counts < sizes * (sizes - 1) / 2)
  if (length(short) > 0) {
    u <- short[1]
    members <- which(cl$unit_of == u)
    in_u <- cl$unit_of[pairs$i] == u
    ex <- .missing_pair(pairs$i[in_u], pairs$j[in_u], members)
    stop(paste0(
      "Missing distances in `", name, "` for ", sizes[u] * (sizes[u] - 1) / 2 - counts[u],
      " pair(s) of subunits in unit: ", cl$unit_ids[u],
      " (e.g. ", ids[ex[1]], "-", ids[ex[2]], ")"
    ))
  }

  pairs
}


# Scales every row of M by v
.row_scale <- function(M, v) {
  if (inherits(M, "Matrix")) {
    return(Matrix::Diagonal(x = v) %*% M)
  }
  M * v
}


# Relative abundances (rows summing to one; empty samples stay at zero)
.relative <- function(ab) {
  .row_scale(ab$A, ifelse(ab$N > 0, 1 / ab$N, 0))
}


# Indicator matrix (subunits x units) of the clustering
.unit_indicator <- function(unit_of, n_units) {
  Matrix::sparseMatrix(
    i = seq_along(unit_of),
    j = unit_of,
    x = 1,
    dims = c(length(unit_of), n_units)
  )
}


# Nonzero entries of a matrix: row index and value
.nonzero <- function(M) {
  if (inherits(M, "sparseMatrix")) {
    s <- Matrix::summary(methods::as(M, "CsparseMatrix"))
    nz <- s$x != 0
    return(list(i = s$i[nz], j = s$j[nz], x = s$x[nz]))
  }
  M <- as.matrix(M)
  idx <- which(M != 0, arr.ind = TRUE)
  list(i = idx[, 1], j = idx[, 2], x = M[idx])
}


# For every sample (row of X and Y), sum over the pairs k of X[, i_k] Y[, j_k] w_k.
# With `group`, the sums are kept apart for every group of pairs (one column per
# group). Pairs are read in chunks, so only a few columns are needed at a time
.pair_sums <- function(X, i, j, w, Y = X, group = NULL, n_groups = 1L) {
  n_rows <- nrow(X)
  out <- matrix(0, nrow = n_rows, ncol = if (is.null(group)) 1L else n_groups)
  m <- length(w)
  if (m == 0 || n_rows == 0) {
    return(out)
  }

  chunk <- max(1L, floor(.chunk_cells / n_rows))
  for (start in seq(1, m, by = chunk)) {
    r <- start:min(m, start + chunk - 1)
    prod <- X[, i[r], drop = FALSE] * Y[, j[r], drop = FALSE]
    if (is.null(group)) {
      out[, 1] <- out[, 1] + as.vector(as.matrix(prod %*% w[r]))
    } else {
      W <- Matrix::sparseMatrix(
        i = seq_along(r),
        j = group[r],
        x = w[r],
        dims = c(length(r), n_groups)
      )
      out <- out + as.matrix(prod %*% W)
    }
  }
  out
}


# Rao's quadratic entropy of every sample, sum_{i != j} p_i p_j d_ij, from the
# relative abundances P and a complete table of pairs
.rao_q <- function(P, pairs) {
  2 * .pair_sums(P, pairs$i, pairs$j, pairs$d)[, 1]
}


# Result named by sample, with NA for empty samples
.by_sample <- function(x, ab) {
  x <- as.numeric(x)
  x[ab$N <= 0] <- NA_real_
  names(x) <- ab$samples
  x
}


.check_sig <- function(sig, name = "sig") {
  if (!is.numeric(sig) || length(sig) != 1 || is.na(sig) || sig <= 0) {
    stop(paste0("`", name, "` must be a single positive numeric value"))
  }
}


.check_q <- function(q) {
  if (!is.numeric(q) || length(q) != 1 || !is.finite(q) || q < 0) {
    stop("`q` must be a single non-negative numeric value")
  }
}


.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(paste0("`", name, "` must be a single TRUE or FALSE value"))
  }
}


# Every pair of units at distance sig
.sigma_unit_pairs <- function(n_units, sig) {
  if (n_units < 2) {
    return(list(i = integer(0), j = integer(0), d = numeric(0)))
  }
  idx <- which(upper.tri(matrix(0, n_units, n_units)), arr.ind = TRUE)
  list(i = idx[, 1], j = idx[, 2], d = rep(sig, nrow(idx)))
}


# Distances between units, aggregated from the distances between their subunits
# while the pairs are read: one counter (min, max, or sum and count) per pair of
# units, merged with every chunk. `finish()` returns the pairs of units (indices
# i < j into the units) with their distances
.unit_distance_counter <- function(unit_of, n_units, method) {
  if (!(method %in% c("min", "max", "average"))) {
    stop(paste0("Cluster distance method: ", method, " not supported"))
  }
  n_units <- as.numeric(n_units)
  key <- numeric(0)
  val <- numeric(0)
  cnt <- numeric(0)

  update <- function(pairs) {
    u1 <- unit_of[pairs$i]
    u2 <- unit_of[pairs$j]
    cross <- u1 != u2
    if (!any(cross)) {
      return(invisible(NULL))
    }
    k <- c(key, pmin(u1, u2)[cross] + (pmax(u1, u2)[cross] - 1) * n_units)
    d <- c(val, pairs$d[cross])

    if (method == "average") {
      ug <- sort(unique(k))
      g <- match(k, ug)
      val <<- as.vector(rowsum(d, g, reorder = TRUE))
      cnt <<- as.vector(rowsum(c(cnt, rep(1, sum(cross))), g, reorder = TRUE))
      key <<- ug
    } else {
      o <- order(k, if (method == "min") d else -d)
      first <- !duplicated(k[o])
      key <<- k[o][first]
      val <<- d[o][first]
    }
    invisible(NULL)
  }

  finish <- function() {
    list(
      i = as.integer((key - 1) %% n_units + 1),
      j = as.integer((key - 1) %/% n_units + 1),
      d = if (method == "average") val / cnt else val
    )
  }

  list(update = update, finish = finish)
}


# Stops unless every pair of units has a distance. With `units` (indices), only
# the pairs of those units are required
.require_unit_pairs <- function(pairs_clust, unit_ids, units = seq_along(unit_ids)) {
  n_units <- length(units)
  expected <- n_units * (n_units - 1) / 2
  keep <- pairs_clust$i %in% units & pairs_clust$j %in% units
  pairs_clust <- .subset_pairs(pairs_clust, keep)
  if (length(pairs_clust$i) < expected) {
    ex <- units[.missing_pair(match(pairs_clust$i, units), match(pairs_clust$j, units), seq_len(n_units))]
    stop(paste0(
      "No distances in `diss` between the subunits of ", expected - length(pairs_clust$i),
      " pair(s) of clusters (e.g. ", unit_ids[ex[1]], "-", unit_ids[ex[2]], ")"
    ))
  }
  invisible(NULL)
}


#' Distances between units
#'
#' Computes the distance between every pair of units (clusters) from the pairwise distances
#' between their subunits, given as a three-column table. The result is itself a three-column
#' table, which can be passed as `diss_clust` to [multiplicity.distance()] with
#' `method = "custom"`.
#'
#' @param diss Data frame with three columns (taken by position): the identifiers of two subunits
#'   and the distance between them, or the path to a file with them (see [distance-files]).
#'   Every pair of subunits of different units must be listed (in either orientation). Ignored
#'   for `method = "sigma"`, where it can be `NULL`.
#' @param clust Named vector of cluster memberships: names are subunit identifiers and values
#'   are the units they belong to. A data frame is also accepted, with two columns (subunit
#'   identifiers, units) or one column of units with the subunit identifiers as row names.
#' @param method Character string specifying the linkage method for aggregating the distances
#'   between the subunits of two units. One of: "average", "min", "max", or "sigma". Sigma sets
#'   all distances between units to `sig`. Default is "average".
#' @param sig Numeric distance between units used by the "sigma" method. Default is 1.
#' @inheritParams raoQuadratic
#'
#' @return Data frame with columns `ID1`, `ID2` and `Distance`, one row per pair of units.
#'
#' @details
#' The distances are read sequentially: every row of `diss` updates a counter (minimum,
#' maximum or sum) of the pair of units of its subunits, so memory depends on the number of
#' pairs of units, not on the number of pairs of subunits. Distances are not capped.
#'
#' Only the pairs of subunits of different units are used: rows of pairs inside a unit are
#' ignored, also when their distance is invalid (NA, infinite or negative). With a distance file, only the
#' pairs of subunits of different units are checked by `check_large_distance_file`.
#'
#' @seealso [multiplicity.distance()], [distance-files]
#'
#' @export
#'
#' @examples
#' clust <- c(a1 = "A", a2 = "A", b1 = "B", c1 = "C")
#' diss <- data.frame(
#'   ID1 = c("a1", "a1", "a1", "a2", "a2", "b1"),
#'   ID2 = c("a2", "b1", "c1", "b1", "c1", "c1"),
#'   Distance = c(0.1, 0.5, 0.9, 0.7, 0.8, 0.6)
#' )
#' unit_distances(diss, clust, method = "average")
#' unit_distances(diss, clust, method = "min")
#'
unit_distances <- function(
  diss,
  clust,
  method = "average",
  sig = 1,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  clust <- .as_clust_vector(clust)
  if (is.null(names(clust))) {
    stop("`clust` must be a named vector (names are subunit identifiers)")
  }
  .check_sig(sig)
  if (!is.character(method) || length(method) != 1 ||
    !(method %in% c("sigma", "min", "max", "average"))) {
    stop(paste0("Cluster distance method: ", method, " not supported"))
  }
  ids <- .as_ids(names(clust))
  cl <- .parse_clust(clust, ids)
  n_units <- length(cl$unit_ids)

  if (method == "sigma") {
    ud <- .sigma_unit_pairs(n_units, sig)
  } else {
    src <- .diss_source(diss, ids, chunk_size, check_large_distance_file, "between", cl, header = header)
    counter <- .unit_distance_counter(cl$unit_of, n_units, method)
    acc <- list(
      scope = "between",
      update = counter$update,
      finish = counter$finish
    )
    ud <- .consume_pairs(src, list(acc), ids, cl)[[1]]
    .require_unit_pairs(ud, cl$unit_ids)
  }

  data.frame(
    ID1 = cl$unit_ids[ud$i],
    ID2 = cl$unit_ids[ud$j],
    Distance = ud$d,
    stringsAsFactors = FALSE
  )
}
