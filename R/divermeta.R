#' Estimate diversity and multiplicity indices across samples
#'
#' Estimate several diversity and multiplicity indices for many samples at once.
#' Abundances are provided with samples in rows and subunits (species/OTUs/genes)
#' in columns. Dissimilarities are provided as a three-column table of pairs of subunits.
#'
#' Indices available (use these labels in `indices`):
#' - "multiplicity_inventory": inventory multiplicity \eqn{^{q}M}{M^q} (order \eqn{q}{q})
#' - "multiplicity_distance": distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} (cutoff `sig`)
#' - "raoQ": Rao quadratic entropy \eqn{Q}{Q}
#' - "FD_sigma": distance-based functional diversity \eqn{\delta D_{\sigma}}{delta D_sigma} (cutoff `sig`)
#' - "FD_q": distance-based functional diversity \eqn{^{q}FD}{FD^q} (order \eqn{q}{q})
#' - "redundancy": functional redundancy \eqn{Re}{Re}
#' - "relative_multiplicity": relative multiplicity \eqn{RM}{RM} (cutoff `sig`)
#' - "average_redundancy": average normalized redundancy \eqn{\overline{Re}}{mean Re}
#'
#' Notes:
#' - Indices that use dissimilarities (`multiplicity_distance`, `raoQ`, `FD_sigma`,
#'   `FD_q`, `redundancy`, `relative_multiplicity`, `average_redundancy`) require `diss`.
#' - Indices that use clustering (`multiplicity_inventory`, `multiplicity_distance`,
#'   `relative_multiplicity`, `average_redundancy`) require `clust`.
#' - `relative_multiplicity` and `average_redundancy` are computed as in
#'   [relative.multiplicity()] and [average.redundancy()] (with `normalize = TRUE`), which also
#'   offer a cutoff per unit, known reference diversities ([relative.multiplicity.ref_div()]) and
#'   the unnormalized redundancy. Subunits listed in a named `clust` that are not in `abund` have zero abundance
#'   in every sample: they count in the units of these two indices and do not change the others.
#' - Every entry of `indices` gives one column, named after the index (aliases such as
#'   `"M_inventory"` are renamed to `"multiplicity_inventory"`). An index requested more than once,
#'   by the same label or through an alias, gives one identical column per request, and the
#'   repeated columns get unique names from [data.frame()] (e.g. `raoQ` and `raoQ.1`).
#'
#' @param abund Numeric matrix, data.frame or sparse `Matrix` of abundances with samples in rows
#'   and subunits in columns. Column names are the subunit identifiers.
#' @param diss Optional data frame with three columns (taken by position): the identifiers of two
#'   subunits and the dissimilarity between them, or the path to a file with them (see
#'   [distance-files]). Every pair of subunits present in at least one sample must be listed
#'   (see [multiplicity.distance()] for the pairs needed with `method = "sigma"`); pairs with a
#'   subunit that has zero abundance in every sample may be missing, with a warning, except when
#'   `relative_multiplicity` is requested with `assume_homogeneous_abundance = TRUE` and
#'   `assume_max_reference_distance = FALSE`, where every subunit counts in the reference. A file
#'   is read only once for all the requested indices.
#' @param indices Character vector of index names to compute (see list above).
#' @param clust Optional vector/factor of cluster memberships for each subunit. If named, names
#'   are subunit identifiers; otherwise it follows the columns of `abund`. A data frame is also accepted, with two columns (subunit identifiers, units) or one
#'   column of units (named by its row names, if set).
#' @param q Numeric order for Hill-number based indices (used by `FD_q` and
#'   `multiplicity_inventory`). Default `1`.
#' @param sig Numeric cutoff `σ` for distance-based measures (used by
#'   `FD_sigma`, `multiplicity_distance` and `relative_multiplicity`). Default `1`.
#' @param method Character string specifying the linkage method for
#'   aggregating distances between clusters. One of: "average", "min",
#'   "max", "sigma" or "custom". Sigma sets all distances between clusters to be
#'    of value `sig`. Default is "sigma".
#' @param diss_clust Optional three-column table of distances between clusters, used by
#'   `multiplicity_distance` when `method` is "custom".
#' @param include_absent Logical (default `FALSE`), used by `relative_multiplicity` and
#'   `average_redundancy`. If `FALSE`, each sample averages only the units present in it. If
#'   `TRUE`, every unit counts in the average, and units absent from the sample contribute 0.
#' @param cap_at_one Logical (default `FALSE`), used by `relative_multiplicity`. If `TRUE`, each
#'   unit's ratio is capped at 1.
#' @param assume_max_reference_distance Logical (default `FALSE`), used by
#'   `relative_multiplicity`. If `TRUE`, every pair of subunits in the reference of a unit is at
#'   distance `sig`, instead of their listed distance.
#' @param assume_homogeneous_abundance Logical (default `FALSE`), used by
#'   `relative_multiplicity`. If `TRUE`, all subunits of a unit have the same abundance in its
#'   reference, instead of their abundances pooled across all samples.
#' @inheritParams raoQuadratic
#' @param normalize Logical indicating whether index values should be normalized to \[0, 1\]
#'   by dividing each column by its maximum value. Useful for comparing indices with
#'   different scales. Default `FALSE`.
#'
#' @return data.frame with one row per sample: a first column `Sample` with the sample identifiers
#'   (the row names of `abund`, or the row numbers when it has none), then one column per requested
#'   index. Empty samples (all zero abundances) return `NA` for all indices.
#'
#' @seealso [multiplicity.inventory()], [multiplicity.distance()], [diversity.functional()],
#'   [diversity.functional.traditional()], [raoQuadratic()], [redundancy()],
#'   [relative.multiplicity()], [average.redundancy()]
#'
#' @examples
#'
#' # Abundance table: rows = samples, columns = subunits
#' abund <- matrix(
#'   c(
#'     10, 0, 4, 6, # Sample1
#'     5,  8, 3, 0, # Sample2
#'     0, 12, 7, 9  # Sample3
#'   ),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(c("Sample1", "Sample2", "Sample3"), c("s1", "s2", "s3", "s4"))
#' )
#'
#' # Dissimilarities among subunits (0-1), one row per pair
#' diss <- data.frame(
#'   ID1 = c("s1", "s1", "s1", "s2", "s2", "s3"),
#'   ID2 = c("s2", "s3", "s4", "s3", "s4", "s4"),
#'   Distance = c(0.2, 0.7, 1.0, 0.6, 0.9, 0.4)
#' )
#'
#' # Clusters per subunit (used by multiplicity_* indices)
#' clust <- c(s1 = "A", s2 = "A", s3 = "B", s4 = "B")
#'
#' ## Run divermeta
#' divermeta(abund,
#'   clust = clust,
#'   diss = diss,
#'   q = 1,
#'   sig = 0.8,
#'   indices = c("multiplicity_inventory", "multiplicity_distance", "raoQ", "redundancy")
#' )
#'
#' ## Indices that average over the clusters
#' divermeta(abund,
#'   clust = clust,
#'   diss = diss,
#'   sig = 0.8,
#'   indices = c("relative_multiplicity", "average_redundancy")
#' )
#'
#' @export
#'
divermeta <- function(
    abund,
    diss = NULL,
    indices = c("multiplicity_inventory"),
    clust = NULL,
    q = 1,
    sig = 1,
    method = "sigma",
    diss_clust = NULL,
    include_absent = FALSE,
    cap_at_one = FALSE,
    assume_max_reference_distance = FALSE,
    assume_homogeneous_abundance = FALSE,
    normalize = FALSE,
    chunk_size = 1e6,
    check_large_distance_file = FALSE,
    header = TRUE) {
  ## Validation
  if (!is.matrix(abund) && !is.data.frame(abund) && !inherits(abund, "Matrix")) {
    stop("abund must be a matrix, data.frame or sparse Matrix")
  }
  if (!is.character(indices) || length(indices) == 0) {
    stop("indices must be a non-empty character vector")
  }
  .check_flag(normalize, "normalize")

  ## Normalize and map index aliases
  idx_map <- list(
    raoQ                   = "raoQ",
    FD_sigma               = "FD_sigma",
    FD_q                   = "FD_q",
    FDq                    = "FD_q",
    redundancy             = "redundancy",
    multiplicity_inventory = "multiplicity_inventory",
    M_inventory            = "multiplicity_inventory",
    multiplicity_distance  = "multiplicity_distance",
    M_distance             = "multiplicity_distance",
    relative_multiplicity  = "relative_multiplicity",
    RM                     = "relative_multiplicity",
    average_redundancy     = "average_redundancy",
    AR                     = "average_redundancy"
  )
  normalized_indices <- vapply(indices, function(x) {
    if (!is.null(idx_map[[x]])) idx_map[[x]] else x
  }, character(1))

  supported <- unname(unlist(idx_map))
  unsupported <- setdiff(normalized_indices, unique(supported))
  if (length(unsupported) > 0) {
    stop(paste0("Unsupported indices: ", paste(unsupported, collapse = ", ")))
  }

  # Indices that average over the units
  by_unit <- c("relative_multiplicity", "average_redundancy")
  with_rm <- "relative_multiplicity" %in% normalized_indices

  needs_diss <- any(normalized_indices %in% c("raoQ", "FD_sigma", "FD_q", "redundancy", "multiplicity_distance", by_unit))
  needs_clust <- any(normalized_indices %in% c("multiplicity_inventory", "multiplicity_distance", by_unit))

  if (needs_diss && is.null(diss)) {
    stop("A dissimilarity table `diss` is required for the requested indices")
  }
  if (needs_clust && is.null(clust)) {
    stop("`clust` must be provided for the requested indices")
  }
  if (any(normalized_indices %in% c("FD_q", "multiplicity_inventory"))) {
    .check_q(q)
  }
  if (any(normalized_indices %in% c("FD_sigma", "multiplicity_distance", "relative_multiplicity"))) {
    .check_sig(sig)
  }
  if (any(normalized_indices %in% by_unit)) {
    .check_flag(include_absent, "include_absent")
  }
  if (with_rm) {
    .check_flag(cap_at_one, "cap_at_one")
    .check_flag(assume_max_reference_distance, "assume_max_reference_distance")
    .check_flag(assume_homogeneous_abundance, "assume_homogeneous_abundance")
  }
  if ("multiplicity_distance" %in% normalized_indices) {
    .check_multiplicity_distance_args(method, sig, diss_clust)
  }

  ## Parse the shared inputs once
  ab <- .parse_abund(abund)
  if (any(normalized_indices %in% by_unit)) {
    # Subunits only listed in `clust` count in the units, with zero abundance
    ab <- .extend_abund(ab, .as_clust_vector(clust))
  }
  cl <- if (needs_clust) .parse_clust(clust, ab$subunits) else NULL

  ## Accumulators of the indices that use distances, fed by a single read of them
  acc_map <- list(
    raoQ = function() .raoQuadratic_acc(ab),
    FD_sigma = function() .diversity.functional_acc(ab, sig),
    FD_q = function() .diversity.functional.traditional_acc(ab, q),
    redundancy = function() .redundancy_acc(ab),
    multiplicity_distance = function() {
      .multiplicity.distance_acc(ab, cl, method, sig, diss_clust)
    },
    relative_multiplicity = function() {
      .relative.multiplicity_acc(
        ab, cl, rep(sig, length(cl$unit_ids)), cap_at_one, include_absent,
        assume_max_reference_distance, assume_homogeneous_abundance
      )
    },
    average_redundancy = function() .average.redundancy_acc(ab, cl, TRUE, include_absent)
  )
  with_diss <- unique(normalized_indices[normalized_indices %in% names(acc_map)])
  res_diss <- list()
  if (length(with_diss) > 0) {
    accs <- lapply(acc_map[with_diss], function(make) make())
    scopes <- vapply(accs, `[[`, "", "scope")
    scope <- if ("all" %in% scopes) "all" else "within"
    # The homogeneous reference with the listed distances uses every subunit,
    # also those with zero abundance in every sample
    need_all <- with_rm && assume_homogeneous_abundance && !assume_max_reference_distance
    src <- .diss_source(
      diss, ab$subunits, chunk_size, check_large_distance_file, scope, cl,
      header = header, need = if (need_all) NULL else .present_subunits(ab)
    )
    res_diss <- .consume_pairs(src, accs, ab$subunits, cl)
    names(res_diss) <- with_diss
  }

  ## Compute every requested index for all samples
  res_list <- lapply(normalized_indices, function(idx) {
    if (idx == "multiplicity_inventory") {
      return(as.numeric(.multiplicity.inventory_core(ab, cl, q)))
    }
    as.numeric(res_diss[[idx]])
  })
  names(res_list) <- normalized_indices

  res <- do.call(cbind, res_list)

  # Normalize if requested
  if (normalize) {
    # Columns with only NA (every sample empty) have no maximum
    col_max <- apply(res, 2, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    # Only normalize columns with non-zero maximum
    nonzero_cols <- col_max != 0 & !is.na(col_max)

    if (any(nonzero_cols)) {
      res[, nonzero_cols] <- sweep(
        res[, nonzero_cols, drop = FALSE],
        2,
        col_max[nonzero_cols], "/"
      )
    }
  }

  res <- data.frame(Sample = ab$samples, res)
  return(res)
}
