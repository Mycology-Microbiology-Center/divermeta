#' Estimate diversity and multiplicity indices across samples
#'
#' Estimate several diversity and multiplicity indices for many samples at once.
#' Abundances are provided with samples in columns and features (species/OTUs/genes)
#' in rows. If supplied, the dissimilarity matrix must correspond to feature names.
#'
#' Indices available (use these labels in `indices`):
#' - "multiplicity_inventory": inventory multiplicity \eqn{^{q}M}{M^q} (order \eqn{q}{q})
#' - "multiplicity_distance": distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} (cutoff `sig`)
#' - "raoQ": Rao quadratic entropy \eqn{Q}{Q}
#' - "FD_sigma": distance-based functional diversity \eqn{\delta D_{\sigma}}{delta D_sigma} (cutoff `sig`)
#' - "FD_q": distance-based functional diversity \eqn{^{q}FD}{FD^q} (order \eqn{q}{q})
#' - "redundancy": functional redundancy \eqn{Re}{Re}
#'
#' Notes:
#' - Indices that use dissimilarities (`multiplicity_distance`, `raoQ`, `FD_sigma`,
#'   `FD_q`, `redundancy`) require `diss`.
#' - Indices that use clustering (`multiplicity_inventory`, `multiplicity_distance`)
#'   require `clusters` (same length/order as feature rows).
#'
#' @param abund Numeric matrix or data.frame of abundances with samples in columns
#'   and features in rows.
#' @param diss Optional numeric dissimilarity matrix or `dist` object among features; square with
#'   row/column names matching `rownames(abund)` when present.
#' @param indices Character vector of index names to compute (see list above).
#' @param clusters Optional vector/factor of cluster memberships for each feature
#'   (row of `abund`).
#' @param q Numeric order for Hill-number based indices (used by `FD_q` and
#'   `multiplicity_inventory`). Default `1`.
#' @param sig Numeric cutoff `σ` for distance-based measures (used by
#'   `FD_sigma` and `multiplicity_distance`). Default `1`.
#' @param normalize boolean value indicating if values should be normalized between 0 and 1. Values will be
#'    normalized by dividing all index columns by their corresponding max value. Default `FALSE`
#'
#' @return data.frame with one row per sample and one column per requested index.
#'
#' @seealso [multiplicity.inventory()], [multiplicity.distance()], [diversity.functional()],
#'   [diversity.functional.traditional()], [raoQuadratic()], [redundancy()]
#'
#' @examples
#'
#' # Abundance matrix: rows = features, columns = samples
#' species <- c("s1", "s2", "s3", "s4")
#' samples <- c("Sample1", "Sample2", "Sample3")
#' abund <- matrix(
#'   c(
#'     10, 5,  0, # s1
#'     0,  8, 12, # s2
#'     4,  3,  7, # s3
#'     6,  0,  9 # s4
#'   ),
#'   nrow = length(species),
#'   ncol = length(samples),
#'   byrow = TRUE,
#'   dimnames = list(species, samples)
#' )
#'
#' # Dissimilarity matrix among features (0-1), names must match abund rownames
#' diss <- matrix(
#'   c(
#'     0.0, 0.2, 0.7, 1.0,
#'     0.2, 0.0, 0.6, 0.9,
#'     0.7, 0.6, 0.0, 0.4,
#'     1.0, 0.9, 0.4, 0.0
#'   ),
#'   nrow = length(species),
#'   ncol = length(species),
#'   byrow = TRUE,
#'   dimnames = list(species, species)
#' )
#'
#' # Clusters per feature (used by multiplicity_* indices)
#' clusters <- c(s1 = "A", s2 = "A", s3 = "B", s4 = "B")
#'
#' ## Run divermeta
#' divermeta(abund,
#'   clusters = clusters,
#'   diss = diss,
#'   q = 1,
#'   sig = 0.8,
#'   indices = c("multiplicity_inventory", "multiplicity_distance", "raoQ", "redundancy")
#' )
#'
#' @export
#'
divermeta <- function(
    abund,
    diss = NULL,
    indices = c("multiplicity_inventory"),
    clusters = NULL,
    q = 1,
    sig = 1,
    normalize = FALSE) {
  ## Validation
  if (!is.matrix(abund) && !is.data.frame(abund)) {
    stop("abund must be a matrix or data.frame")
  }
  abund <- as.matrix(abund)
  if (!is.numeric(abund)) {
    stop("abund must be numeric")
  }
  if (any(is.na(abund))) {
    stop("abund contains NA values; please impute or remove them")
  }

  if (!is.character(indices) || length(indices) == 0) {
    stop("indices must be a non-empty character vector")
  }

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
    M_distance             = "multiplicity_distance"
  )
  normalized_indices <- vapply(indices, function(x) {
    if (!is.null(idx_map[[x]])) idx_map[[x]] else x
  }, character(1))

  supported <- unname(unlist(idx_map))
  unsupported <- setdiff(normalized_indices, unique(supported))
  if (length(unsupported) > 0) {
    stop(paste0("Unsupported indices: ", paste(unsupported, collapse = ", ")))
  }

  needs_diss <- any(normalized_indices %in% c("raoQ", "FD_sigma", "FD_q", "redundancy", "multiplicity_distance"))
  needs_clusters <- any(normalized_indices %in% c("multiplicity_inventory", "multiplicity_distance"))

  ## Dissimilarity alignment
  diss_mtx <- NULL
  diss_frame <- NULL
  feature_names <- rownames(abund)

  if (needs_diss) {
    if (is.null(diss)) {
      stop("A dissimilarity matrix `diss` is required for the requested indices")
    }
    if (!is.numeric(diss)) {
      stop("`diss` must be numeric")
    }

    ## Try aligning by rownames if available
    if (!is.null(rownames(diss)) && !is.null(colnames(diss))) {
      if (!identical(sort(rownames(diss)), sort(colnames(diss)))) {
        stop("`diss` row and column names must match")
      }

      if (!is.null(feature_names)) {
        common <- intersect(feature_names, rownames(diss))
        if (length(common) == 0) {
          stop("No overlap between `rownames(abund)` and `dimnames(diss)`")
        }

        if (length(common) < nrow(abund)) {
          warning("Some features in `abund` are missing from `diss`; dropping them")
        }
        if (length(common) < nrow(diss)) {
          warning("Some features in `diss` are missing from `abund`; subsetting dissimilarities")
        }

        ## Subset and reorder both to common
        abund <- abund[common, , drop = FALSE]
        diss_mtx <- diss[common, common, drop = FALSE]
        feature_names <- common
      } else {
        ## No rownames in abund: enforce dimension match
        if (nrow(abund) != nrow(diss)) {
          stop("When `abund` has no rownames, `diss` dimensions must match it exactly")
        }
        diss_mtx <- diss
      }
    } else {
      ## No dimnames on diss; assume it matches abund as-is
      # If dist object
      if (inherits(diss, "dist")) {
        n <- attr(diss, "Size")
        if (n != nrow(abund)) {
          stop(paste0(
            "Lack dimnames for alignment so abundance vector and matrix must have compatible sizes. Matrix: ",
            n, "x", n, ". Vector: ", nrow(abund)
          ))
        }
      } else {
        dims <- dim(diss)
        n <- dims[1]
        if (dims[1] != dims[2]) {
          stop(paste0(
            "Lack dimnames for alignment so distance matrix must be square. Matrix: ",
            n, "x", n, "."
          ))
        }

        if (dims[1] != nrow(abund)) {
          stop(paste0(
            "Lack dimnames for alignment so abundance vector and matrix must have compatible sizes. Matrix: ",
            dims[1], "x", dims[2], ". Vector: ", nrow(abund)
          ))
        }
      }
      diss_mtx <- diss
    }



    ## Build an upper-triangle three-column frame for multiplicity_distance
    if (any(normalized_indices == "multiplicity_distance")) {
      if (inherits(diss_mtx, "dist")) {
        # For dist object, use combn to get all pairs
        ut_idx <- combn(n, 2)
        diss_frame <- data.frame(
          ID1 = feature_names[ut_idx[1, ]],
          ID2 = feature_names[ut_idx[2, ]],
          Distance = as.vector(diss_mtx),
          stringsAsFactors = FALSE
        )
      } else {
        ut_idx <- which(upper.tri(diss_mtx), arr.ind = TRUE)
        diss_frame <- data.frame(
          ID1 = feature_names[ut_idx[, 1]],
          ID2 = feature_names[ut_idx[, 2]],
          Distance = diss_mtx[upper.tri(diss_mtx)],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  ## Cluster alignment
  clust_vec <- NULL
  if (needs_clusters) {
    if (is.null(clusters)) {
      stop("`clusters` must be provided for the requested indices")
    }
    if (length(clusters) != nrow(as.matrix(clusters))) {
      ## Protect against matrices accidentally passed
      clusters <- as.vector(clusters)
    }
    if (!is.null(names(clusters)) && !is.null(feature_names)) {
      missing <- setdiff(feature_names, names(clusters))
      if (length(missing) > 0) {
        stop("`clusters` is named but missing some features present in `abund`/`diss`")
      }
      clust_vec <- clusters[feature_names]
    } else {
      if (length(clusters) != nrow(abund)) {
        stop("`clusters` length must equal number of features (rows of `abund`)")
      }
      clust_vec <- clusters
    }
  }

  ## Helper to guard against empty (all-zero) samples
  guard_nonempty <- function(f) {
    function(ab_vec) {
      if (sum(ab_vec) <= 0) {
        return(NA_real_)
      }
      f(ab_vec)
    }
  }

  ## Build compute functions per requested index
  compute_map <- list()

  if ("raoQ" %in% normalized_indices) {
    compute_map$raoQ <- guard_nonempty(function(ab_vec) raoQuadratic(ab_vec, diss_mtx))
  }
  if ("FD_sigma" %in% normalized_indices) {
    compute_map$FD_sigma <- guard_nonempty(function(ab_vec) diversity.functional(ab_vec, diss_mtx, sig = sig))
  }
  if ("FD_q" %in% normalized_indices) {
    compute_map$FD_q <- guard_nonempty(function(ab_vec) diversity.functional.traditional(ab_vec, diss_mtx, q = q))
  }
  if ("redundancy" %in% normalized_indices) {
    compute_map$redundancy <- guard_nonempty(function(ab_vec) redundancy(ab_vec, diss_mtx))
  }
  if ("multiplicity_inventory" %in% normalized_indices) {
    compute_map$multiplicity_inventory <- guard_nonempty(function(ab_vec) multiplicity.inventory(ab_vec, clust_vec, q = q))
  }
  if ("multiplicity_distance" %in% normalized_indices) {
    ids <- feature_names
    compute_map$multiplicity_distance <- guard_nonempty(function(ab_vec) {
      multiplicity.distance.by_blocks(
        ids = ids,
        ab = ab_vec,
        diss_frame = diss_frame,
        clust = clust_vec,
        sigma = sig
      )
    })
  }

  ## Evaluate per index using
  res_list <- lapply(normalized_indices, function(idx) {
    fn <- compute_map[[idx]]
    if (is.null(fn)) stop(paste0("Unexpected index label: ", idx))
    vals <- apply(abund, 2, fn)
    as.numeric(vals)
  })
  names(res_list) <- normalized_indices

  res <- do.call(cbind, res_list)

  # Checks if needs normalization
  if (normalize == TRUE) {
    col_max <- apply(res, 2, max, na.rm = TRUE)
    # If columns are all zero will not divide
    nonzero_cols <- col_max != 0 & !is.na(col_max)

    res[, nonzero_cols] <- sweep(
      res[, nonzero_cols, drop = FALSE],
      2,
      col_max[nonzero_cols], "/"
    )

    res
  }

  res <- data.frame(Sample = colnames(abund), res)
  return(res)
}
