divermeta <- function(abund,
  diss = NULL,
  indices = c("multiplicity_inventory"),
  clusters = NULL,
  q = 1,
  sig = 1) {

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
    if (!is.matrix(diss) || nrow(diss) != ncol(diss)) {
      stop("`diss` must be a square numeric matrix")
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
      if (nrow(abund) != nrow(diss)) {
        stop("`diss` dimensions do not match `abund` and lack dimnames for alignment")
      }
      diss_mtx <- diss
    }

    ## Build an upper-triangle three-column frame for multiplicity_distance
    if (any(normalized_indices == "multiplicity_distance")) {
      ut_idx <- which(upper.tri(diss_mtx), arr.ind = TRUE)
      diss_frame <- data.frame(
        ID1 = feature_names[ut_idx[, 1]],
        ID2 = feature_names[ut_idx[, 2]],
        Distance = diss_mtx[upper.tri(diss_mtx)],
        stringsAsFactors = FALSE
      )
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
      if (sum(ab_vec) <= 0) return(NA_real_)
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
    compute_map$multiplicity_distance <- guard_nonempty(function(ab_vec) multiplicity.distance.by_blocks(
      ids = ids,
      ab = ab_vec,
      diss_frame = diss_frame,
      clust = clust_vec,
      sigma = sig
    ))
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
  res <- data.frame(Sample = colnames(abund), res)
  return(res)
}

