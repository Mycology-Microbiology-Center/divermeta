#' Distance-based multiplicity
#'
#' Computes distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} as the ratio of distance-based
#' diversity before vs after clustering, using a cutoff \eqn{\sigma}{sigma} to cap pairwise distances.
#' This index quantifies the diversity lost when clustering elements into operational units, incorporating
#' pairwise dissimilarities (e.g., genetic, functional, or phylogenetic distances) among elements.
#'
#' @param ab Numeric vector of element abundances before clustering.
#' @param diss Numeric matrix or `dist` object of pairwise dissimilarities among elements
#'   before clustering. Must be square and match the length of `ab`.
#' @param clust Vector or factor of cluster memberships for each element of `ab`.
#'   Must have the same length as `ab`.
#' @param method Character string specifying the linkage method for
#'   aggregating distances between clusters. One of: "average", "min",
#'   "max", or "sigma". Sigma sets all distances between clusters to be
#'    of value `sig`. Default is "sigma".
#' @param sig Numeric cutoff \eqn{\sigma}{sigma} (default `1`) defining the threshold distance
#'   at which two units are considered maximally different. All distances greater than `sig`
#'   are capped at `sig`. This parameter should be chosen based on the biological meaning
#'   of distances in your dataset (e.g., maximum expected genetic distance, functional
#'   dissimilarity threshold).
#' @param clust_ids_order Integer or factor vector specifying the order of
#'   cluster IDs as they appear in the output. If `NULL` (default), uses the natural
#'   order of unique cluster IDs found in `clust`.
#' @param diss_clust Numeric matrix or `dist` object of pairwise dissimilarities among clusters
#'   after clustering. Only used when `method` is "custom". Default is `NULL`.
#'
#' @return Numeric scalar, the distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma}.
#'   This value represents the effective number of functionally distinct elements per cluster.
#'   A value of 1 indicates minimal functional diversity within clusters, while higher values
#'   indicate greater functional diversity lost through clustering.
#'
#' @details
#' Unlike inventory multiplicity, distance-based multiplicity is only defined for \eqn{q = 1}{q = 1}
#' (Shannon-type weighting), meaning all elements are weighted proportionally to their abundance.
#'
#' The function supports several methods for computing inter-cluster distances:
#' \itemize{
#'   \item \code{"sigma"}: All inter-cluster distances are set to \code{sig}, creating maximally
#'     different clusters (default)
#'   \item \code{"average"}: Mean distance between all cross-cluster element pairs
#'   \item \code{"min"}: Minimum distance between any cross-cluster element pair
#'   \item \code{"max"}: Maximum distance between any cross-cluster element pair
#'   \item \code{"custom"}: Use a pre-computed cluster distance matrix provided via \code{diss_clust}
#' }
#'
#' When using the "custom" method, you must provide \code{diss_clust} and ensure its dimensions
#' match the number of clusters. The \code{clust_ids_order} parameter determines how cluster
#' abundances are aligned with the distance matrix.
#'
#' @seealso [multiplicity.distance.by_blocks()] for computing from a compact distance table,
#'   [raoQuadratic()] for Rao's quadratic entropy,
#'   [diversity.functional()] for distance-based functional diversity,
#'   [cluster_distance_matrix()] for computing cluster distance matrices
#'
#' @export
#'
#' @examples
#' # Example 1: Using sigma method (default)
#' # Two clusters with maximally different clusters
#' ab <- c(10, 10, 10, 10)  # 4 elements
#' clust <- c(1, 1, 2, 2)    # 2 clusters
#' diss <- matrix(c(
#'   0.0, 0.3, 1.0, 1.0,
#'   0.3, 0.0, 1.0, 1.0,
#'   1.0, 1.0, 0.0, 0.4,
#'   1.0, 1.0, 0.4, 0.0
#' ), nrow = 4, byrow = TRUE)
#'
#' multiplicity.distance(ab, diss, clust, method = "sigma", sig = 1)
#'
#' # Example 2: Using average linkage
#' # Compute cluster distances as mean of cross-cluster element distances
#' multiplicity.distance(ab, diss, clust, method = "average", sig = 1)
#'
#' # Example 3: Using custom cluster distances
#' # Provide pre-computed cluster distance matrix
#' ab <- c(10, 10, 10, 10)
#' clust <- c("A", "A", "B", "B")
#' clust_ids <- c("A", "B")
#' diss <- dist(matrix(rnorm(4 * 3), ncol = 3))
#' diss_clust <- matrix(c(0, 0.8, 0.8, 0), nrow = 2)
#'
#' multiplicity.distance(
#'   ab, diss, clust,
#'   method = "custom",
#'   diss_clust = diss_clust,
#'   clust_ids_order = clust_ids
#' )
#'
#' # Example 4: Low diversity within clusters
#' # Elements within clusters are functionally similar (distance ~0)
#' diss_low <- matrix(c(
#'   0.0, 0.05, 1.0, 1.0,
#'   0.05, 0.0, 1.0, 1.0,
#'   1.0, 1.0, 0.0, 0.05,
#'   1.0, 1.0, 0.05, 0.0
#' ), nrow = 4, byrow = TRUE)
#'
#' multiplicity.distance(ab, diss_low, clust, method = "sigma", sig = 1)
#'
multiplicity.distance <- function(
  ab,
  diss,
  clust,
  method = "sigma",
  sig = 1,
  clust_ids_order = NULL,
  diss_clust = NULL
) {
  # Input validation
  if (!is.numeric(ab)) {
    stop("Abundance vector must be numeric")
  }
  if (!is.numeric(sig) || length(sig) != 1 || sig <= 0) {
    stop("`sig` must be a single positive numeric value")
  }

  if (sum(ab) == 0) {
    stop("Total abundance cannot be zero")
  }

  if (length(ab) != length(clust)) {
    stop("`ab` and `clust` must have the same length")
  }

  #  Computes the cluster ids
  if (method != "custom") {
    clust_ids_order <- unique(clust)
  }

  n_clust <- length(clust_ids_order)

  # Checks custom
  if (method == "custom") {
    if (is.null(diss_clust)) {
      stop("diss_clust cannot be NULL if method is 'custom'")
    }

    # If dist object
    if (inherits(diss_clust, "dist")) {
      n <- attr(diss_clust, "Size")
      if (n != length(clust_ids_order)) {
        stop(paste0(
          "Abundance vector and matrix must have compatible sizes. Matrix: ",
          n, "x", n, ". Vector: ", length(ab)
        ))
      }
    } else {
      dims <- dim(diss_clust)
      if (dims[1] != dims[2]) {
        stop(paste0(
          "Distance matrix must be square. Matrix: ",
          n_clust, "x", n_clust, "."
        ))
      }

      if (dims[1] != length(clust_ids_order)) {
        stop(paste0(
          "Abundance vector and matrix must have compatible sizes. Matrix: ",
          dims[1], "x", dims[2], ". Vector: ", length(clust_ids_order)
        ))
      }
    }


    # Converts to dist object
    if (!inherits(diss_clust, "dist")) {
      diss_clust <- as.dist(diss_clust)
    }

    n <- attr(diss_clust, "Size")

    if (n != n_clust) {
      stop(paste0(
        "Custom cluster distance matrix is not compatible with clusters. Matrix: ",
        n,
        "x",
        n,
        ". Cluster IDs order: ",
        n_clust
      ))
    }
  }

  # Check for empty input
  if (length(ab) == 0 || length(clust) == 0) {
    return(0)
  }

  #  Computes the cluster distances
  if (method != "custom") {
    diss_clust <- cluster_distance_matrix(
      diss = diss,
      clust = clust,
      clust_ids_order = clust_ids_order,
      method = method,
      sig = sig
    )
  }

  #  Computes the ab_clust
  ab_clust <- tapply(ab, clust, sum)

  # Reorder to match clust_ids_order
  ab_clust <- ab_clust[clust_ids_order]

  # Cap distances at sigma (modify copies to avoid side effects)
  diss_clust[diss_clust > sig] <- sig
  diss[diss > sig] <- sig

  # Compute Rao's quadratic entropy
  Q_before <- raoQuadratic(ab, diss)
  Q_after <- raoQuadratic(ab_clust, diss_clust)

  # Check for edge case where denominator would be zero
  if (abs(sig - Q_before) < .Machine$double.eps) {
    stop(
      "Cannot compute multiplicity: denominator (sigma - Q_before) is effectively zero"
    )
  }

  # Compute multiplicity
  (sig - Q_after) / (sig - Q_before)
}


#' Distance-based multiplicity by blocks
#'
#' Computes distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} from a compact
#' three-column distance table. This is an efficient implementation for large datasets where
#' storing the full distance matrix would be memory-intensive. Distances between elements from
#' different clusters are assumed to be equal to \eqn{\sigma}{sigma} (maximally different).
#'
#' @param ids Character or integer vector of element identifiers (length `n`). Must match
#'   the identifiers used in `diss_frame`.
#' @param ab Numeric vector of pre-clustering abundances (length `n`, same order as `ids`).
#' @param diss_frame Data frame with columns `ID1`, `ID2`, `Distance` containing pairwise
#'   distances for unique pairs of elements. Should only include within-cluster pairs or
#'   pairs where distances are less than `sigma`; cross-cluster distances are automatically
#'   set to `sigma`.
#' @param clust Vector or factor of cluster memberships for each element (length `n`, same order as `ids`).
#' @param sigma Numeric cutoff \eqn{\sigma}{sigma} (default `1`) at which two units are
#'   considered maximally different. All distances in `diss_frame` greater than `sigma` are
#'   capped, and distances between elements from different clusters are set to `sigma`.
#'
#' @return Numeric scalar, the distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma}.
#'
#' @seealso [multiplicity.distance()] for the standard implementation using full matrices,
#'   [raoQuadratic()] for Rao's quadratic entropy,
#'   [diversity.functional()] for distance-based functional diversity
#'
#' @export
#'
#' @examples
#' # Example: Compute multiplicity from a distance table
#' ids <- c("elem1", "elem2", "elem3", "elem4")
#' ab <- c(10, 15, 20, 25)
#' clust <- c(1, 1, 2, 2)
#'
#' # Distance table: only within-cluster pairs (cross-cluster assumed = sigma)
#' diss_frame <- data.frame(
#'   ID1 = c("elem1", "elem3"),
#'   ID2 = c("elem2", "elem4"),
#'   Distance = c(0.3, 0.4),
#'   stringsAsFactors = FALSE
#' )
#'
#' multiplicity.distance.by_blocks(ids, ab, diss_frame, clust, sigma = 1)
#'
multiplicity.distance.by_blocks <- function(
  ids,
  ab,
  diss_frame,
  clust,
  sigma = 1
) {
  # Input validation
  if (length(ids) == 0) {
    return(0)
  }
  if (length(ids) != length(ab) || length(ids) != length(clust)) {
    stop("`ids`, `ab`, and `clust` must have the same length")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
    stop("`sigma` must be a single positive numeric value")
  }
  if (
    !is.data.frame(diss_frame) ||
      !all(c("ID1", "ID2", "Distance") %in% colnames(diss_frame))
  ) {
    # Try to rename if columns exist but have different names
    if (ncol(diss_frame) == 3) {
      colnames(diss_frame) <- c("ID1", "ID2", "Distance")
    } else {
      stop(
        "`diss_frame` must be a data frame with columns `ID1`, `ID2`, `Distance`"
      )
    }
  }

  # Compute clustered abundances
  ab_clust <- tapply(ab, clust, sum)
  p_clust <- ab_clust / sum(ab_clust)

  # Normalize abundances
  p <- ab / sum(ab)
  names(p) <- ids

  # Cap distances at sigma
  diss_frame$Distance[diss_frame$Distance > sigma] <- sigma

  # Builds matrix multiplication by join
  diss_block <- diss_frame[diss_frame$ID1 != diss_frame$ID2, ]
  diss_block <- merge(
    diss_block,
    data.frame(ID1 = ids, Abundance1 = p),
    by = "ID1",
    all = FALSE
  )
  diss_block <- merge(
    diss_block,
    data.frame(ID2 = ids, Abundance2 = p),
    by = "ID2",
    all = FALSE
  )

  # Compute Rao's quadratic entropy before clustering
  # This accounts for within-cluster distances from diss_frame and assumes cross-cluster distances = sigma
  raoQ_before <- sigma +
    2 *
      sum(diss_block$Abundance1 * diss_block$Distance * diss_block$Abundance2) -
    sigma * sum(p^2) -
    2 * sigma * sum(diss_block$Abundance1 * diss_block$Abundance2)

  # Compute Rao's quadratic entropy after clustering
  raoQ_after <- sigma - sigma * sum(p_clust^2)

  # Check for edge case
  if (abs(sigma - raoQ_before) < .Machine$double.eps) {
    stop(
      "Cannot compute multiplicity: denominator (sigma - Q_before) is effectively zero"
    )
  }

  # Compute multiplicity
  (sigma - raoQ_after) / (sigma - raoQ_before)
}
