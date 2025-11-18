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
#' @param ab_clust Numeric vector of cluster abundances after clustering (sum of element
#'   abundances within each cluster).
#' @param diss_clust Numeric matrix or `dist` object of pairwise dissimilarities among clusters
#'   after clustering. Must be square and match the length of `ab_clust`.
#' @param sig Numeric cutoff \eqn{\sigma}{sigma} (default `1`) defining the threshold distance
#'   at which two units are considered maximally different. All distances greater than `sig`
#'   are capped at `sig`. This parameter should be chosen based on the biological meaning
#'   of distances in your dataset (e.g., maximum expected genetic distance, functional
#'   dissimilarity threshold).
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
#' @seealso [multiplicity.distance.by_blocks()] for computing from a compact distance table,
#'   [raoQuadratic()] for Rao's quadratic entropy,
#'   [diversity.functional()] for distance-based functional diversity
#'
#' @export
#'
#' @examples
#' # Example 1: High diversity within clusters
#' # Two clusters, each with functionally diverse elements
#' ab <- c(10, 10, 10, 10)  # 4 elements
#' clust <- c(1, 1, 2, 2)    # 2 clusters
#'
#' # Distance matrix: elements within clusters are dissimilar (0.3, 0.4)
#' # elements between clusters are maximally dissimilar (1.0)
#' diss <- matrix(c(
#'   0.0, 0.3, 1.0, 1.0,
#'   0.3, 0.0, 1.0, 1.0,
#'   1.0, 1.0, 0.0, 0.4,
#'   1.0, 1.0, 0.4, 0.0
#' ), nrow = 4, byrow = TRUE)
#'
#' # Clustered representation
#' ab_clust <- c(20, 20)  # Two clusters
#' diss_clust <- matrix(c(0, 1, 1, 0), nrow = 2)  # Clusters maximally different
#'
#' multiplicity.distance(ab, diss, ab_clust, diss_clust, sig = 1)
#'
#' # Example 2: Low diversity within clusters
#' # Elements within clusters are functionally similar (distance ~0)
#' diss_low <- matrix(c(
#'   0.0, 0.05, 1.0, 1.0,
#'   0.05, 0.0, 1.0, 1.0,
#'   1.0, 1.0, 0.0, 0.05,
#'   1.0, 1.0, 0.05, 0.0
#' ), nrow = 4, byrow = TRUE)
#'
#' multiplicity.distance(ab, diss_low, ab_clust, diss_clust, sig = 1)
#'
multiplicity.distance <- function(ab, diss, ab_clust, diss_clust, sig = 1) {

  # Input validation
  if (!is.numeric(ab) || !is.numeric(ab_clust)) {
    stop("Abundance vectors must be numeric")
  }
  if (!is.numeric(sig) || length(sig) != 1 || sig <= 0) {
    stop("`sig` must be a single positive numeric value")
  }
  if (sum(ab) == 0 || sum(ab_clust) == 0) {
    stop("Total abundance cannot be zero")
  }

  if (inherits(diss, "dist")) {
    diss <- as.matrix(diss)
  }
  if (inherits(diss_clust, "dist")) {
    diss_clust <- as.matrix(diss_clust)
  }

  # Cap distances at sigma (modify copies to avoid side effects)
  diss_clust[diss_clust > sig] <- sig
  diss[diss > sig] <- sig

  # Compute Rao's quadratic entropy
  Q_before <- raoQuadratic(ab, diss)
  Q_after <- raoQuadratic(ab_clust, diss_clust)

  # Check for edge case where denominator would be zero
  if (abs(sig - Q_before) < .Machine$double.eps) {
    stop("Cannot compute multiplicity: denominator (sigma - Q_before) is effectively zero")
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
multiplicity.distance.by_blocks <- function(ids, ab, diss_frame, clust, sigma = 1) {

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
  if (!is.data.frame(diss_frame) || !all(c("ID1", "ID2", "Distance") %in% colnames(diss_frame))) {
    # Try to rename if columns exist but have different names
    if (ncol(diss_frame) == 3) {
      colnames(diss_frame) <- c("ID1", "ID2", "Distance")
    } else {
      stop("`diss_frame` must be a data frame with columns `ID1`, `ID2`, `Distance`")
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
  diss_block <- merge(diss_block, data.frame(ID1 = ids, Abundance1 = p), by = "ID1", all = FALSE)
  diss_block <- merge(diss_block, data.frame(ID2 = ids, Abundance2 = p), by = "ID2", all = FALSE)

  # Compute Rao's quadratic entropy before clustering
  # This accounts for within-cluster distances from diss_frame and assumes cross-cluster distances = sigma
  raoQ_before <- sigma +
    2 * sum(diss_block$Abundance1 * diss_block$Distance * diss_block$Abundance2) -
    sigma * sum(p^2) -
    2 * sigma * sum(diss_block$Abundance1 * diss_block$Abundance2)

  # Compute Rao's quadratic entropy after clustering
  raoQ_after <- sigma - sigma * sum(p_clust^2)

  # Check for edge case
  if (abs(sigma - raoQ_before) < .Machine$double.eps) {
    stop("Cannot compute multiplicity: denominator (sigma - Q_before) is effectively zero")
  }

  # Compute multiplicity
  (sigma - raoQ_after) / (sigma - raoQ_before)
}
