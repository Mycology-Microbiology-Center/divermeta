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
multiplicity.distance <- function(ab, diss, ab_clust, diss_clust, sig = 1) {
  diss_clust[diss_clust > sig] <- sig
  diss[diss > sig] <- sig

  (sig - raoQuadratic(ab_clust, diss_clust)) / (sig - raoQuadratic(ab, diss))
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
multiplicity.distance.by_blocks <- function(ids, ab, diss_frame, clust, sigma = 1) {
  # Checks
  if (length(ids) == 0) {
    return(0)
  }

  # Renames
  colnames(diss_frame) <- c("ID1", "ID2", "Distance")

  # Computes post Clustered
  ab_clust <- tapply(ab, clust, sum)
  p_clust <- ab_clust / sum(ab_clust)

  # Normalizes abundance
  p <- ab / sum(ab)

  # Cuts off at sigma
  diss_frame$Distance[diss_frame$Distance > sigma] <- sigma

  # Builds matrix multiplication by join
  diss_block <- diss_frame[diss_frame$ID1 != diss_frame$ID2, ]
  diss_block <- merge(diss_block, data.frame(ID1 = ids, Abundance1 = p), by = "ID1", all = FALSE)
  diss_block <- merge(diss_block, data.frame(ID2 = ids, Abundance2 = p), by = "ID2", all = FALSE)


  # Computes RaoQ before and after Clustering
  raoQ_before <- sigma + 2 * sum(diss_block$Abundance1 * diss_block$Distance * diss_block$Abundance2) - sigma * sum(p**2) - 2 * sigma * sum(diss_block$Abundance1 * diss_block$Abundance2)
  raoQ_after <- sigma - sigma * sum(p_clust**2)

  # Computes multiplicity
  (sigma - raoQ_after) / (sigma - raoQ_before)
}
