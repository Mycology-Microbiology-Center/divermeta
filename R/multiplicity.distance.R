#' Distance-based multiplicity
#'
#' Computes distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} as the ratio of distance-based
#' functional diversity before vs after clustering, using a cutoff \eqn{\sigma}{sigma} to cap
#' pairwise distances. Provide abundances and dissimilarities both for the
#' unclustered elements and for the clustered representation.
#'
#' @param ab Numeric vector of element abundances before clustering.
#' @param diss Numeric matrix or `dist` of pairwise dissimilarities among elements before clustering.
#' @param ab_clust Numeric vector of element abundances after clustering.
#' @param diss_clust Numeric matrix or `dist` object of pairwise dissimilarities among clusters after clustering.
#' @param sig Numeric cutoff \eqn{\sigma}{sigma} at which two units are considered different (default `1`).
#'
#' @return Numeric scalar, the distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma}.
#' @seealso [raoQuadratic()], [diversity.functional()], [multiplicity.distance.by_blocks()]
#' @export
multiplicity.distance <- function(ab, diss, ab_clust, diss_clust, sig = 1) {
  diss_clust[diss_clust > sig] <- sig
  diss[diss > sig] <- sig

  (sig - raoQuadratic(ab_clust, diss_clust)) / (sig - raoQuadratic(ab, diss))
}



#' Distance-based multiplicity by blocks
#'
#' Computes distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} from a compact three-column
#' distance table. Distances are capped at \eqn{\sigma}{sigma} and distances between elements
#' from different clusters are assumed to be equal to \eqn{\sigma}{sigma}.
#'
#' @param ids Character or integer vector of element identifiers (length `n`).
#' @param ab Numeric vector of pre-clustering abundances (length `n`, same order as `ids`).
#' @param diss_frame Data frame with columns `ID1`, `ID2`, `Distance` for unique pairs of elements.
#' @param clust Vector or factor of cluster memberships for each element (length `n`).
#' @param sigma Numeric cutoff \eqn{\sigma}{sigma} at which two units are considered different (default `1`).
#'
#' @return Numeric scalar, the distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma}.
#' @seealso [multiplicity.distance()], [raoQuadratic()], [diversity.functional()]
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
