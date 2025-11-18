#' Inventory multiplicity
#'
#' Computes inventory multiplicity \eqn{^{q}M}{M^q}, the within-cluster diversity component
#' under Hill-number partitioning. It summarizes the average diversity inside
#' clusters given element abundances and their cluster memberships. Use `q`
#' to control abundance weighting (e.g., `q = 0` richness-like, `q = 1`
#' Shannon-type).
#'
#' @param ab Numeric vector of element (subunit) abundances. Elements with zero abundance
#'   are automatically removed before computation.
#' @param clust Vector or factor of cluster memberships for each element of `ab`.
#'   Must have the same length as `ab` (after removing zeros).
#' @param q Numeric order of the Hill number (default `1`). Controls abundance weighting:
#'   \describe{
#'     \item{`q = 0`}{Richness-like weighting (all elements weighted equally)}
#'     \item{`q = 1`}{Shannon-type weighting (default, proportional to abundance)}
#'     \item{`q = 2`}{Simpson-type weighting (emphasizes abundant elements)}
#'   }
#'   Must be non-negative.
#'
#' @return Numeric scalar, the inventory multiplicity \eqn{^{q}M}{M^q}. This value can be
#'   interpreted as the effective number of equally abundant elements per cluster. A value of 1
#'   indicates each cluster contains essentially one element (no diversity lost), while higher
#'   values indicate greater intra-cluster diversity.
#'
#' @references
#' \itemize{
#' \item Hill MO (1973) Diversity and evenness: a unifying notation and its consequences.
#'   Ecology 54(2):427-432. \doi{10.2307/1934352}
#' \item Jost L (2007) Partitioning diversity into independent alpha and beta components.
#'   Ecology 88(10):2427-2439. \doi{10.1890/06-1736.1}
#' }
#'
#' @seealso [multiplicity.distance()] for distance-based multiplicity,
#'   [diversity.functional()] for functional diversity indices
#'
#' @export
#'
#' @examples
#' ab <- rep(10, 9)
#' clust <- c(rep(1, 3), rep(2, 3), rep(3, 3))
#' multiplicity.inventory(ab, clust)
#'
multiplicity.inventory <- function(ab, clust, q = 1) {

  # Remove elements with zero abundance
  zeros <- ab == 0
  if (sum(zeros) > 0) {
    ab <- ab[!zeros]
    clust <- clust[!zeros]
  }

  # Checks
  if(length(ab) == 0 || length(clust) == 0)
    return(0)

  # Pre clust
  p <- ab /sum(ab)
  if(q == 1)
    div <- exp(-sum(p*log(p)))
  else
    div <- (sum(p^q))^(1/(1-q))

  # Post Clustering
  ab_clust <- tapply(ab, clust, sum)
  p_clust <- ab_clust /sum(ab_clust)

  if(q == 1)
    div_clust <- exp(-sum(p_clust*log(p_clust)))
  else
    div_clust <- (sum(p_clust^q))^(1/(1-q))

  return( div / div_clust)

}
