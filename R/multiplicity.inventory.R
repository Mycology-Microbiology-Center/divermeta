#' Inventory multiplicity
#'
#' Computes inventory multiplicity \eqn{^{q}M}{M^q}, the within-cluster diversity component
#' under Hill-number partitioning. It summarizes the average diversity inside
#' clusters given element abundances and their cluster memberships. Use `q`
#' to control abundance weighting (e.g., `q = 0` richness-like, `q = 1`
#' Shannon-type).
#'
#' @param ab Numeric vector of subunit abundances.
#' @param clust Vector or factor of cluster memberships for each element of `ab`.
#' @param q Numeric order of the Hill number (default `1`).
#'
#' @return Numeric scalar, the inventory multiplicity \eqn{^{q}M}{M^q} for the input.
#' @seealso [diversity.functional()], [multiplicity.distance()]
#' @export
#'
#' @examples
#' ab <- rep(10, 9)
#' clust <- c(rep(1, 3), rep(2, 3), rep(3, 3))
#' multiplicity.inventory(ab, clust)
#'
multiplicity.inventory <- function(ab, clust, q = 1) {

  # Removes any elements with zeros abundance
  zeros <- ab == 0
  if(sum(zeros) > 0){
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
