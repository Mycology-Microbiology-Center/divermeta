#' Inventory Multiplicity
#'
#' Computes the inventory multiplicity, derived from the factor of the true diversity of unclustered abundances and the clustered abundances.
#'
#' @param ab A numeric vector or list of element abundances.
#' @param clust A numeric vector or list indicating the cluster assignment for each element.
#' @param q A numeric parameter for the Hill number, which determines the diversity order.
#'
#' @return A numeric value representing the inventory multiplicity, \eqn{^qM}.
#' @export
#'
#' @examples
#' ab <- rep(10,9)
#' clust <- c(rep(1,3), rep(2,3), rep(3,3))
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
