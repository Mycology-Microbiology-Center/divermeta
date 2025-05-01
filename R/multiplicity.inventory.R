#' Inventory Multiplicity
#'
#' Inventory multiplicity directly from the factor of the true diversity
#' unclustered abundances and the clustered abundances.
#'
#' @param ab List or vector element of size n with the unit's abundances.
#' @param clust List or vector  element of size n with the unit's corresponding
#' cluster
#' @param q q parameter for the Hill number
#'
#' @return \eqn{^qM}
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





