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
#' # Example 1: High multiplicity (many diverse elements per cluster)
#' # Three clusters, each with 3 equally abundant elements
#' ab <- rep(10, 9)
#' clust <- c(rep(1, 3), rep(2, 3), rep(3, 3))
#' multiplicity.inventory(ab, clust)
#' # Result: 3 (each cluster has effective diversity of 3)
#'
#' # Example 2: Low multiplicity (one element per cluster)
#' # Three clusters, each with one element
#' ab <- c(10, 20, 30)
#' clust <- c(1, 2, 3)
#' multiplicity.inventory(ab, clust)
#' # Result: 1 (no diversity lost, each cluster is a single element)
#'
#' # Example 3: Unequal abundances within clusters
#' ab <- c(10, 5, 2,  # cluster 1: unequal abundances
#'         8, 8, 8,    # cluster 2: equal abundances
#'         20, 1)      # cluster 3: very unequal abundances
#' clust <- c(rep(1, 3), rep(2, 3), rep(3, 2))
#' multiplicity.inventory(ab, clust, q = 1)
#'
#' # Example 4: Using different q values
#' ab <- rep(1:12)
#' clust <- c(rep(1, 4), rep(2, 4), rep(3, 4))
#' multiplicity.inventory(ab, clust, q = 0)  # Richness-like
#' multiplicity.inventory(ab, clust, q = 1)  # Shannon-type (default)
#' multiplicity.inventory(ab, clust, q = 2)  # Simpson-type
#'
multiplicity.inventory <- function(ab, clust, q = 1) {

  # Input validation
  if (!is.numeric(ab)) {
    stop("`ab` must be a numeric vector")
  }
  if (!is.numeric(q) || length(q) != 1 || q < 0) {
    stop("`q` must be a single non-negative numeric value")
  }
  if (length(ab) != length(clust)) {
    stop("`ab` and `clust` must have the same length")
  }

  # Remove elements with zero abundance
  zeros <- ab == 0
  if (sum(zeros) > 0) {
    ab <- ab[!zeros]
    clust <- clust[!zeros]
  }

  # Check for empty input
  if (length(ab) == 0 || length(clust) == 0) {
    return(0)
  }

  # Check for single cluster case
  if (length(unique(clust)) == 1) {
    # Single cluster: multiplicity equals overall diversity
    p <- ab / sum(ab)
    if (q == 1) {
      # Handle log(0) case
      p_nonzero <- p[p > 0]
      if (length(p_nonzero) == 0) {
        return(0)
      }
      div <- exp(-sum(p_nonzero * log(p_nonzero)))
    } else {
      div <- (sum(p^q))^(1 / (1 - q))
    }
    return(div)
  }

  # Compute true diversity before clustering
  p <- ab / sum(ab)
  if (q == 1) {
    # Handle log(0) case
    p_nonzero <- p[p > 0]
    if (length(p_nonzero) == 0) {
      return(0)
    }
    div_gamma <- exp(-sum(p_nonzero * log(p_nonzero)))
  } else {
    div_gamma <- (sum(p^q))^(1 / (1 - q))
  }

  # Compute true diversity after clustering
  ab_clust <- tapply(ab, clust, sum)
  p_clust <- ab_clust / sum(ab_clust)

  if (q == 1) {
    p_clust_nonzero <- p_clust[p_clust > 0]
    if (length(p_clust_nonzero) == 0) {
      return(0)
    }
    div_beta <- exp(-sum(p_clust_nonzero * log(p_clust_nonzero)))
  } else {
    div_beta <- (sum(p_clust^q))^(1 / (1 - q))
  }

  # Multiplicity
  return(div_gamma / div_beta)
}
