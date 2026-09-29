#' Inventory multiplicity
#'
#' Computes inventory multiplicity \eqn{^{q}M}{M^q}, the within-cluster diversity component
#' under Hill-number partitioning, for every sample. It summarizes the average diversity inside
#' clusters given element abundances and their cluster memberships. Use `q`
#' to control abundance weighting (e.g., `q = 0` richness-like, `q = 1`
#' Shannon-type).
#'
#' @inheritParams raoQuadratic
#' @param clust Vector or factor of cluster memberships for each subunit. If named, names are
#'   subunit identifiers and every column of `abund` must be listed; otherwise it must follow the
#'   columns of `abund`. A data frame is also accepted, with two columns (subunit identifiers, units) or one
#'   column of units (named by its row names, if set).
#' @param q Numeric order of the Hill number (default `1`). Controls abundance weighting:
#'   \describe{
#'     \item{`q = 0`}{Richness-like weighting (all elements weighted equally)}
#'     \item{`q = 1`}{Shannon-type weighting (default, proportional to abundance)}
#'     \item{`q = 2`}{Simpson-type weighting (emphasizes abundant elements)}
#'   }
#'   Must be non-negative.
#'
#' @return Named numeric vector with the inventory multiplicity \eqn{^{q}M}{M^q} of every sample.
#'   This value can be interpreted as the effective number of equally abundant elements per cluster.
#'   A value of 1 indicates each cluster contains essentially one element (no diversity lost), while
#'   higher values indicate greater intra-cluster diversity. Samples with zero total abundance get `NA`.
#'
#' @details
#' Subunits with zero abundance in a sample are ignored for that sample.
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
#' # Three clusters with three subunits each (a, b, c)
#' clust <- c(a1 = "a", a2 = "a", a3 = "a", b1 = "b", b2 = "b", b3 = "b",
#'            c1 = "c", c2 = "c", c3 = "c")
#' abund <- rbind(
#'   High = rep(10, 9),                     # every subunit equally abundant
#'   Low = c(10, 0, 0, 20, 0, 0, 30, 0, 0), # one subunit per cluster
#'   Uneven = c(10, 5, 2, 8, 8, 8, 20, 1, 0)
#' )
#' colnames(abund) <- names(clust)
#'
#' multiplicity.inventory(abund, clust)          # High: 3, Low: 1
#' multiplicity.inventory(abund, clust, q = 0)  # Richness-like
#' multiplicity.inventory(abund, clust, q = 2)  # Simpson-type
#'
multiplicity.inventory <- function(abund, clust, q = 1) {
  # Input validation
  .check_q(q)
  ab <- .parse_abund(abund)
  cl <- .parse_clust(clust, ab$subunits)
  .multiplicity.inventory_core(ab, cl, q)
}


.multiplicity.inventory_core <- function(ab, cl, q) {
  # Abundances of every unit
  U <- ab$A %*% .unit_indicator(cl$unit_of, length(cl$unit_ids))

  # True diversity before and after clustering
  div_gamma <- .hill_numbers(ab$A, ab$N, q)
  div_beta <- .hill_numbers(U, ab$N, q)

  # Multiplicity
  .by_sample(div_gamma / div_beta, ab)
}


# Hill number of order q of every row of M (row totals N), over its nonzero entries
.hill_numbers <- function(M, N, q) {
  nz <- .nonzero(M)
  n_rows <- length(N)
  p <- nz$x / N[nz$i]
  rows <- factor(nz$i, levels = seq_len(n_rows))

  if (abs(q - 1) < .q_one_tol) {
    H <- -as.vector(tapply(p * log(p), rows, sum, default = 0))
    return(exp(H))
  }
  as.vector(tapply(p^q, rows, sum, default = 0))^(1 / (1 - q))
}
