#' Quadratic Form for a Distance Object
#'
#' Computes the quadratic form \( \sum_{i,j} 2 \, p_i \, p_j \, m_{ij} \)
#' for a vector `p` and a `dist` object `m`. This is equivalent to
#' \( P^\top M P \) where `M` is the full distance matrix, without
#' explicitly constructing the full matrix.
#'
#' @param p Numeric vector of length \( n \), representing the weights or abundances.
#' @param diss A `dist` object of size \( n \) (as returned by `vegdist` or `dist`),
#'   containing the pairwise distances between observations.
#'
#' @return A numeric scalar, the value of the quadratic form \( P^\top M P \).
#'
#' @details
#' The function uses the fact that a `dist` object stores only the lower triangle of
#' the distance matrix in column-major order. It constructs all pairs of indices
#' corresponding to this storage order and computes the sum efficiently.
#'
#' @examples
#' library(vegan)
#' data(dune)
#' diss <- vegdist(dune)
#' p <- runif(nrow(dune))
#' dist_quadratic_form(p, diss)
#'
#' @export
dist_quadratic_form <- function(p, diss) {
    # Get row/column indices for lower-triangle elements in dist
    # This follows the order used in dist objects
    idx <- combn(length(p), 2)
    sum(2 * p[idx[1, ]] * p[idx[2, ]] * diss)
}
