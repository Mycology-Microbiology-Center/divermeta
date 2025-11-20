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
#' the distance matrix in column-major order. It removes zero entries and then constructs 
#' all pairs of indices corresponding to this storage order and computes the sum efficiently.
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
  # Get indices of non-zero probabilities
  non_zero <- which(p > 0)
  n <- length(non_zero)

  if (n < 2) return(0)  # Need at least 2 non-zero entries

  # Only generate combinations for non-zero entries
  idx <- combn(non_zero, 2)

  # Convert matrix indices to vector index for dist object
  # Using the formula: k = n*(i-1) - i*(i-1)/2 + j-i  for i < j
  n_total <- length(p)
  vec_indices <- (idx[1, ] - 1) * n_total - (idx[1, ] - 1) * idx[1, ] / 2 + (idx[2, ] - idx[1, ])

  sum(2 * p[idx[1, ]] * p[idx[2, ]] * diss[vec_indices])
}
