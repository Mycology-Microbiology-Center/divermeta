#' Rao Quadratic
#'
#' TODO
#'
#' @param ab 1 x n vector of abundances.
#' @param diss n x n matrix with dissimilarities / distances
#'
#' @return Q.
#' @export
raoQuadratic <- function(ab, diss) {
  P <- as.matrix(ab / sum(ab))
  return(sum(P %*% diss %*% t(P)))
}
