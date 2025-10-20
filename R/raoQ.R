#' Rao quadratic entropy (Q) (Rao 1982)
#'
#' Computes Rao quadratic entropy \eqn{Q}{Q}, a diversity measure that accounts for
#' element abundances and their pairwise dissimilarities.
#'
#' @param ab Numeric vector of element abundances.
#' @param diss Numeric square matrix of pairwise dissimilarities among elements.
#'
#' @return Numeric scalar, Rao quadratic entropy \eqn{Q}{Q}.
#' @references
#' \itemize{
#' \item Rao CR (1982) Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21. \doi{10.1016/0040-5809(82)90004-1}.
#' }
#' @export
raoQuadratic <- function(ab, diss) {

  # Validate inputs
  if (any(is.na(ab)) || any(is.na(diss))){
    stop("Input contains NA values\n")
  }
  if (any(ab < 0)){
    stop("Abundances must be non-negative\n")
  }

  total_ab <- sum(ab, na.rm = TRUE)
  if (total_ab == 0){
    stop("Total abundance cannot be zero")
  }

  P <- as.vector(ab / total_ab)


  return(sum(t(P) %*% diss %*% P))
}
