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
