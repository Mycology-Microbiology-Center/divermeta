#' Rao's Quadratic Entropy
#'
#' Computes Rao's Quadratic Entropy, 
#' a measure of diversity that accounts for both element abundances and their dissimilarities.
#'
#' @param ab A numeric vector of element abundances. Each element represents the abundance of an entity.
#' @param diss A numeric matrix representing the dissimilarities or distances between entities. It should be a square matrix with dimensions equal to the length of `ab`.
#'
#' @return A numeric value representing Rao's Quadratic Entropy.
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

  P <- as.matrix(ab / total_ab)
  return(sum(P %*% diss %*% t(P)))
}
