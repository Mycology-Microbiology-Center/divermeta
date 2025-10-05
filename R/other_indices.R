
#' Functional diversity
#'
#' Computes the functional diversity, derived from Chiu & Chao (2014), assuming no intra-specific variation and using a default q value of 1.
#'
#' @param ab A numeric vector of element abundances.
#' @param diss A numeric matrix representing the dissimilarities or distances between elements.
#' @param sig A numeric value determining the threshold (sigma) at which two units are considered different.
#'
#' @return A numeric value representing the functional diversity, \eqn{FD_{\sigma}}.
#' @export
diversity.functional <- function(ab, diss, sig = 1) {
  diss[diss > sig] <- sig
  vals <- 1 / (1 - raoQuadratic(ab, diss) / sig)

  return(vals)
}



#' Functional diversity (as per Chiu & Chao 2014 )
#'
#' Computes the functional diversity, derived from Chiu & Chao (2014). 
#' This corresponds to D(Q) in their paper and \eqn{^qFD} in the multiplicity manuscript.
#'
# TODO: create the case for q = 1
#' 
#' @param ab A numeric vector of element abundances.
#' @param diss A numeric matrix representing the dissimilarities or distances between elements.
#' @param q A numeric parameter for the Hill number, which determines the diversity order.
#'
#' @return A numeric value representing the functional diversity, \eqn{^qFD}.
#' @export
diversity.functional.traditional <- function(ab, diss, q = 1) {

  # Validate inputs
  if (!is.numeric(ab)) {
    stop("Abundance vector must be numeric")
  }
  if (!is.matrix(diss)) {
    stop("Dissimilarity input must be a matrix")
  }
  if (length(ab) != nrow(diss)) {
    stop("Abundance vector length must match dissimilarity matrix dimensions")
  }
  if (!is.numeric(q) || length(q) != 1) {
    stop("q parameter must be a single numeric value")
  }
  if (q < 0) {
    stop("q parameter must be positive")
  }

  if(q == 1)
    q = q + 10e-12

  P <- as.matrix(ab / sum(ab))
  Q <- raoQuadratic(ab, diss)
  Pq <- P^q

  vals <- sum(Pq %*% diss %*% t(Pq))

  vals <- vals / Q

  vals <- vals^(1 / (1 - q))

  return(sqrt(vals))
}

