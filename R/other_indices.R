#' Functional diversity (Chiu & Chao 2014)
#'
#' Computes the functional diversity, derived from Chiu & Chao (2014), 
#' assuming no intra-specific variation and using a default q value of 1.
#'
#' @param ab A numeric vector of element abundances.
#' @param diss A numeric matrix representing the dissimilarities or distances between elements.
#' @param sig A numeric value determining the threshold (sigma) at which two units are considered different.
#'
#' @return A numeric value representing the functional diversity, \eqn{FD_{\sigma}}.
#' @references
#' Chiu CH, Chao A (2014) Distance-based functional diversity measures and their
#' decomposition: A framework based on Hill numbers. PLOS ONE 9(7).
#' DOI:10.1371/journal.pone.0100014. URL: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014
#' @export
diversity.functional <- function(ab, diss, sig = 1) {
  diss[diss > sig] <- sig
  vals <- 1 / (1 - raoQuadratic(ab, diss) / sig)

  return(vals)
}



#' Functional diversity (Chiu & Chao 2014)
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
#' @references
#' Chiu CH, Chao A (2014) Distance-based functional diversity measures and their
#' decomposition: A framework based on Hill numbers. PLOS ONE 9(7).
#' DOI:10.1371/journal.pone.0100014. URL: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014
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



#' Redundancy index (Ricotta & Pavoine 2025)
#'
#' Computes functional redundancy as the degree of functional similarity among distinct species.
#'
#' Formally: \eqn{Re = D - Q}, where \eqn{D = 1 - \sum_i p_i^2} is Simpson's
#' index and \eqn{Q = \sum_{i,j} p_i \, p_j \, \delta_{ij}} is Rao's quadratic
#' entropy. This redundancy measure is mathematically fixed to \eqn{q = 2}
#' because both components (Simpson and Rao) belong to the Simpson–Rao family
#' of quadratic diversity indices.
#'
#' @param ab A numeric vector of element abundances.
#' @param diss A numeric matrix of pairwise dissimilarities among elements
#'   scaled in \[0, 1\]. Must be square with dimension equal to length of `ab`.
#'
#' @return A numeric value representing functional redundancy, \eqn{Re = D - Q}.
#' @references
#' \itemize{
#' \item Ricotta C, Pavoine S (2025) What do functional diversity, redundancy, rarity, and originality actually measure? A theoretical guide for ecologists and conservationists. Ecological Complexity 61. \doi{10.1016/j.ecocom.2025.101116}. \url{https://www.sciencedirect.com/science/article/pii/S1476945X25000017}
#' \item Rao CR (1982) Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology 21. \doi{10.1016/0040-5809(82)90004-1}.
#' }
#' @export
#' 
redundancy <- function(ab, diss) {

  # Validate inputs
  if (!is.numeric(ab)) {
    stop("Abundance vector must be numeric")
  }
  if (!is.matrix(diss)) {
    stop("Dissimilarity input must be a matrix")
  }
  if (length(ab) != nrow(diss) || nrow(diss) != ncol(diss)) {
    stop("Abundance length must match square dissimilarity matrix dimensions")
  }
  if (any(is.na(ab)) || any(is.na(diss))) {
    stop("Input contains NA values")
  }
  if (any(ab < 0)) {
    stop("Abundances must be non-negative")
  }

  total_ab <- sum(ab)
  if (total_ab == 0) {
    stop("Total abundance cannot be zero")
  }

  p <- as.vector(ab / total_ab)

  # Simpson's diversity (q = 2): D = 1 - sum p_i^2
  D <- 1 - sum(p^2)

  # Rao's quadratic entropy
  Q <- raoQuadratic(ab, diss)

  # Functional redundancy
  res <- D - Q
  return(res)
}

