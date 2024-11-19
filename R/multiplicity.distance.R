#' Distance Multiplicity
#'
#' Distance Multiplicity defined as the ratio between the functional
#' diversities of the unclustered set and the cluster set
#'
#' @param ab 1 x n vector of abundances before clustering.
#' @param diss n x n matrix with dissimilarities / distances before clustering
#' @param ab_clust 1 x n vector of abundances after clustering.
#' @param diss_clust n x n matrix with dissimilarities / distances before
#' clustering
#' @param sig float value determining that two units are different
#'
#' @return \eqn{FM_{\sigma}}.
#' @export
multiplicity.distance <- function(ab, diss, ab_clust, diss_clust, sig = 1) {

  diss_clust[diss_clust > sig] = sig
  diss[diss > sig] = sig

  return((sig - raoQuadratic(ab_clust, diss_clust)) / (sig - raoQuadratic(ab, diss)))
}

#' Functional diversity
#'
#' Functional diversity derived from Chiu & Chao 2014, with q = 1 ans assuming
#' no intra specific variation. FD_sigma in the manuscript.
#'
#' @param ab 1 x n vector of abundances.
#' @param diss n x n matrix with dissimilarities / distances
#' @param sig float value determining that two units are different
#'
#' @return FD_sigma.
#' @export
diversity.functional <- function(ab, diss, sig = 1) {
  diss[diss > sig] <- sig
  vals <- 1 / (1 - raoQuadratic(ab, diss) / sig)

  return(vals)
}

#' Functional diversity (From Chiu & Chao 2014 )
#'
#' Functional diversity derived from Chiu & Chao 2014. This corresponds to
#' D(Q) in their paper and ^qFD in the manusctipt
#'
#' TODO: create the case for q = 1
#'
#' @param ab 1 x n vector of abundances.
#' @param diss n x n matrix with dissimilarities / distances
#' @param q q parameter for the Hill number
#'
#' @return ^qFD
#' @export
diversity.functional.traditional <- function(ab, diss, q = 1 + 10e-8) {
  P <- as.matrix(ab / sum(ab))
  Q <- raoQuadratic(ab, diss)
  Pq <- P^q

  vals <- sum(Pq %*% diss %*% t(Pq))

  vals <- vals / Q

  vals <- vals^(1 / (1 - q))

  return(sqrt(vals))
}




