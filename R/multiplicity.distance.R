#' Distance based Multiplicity
#'
#' Distance Multiplicity defined as the ratio between the functional
#' diversities of the unclustered set and the cluster set
#'
#' @param ab 1 x n vector of abundances before clustering.
#' @param diss n x n matrix with dissimilarities / distances before clustering
#' @param ab_clust 1 x k vector of abundances after clustering.
#' @param diss_clust k x k matrix with dissimilarities / distances before
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


#' Distance based Multiplicity by Blocks
#'
#' Distance Multiplicity defined as the ratio between the functional
#' diversities of the unclustered set and the cluster set. Computed over a
#' block distance matrix. Here we assume that the distance between elements of
#' different clusters is maximum (i.e. sigma)
#'
#'
#' @param ids 1 x n vector of the ids of the elements.
#' @param ab 1 x n vector of abundances before clustering, in the same order
#' as the ids.
#' @param diss_frame a three columns data frame, with the distances between
#' elements. Columns are assumed to be in the order: ID1, ID2 and distance.
#' @param clust List or vector  element of size n with the unit's corresponding
#' cluster
#' @param sig float value determining that two units are different
#'
#' @return \eqn{FM_{\sigma}}.
#' @export
multiplicity.distance.by_blocks <- function(ids, ab, diss_frame, clust, sigma = 1)
{
  # Renames
  colnames(diss_frame) <- c("ID1","ID2","Distance")

  # Computes post Clustered
  ab_clust <- tapply(ab, clust, sum)
  p_clust <- ab_clust /sum(ab_clust)

  # Normalizes abundance
  p <- ab /sum(ab)

  # Cuts off at sigma
  diss_frame$Distance[diss_frame$Distance > sigma] <- sigma

  # Builds matrix multiplication by join
  diss_block <- diss_frame[diss_frame$ID1 != diss_frame$ID2,]
  diss_block <- merge(diss_block, data.frame(ID1 = ids, Abundance1 = p), by = "ID1", all = FALSE)
  diss_block <- merge(diss_block, data.frame(ID2 = ids, Abundance2 = p), by = "ID2", all = FALSE)


  # Computes RaoQ before and after Clustering
  raoQ_before <- sigma + 2*sum(diss_block$Abundance1 * diss_block$Distance * diss_block$Abundance2) - sigma*sum(p**2) - 2*sigma*sum(diss_block$Abundance1 * diss_block$Abundance2)
  raoQ_after <- sigma - sigma*sum(p_clust**2)

  # Computes multiplicity
  return((sigma - raoQ_after) / (sigma - raoQ_before))

}



#' Functional diversity
#'
#' Functional diversity derived from Chiu & Chao 2014, with q = 1 and assuming
#' no intra specific variation. FD_sigma in the manuscript. Function used for
#' testing.
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
#'
#' @param ab 1 x n vector of abundances.
#' @param diss n x n matrix with dissimilarities / distances
#' @param q q parameter for the Hill number
#'
#' @return ^qFD
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




