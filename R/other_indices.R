#' Distance-based functional diversity (q = 1) (Chiu & Chao 2014)
#'
#' Computes distance-based functional diversity \eqn{\delta D_{\sigma}}{delta D_sigma} following Chiu & Chao (2014)
#' with Shannon-type weighting (order \eqn{q = 1}{q = 1}). Pairwise distances are capped at
#' the cutoff \eqn{\sigma}{sigma}.
#'
#' @param ab Numeric vector of element abundances.
#' @param diss Numeric matrix or `dist` object of pairwise dissimilarities among elements.
#' @param sig Numeric cutoff \eqn{\sigma}{sigma} at which two units are considered different (default `1`).
#'
#' @return Numeric scalar, the distance-based functional diversity \eqn{\delta D_{\sigma}}{delta D_sigma}.
#' @references
#' \itemize{
#' \item Chiu CH, Chao A (2014) Distance-based functional diversity measures and their decomposition: A framework based on Hill numbers. PLOS ONE 9(7). \doi{10.1371/journal.pone.0100014}. \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014}
#' }
#' @seealso [raoQuadratic()], [diversity.functional.traditional()]
#' @export
diversity.functional <- function(ab, diss, sig = 1) {
  diss[diss > sig] <- sig
  vals <- 1 / (1 - raoQuadratic(ab, diss) / sig)

  return(vals)
}



#' Distance-based functional diversity (order q) (Chiu & Chao 2014)
#'
#' Computes distance-based functional diversity \eqn{^{q}FD}{FD^q} following Chiu & Chao (2014).
#' This corresponds to D(Q) in their paper and \eqn{^{q}FD}{FD^q} in the divermeta manuscript.
#' For \eqn{q = 1}{q = 1}, an approximation is used internally.
#'
#' @param ab Numeric vector of element abundances.
#' @param diss Numeric matrix of pairwise dissimilarities among elements.
#' @param q Numeric order of the Hill number (non-negative).
#'
#' @return Numeric scalar, the distance-based functional diversity \eqn{^{q}FD}{FD^q}.
#' @references
#' \itemize{
#' \item Chiu CH, Chao A (2014) Distance-based functional diversity measures and their decomposition: A framework based on Hill numbers. PLOS ONE 9(7). \doi{10.1371/journal.pone.0100014}. \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014}
#' }
#' @seealso [diversity.functional()], [raoQuadratic()]
#' @export
diversity.functional.traditional <- function(ab, diss, q = 1) {
  # Validate inputs
  if (!is.numeric(ab)) {
    stop("Abundance vector must be numeric")
  }
  if (!is.numeric(q) || length(q) != 1) {
    stop("q parameter must be a single numeric value")
  }
  if (q < 0) {
    stop("q parameter must be positive")
  }

  # If dist object
  if (inherits(diss, "dist")) {
    n <- attr(diss, "Size")
    if (n != length(ab)) {
      stop(paste0(
        "Abundance vector and matrix must have compatible sizes. Matrix: ",
        n, "x", n, ". Vector: ", length(ab)
      ))
    }
  } else {
    dims <- dim(diss)
    if (dims[1] != dims[2]) {
      stop(paste0(
        "Distance matrix must be square. Matrix: ",
        dims[1], "x", dims[2], "."
      ))
    }

    if (dims[1] != length(ab)) {
      stop(paste0(
        "Abundance vector and matrix must have compatible sizes. Matrix: ",
        dims[1], "x", dims[2], ". Vector: ", length(ab)
      ))
    }
  }

  P <- as.vector(ab / sum(ab))
  Q <- raoQuadratic(ab, diss)

  # Handle q = 1 with the analytic limit q -> 1:
  # exp(-sum_ij d_ij p_i p_j ln(p_i) / Q) (square root already included)
  if (abs(q - 1) < .Machine$double.eps) {
    w <- ifelse(P > 0, P * log(P), 0)
    D <- if (inherits(diss, "dist")) as.matrix(diss) else diss
    return(exp(-as.numeric(t(w) %*% D %*% P) / Q))
  }

  Pq <- P^q

  if (inherits(diss, "dist")) {
    vals <- dist_quadratic_form(Pq, diss)
  } else {
    vals <- as.numeric(t(Pq) %*% diss %*% Pq)
  }

  vals <- vals / Q

  vals <- vals^(1 / (1 - q))

  return(sqrt(vals))
}



#' Functional redundancy (Re) (Ricotta & Pavoine 2025)
#'
#' Computes functional redundancy \eqn{Re}{Re}, a measure of the degree to which
#' distinct elements are functionally similar given their abundances and
#' pairwise dissimilarities. This implementation follows the Simpson–Rao
#' family and corresponds to the `q = 2` case.
#'
#' @param ab Numeric vector of element abundances.
#' @param diss Numeric square matrix or `dist` object of pairwise dissimilarities among elements
#'   scaled to the range \[0, 1\].
#'
#' @return Numeric scalar, functional redundancy `Re`.
#' @references
#' \itemize{
#' \item Ricotta C, Pavoine S (2025) What do functional diversity, redundancy, rarity, and originality actually measure? A theoretical guide for ecologists and conservationists. Ecological Complexity 61. \doi{10.1016/j.ecocom.2025.101116}. \url{https://www.sciencedirect.com/science/article/pii/S1476945X25000017}
#' \item Rao CR (1982) Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology 21. \doi{10.1016/0040-5809(82)90004-1}.
#' }
#' @seealso [raoQuadratic()], [redundancy.by_blocks()] for computing from a compact distance table
#' @export
#'
redundancy <- function(ab, diss) {


  # Validate inputs
  if (!is.numeric(ab)) {
    stop("Abundance vector must be numeric")
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

  # Edge case
  if (length(ab) == 1) {
    return(0)
  }

  # If dist object
  if (inherits(diss, "dist")) {
    n <- attr(diss, "Size")
    if (n != length(ab)) {
      stop(paste0(
        "Abundance vector and matrix must have compatible sizes. Matrix: ",
        n, "x", n, ". Vector: ", length(ab)
      ))
    }
  } else {
    dims <- dim(diss)
    if (dims[1] != dims[2]) {
      stop(paste0(
        "Distance matrix must be square. Matrix: ",
        dims[1], "x", dims[2], "."
      ))
    }

    if (dims[1] != length(ab)) {
      stop(paste0(
        "Abundance vector and matrix must have compatible sizes. Matrix: ",
        dims[1], "x", dims[2], ". Vector: ", length(ab)
      ))
    }
  }

  p <- as.vector(ab / total_ab)

  # Simpson's diversity (q = 2): D = 1 - sum p_i^2
  d <- 1 - sum(p^2)

  # Rao's quadratic entropy
  q <- raoQuadratic(ab, diss)

  # Functional redundancy
  res <- d - q

  res
}


#' Functional redundancy (Re) by blocks
#'
#' Computes functional redundancy \eqn{Re}{Re} from a compact three-column distance
#' table. This is an efficient implementation of [redundancy()] for large datasets where
#' storing the full distance matrix would be memory-intensive. Only the distances between
#' elements inside the same unit (block) need to be provided; every pair of elements not
#' listed in the table is assumed to be maximally different (distance equal to 1).
#'
#' @param ids Character or integer vector of element identifiers (length `n`). Must match
#'   the identifiers used in `diss_frame`.
#' @param ab Numeric vector of element abundances (length `n`, same order as `ids`).
#' @param diss_frame Data frame with columns `ID1`, `ID2`, `Distance` containing pairwise
#'   dissimilarities, scaled to the range \[0, 1\], for unique pairs of elements. Each
#'   unordered pair must be listed only once. Should only include within-unit pairs or
#'   pairs where distances are less than 1; unlisted pairs are automatically set to 1 and
#'   distances greater than 1 are capped at 1.
#'
#' @return Numeric scalar, functional redundancy `Re`.
#'
#' @details
#' Let \eqn{p_i}{p_i} be the relative abundance of element \eqn{i}{i} and \eqn{P}{P} the set
#' of unordered pairs \eqn{(i, j)}{(i, j)}, \eqn{i \neq j}{i != j}, listed in `diss_frame`.
#' Assuming every pair not in \eqn{P}{P} is at distance 1, Rao's quadratic entropy is
#' \deqn{Q = \left(1 - \sum_i p_i^2\right) - 2 \sum_{(i,j) \in P} p_i p_j + 2 \sum_{(i,j) \in P} p_i p_j d_{ij},}{Q = (1 - sum_i p_i^2) - 2 sum_P p_i p_j + 2 sum_P p_i p_j d_ij,}
#' so that the redundancy \eqn{Re = \left(1 - \sum_i p_i^2\right) - Q}{Re = (1 - sum_i p_i^2) - Q} reduces to
#' \deqn{Re = 2 \sum_{(i,j) \in P} p_i p_j \left(1 - d_{ij}\right).}{Re = 2 sum_P p_i p_j (1 - d_ij).}
#' Only the listed (within-unit) pairs contribute, so the computation scales with the
#' number of rows of `diss_frame` instead of \eqn{n^2}{n^2}.
#'
#' @references
#' \itemize{
#' \item Ricotta C, Pavoine S (2025) What do functional diversity, redundancy, rarity, and originality actually measure? A theoretical guide for ecologists and conservationists. Ecological Complexity 61. \doi{10.1016/j.ecocom.2025.101116}. \url{https://www.sciencedirect.com/science/article/pii/S1476945X25000017}
#' \item Rao CR (1982) Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology 21. \doi{10.1016/0040-5809(82)90004-1}.
#' }
#' @seealso [redundancy()] for the standard implementation using full matrices,
#'   [multiplicity.distance.by_blocks()] for distance-based multiplicity by blocks,
#'   [raoQuadratic()] for Rao's quadratic entropy
#' @export
#'
#' @examples
#' # Example: Compute redundancy from a distance table
#' ids <- c("elem1", "elem2", "elem3", "elem4")
#' ab <- c(10, 15, 20, 25)
#'
#' # Distance table: only within-unit pairs (other pairs assumed = 1)
#' diss_frame <- data.frame(
#'   ID1 = c("elem1", "elem3"),
#'   ID2 = c("elem2", "elem4"),
#'   Distance = c(0.3, 0.4),
#'   stringsAsFactors = FALSE
#' )
#'
#' redundancy.by_blocks(ids, ab, diss_frame)
#'
redundancy.by_blocks <- function(ids, ab, diss_frame) {
  # Validate inputs
  if (length(ids) == 0) {
    return(0)
  }
  if (!is.numeric(ab)) {
    stop("Abundance vector must be numeric")
  }
  if (length(ids) != length(ab)) {
    stop("`ids` and `ab` must have the same length")
  }
  if (any(is.na(ab))) {
    stop("Input contains NA values")
  }
  if (any(ab < 0)) {
    stop("Abundances must be non-negative")
  }

  total_ab <- sum(ab)
  if (total_ab == 0) {
    stop("Total abundance cannot be zero")
  }

  # Edge case
  if (length(ab) == 1) {
    return(0)
  }

  p <- as.vector(ab / total_ab)

  # Simpson's diversity (q = 2): D = 1 - sum p_i^2
  d <- 1 - sum(p^2)

  # Rao's quadratic entropy (unlisted pairs are at distance 1)
  q <- raoQuadratic.by_blocks(ids, ab, diss_frame, sigma = 1)

  # Functional redundancy
  res <- d - q

  res
}


#' Metagenomic Alpha-Diversity Index (MAD) (Finn 2024)
#'
#' Computes the Metagenomic Alpha-Diversity Index (MAD), a metric that measures the
#' average dissimilarity of elements (e.g., protein-encoding genes) within clusters
#' relative to cluster representatives. Unlike multiplicity, MAD does not account for
#' element abundances and decreases as the number of elements per cluster increases.
#'
#' @param clust Vector or factor of cluster memberships for each element (e.g., gene).
#'   Must have the same length as the number of rows/columns in `diss`.
#' @param diss Numeric square matrix of pairwise dissimilarities among elements.
#'   Should be scaled to the range \[0, 1\], where 0 indicates identical elements
#'   and 1 indicates maximally different elements.
#' @param representatives Optional named vector mapping cluster names to the index
#'   (position) of the representative element for each cluster. If `NULL` (default),
#'   the first element in each cluster is used as the representative.
#'
#' @return Numeric scalar, the Metagenomic Alpha-Diversity Index (MAD). Higher values
#'   indicate greater average dissimilarity within clusters.
#'
#' @details
#' Note: Unlike multiplicity indices, MAD does not incorporate element abundances and
#' decreases as cluster size increases, which may not reflect biological complexity.
#'
#' @references
#' \itemize{
#' \item Finn DR (2024) A metagenomic alpha-diversity index for microbial functional
#'   biodiversity. FEMS Microbiology Ecology 100(3):fiae019.
#'   \doi{10.1093/femsec/fiae019}
#' }
#'
#' @seealso [multiplicity.distance()] for abundance-weighted distance-based multiplicity,
#'   [multiplicity.inventory()] for inventory multiplicity
#'
#' @export
#'
#' @examples
#' # Example: Compute MAD for gene clusters
#' clust <- c(1, 1, 2, 2, 2, 3)
#' diss <- matrix(c(
#'   0.0, 0.2, 0.8, 0.9, 0.9, 0.9,
#'   0.2, 0.0, 0.8, 0.9, 0.9, 0.9,
#'   0.8, 0.8, 0.0, 0.1, 0.15, 0.9,
#'   0.9, 0.9, 0.1, 0.0, 0.12, 0.9,
#'   0.9, 0.9, 0.15, 0.12, 0.0, 0.9,
#'   0.9, 0.9, 0.9, 0.9, 0.9, 0.0
#' ), nrow = 6, byrow = TRUE)
#'
#' # Use default (first element as representative)
#' metagenomic.alpha.index(clust, diss)
#'
#' # Specify representatives explicitly
#' reps <- c("1" = 1, "2" = 3, "3" = 6)
#' metagenomic.alpha.index(clust, diss, representatives = reps)
#'
metagenomic.alpha.index <- function(clust, diss, representatives = NULL) {
  if (!is.matrix(diss)) {
    stop("Dissimilarity input must be a matrix")
  }
  if (length(clust) != nrow(diss) || nrow(diss) != ncol(diss)) {
    stop("Cluster memberships length must match square dissimilarity matrix dimensions")
  }
  if (any(is.na(clust)) || any(is.na(diss))) {
    stop("Input contains NA values")
  }

  # Cluster names
  cluster_names <- unique(clust)

  # Validate representatives
  if (!is.null(representatives)) {
    if (!is.numeric(representatives)) {
      stop("`representatives` must be a numeric vector (indices)")
    }
    if (any(representatives < 1) || any(representatives > length(clust))) {
      stop(paste0(
        "Representative indices must be between 1 and ",
        length(clust),
        "."
      ))
    }
    if (!all(cluster_names %in% names(representatives))) {
      missing <- cluster_names[!(cluster_names %in% names(representatives))]
      stop(paste0(
        "Missing representatives for clusters: ",
        paste(missing, collapse = ", ")
      ))
    }
  }

  # Summand function
  summand <- function(cluster_name) {
    ids <- which(clust == cluster_name)
    rep <- ids[1]
    if (!is.null(representatives)) {
      rep <- representatives[[cluster_name]]
    }

    1 + mean(diss[rep, ids])
  }

  # Sum and divide
  sum(sapply(cluster_names, summand)) / length(clust)
}
