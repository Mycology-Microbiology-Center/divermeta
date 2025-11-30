#' Quadratic Form for a Distance Object
#'
#' Computes the quadratic form \eqn{\sum_{i,j} 2 p_i p_j m_{ij}} 
#' for a vector \code{p} and a \code{dist} object \code{m}. This is equivalent to
#' \eqn{P^T M P}, where \code{M} is the full distance matrix, without
#' explicitly constructing the full matrix.
#'
#' @param p Numeric vector of length \( n \), representing the weights or abundances.
#' @param diss A `dist` object of size \( n \) (as returned by `vegdist` or `dist`),
#'   containing the pairwise distances between observations.
#'
#' @return A numeric scalar, the value of the quadratic form \eqn{P^T M P}.
#'
#' @details
#' The function uses the fact that a `dist` object stores only the lower triangle of
#' the distance matrix in column-major order. It removes zero entries and then constructs
#' all pairs of indices corresponding to this storage order and computes the sum efficiently.
#'
#' @examples
#' set.seed(123)
#' p <- runif(5)  # abundances
#' diss <- stats::dist(matrix(runif(5 * 3), ncol = 3)) # distance matrix
#' dist_quadratic_form(p, diss)
#'
#' @export
dist_quadratic_form <- function(p, diss) {
  # Get indices of non-zero probabilities
  non_zero <- which(p > 0)
  n <- length(non_zero)

  if (n < 2) {
    return(0)
  } # Need at least 2 non-zero entries

  # Only generate combinations for non-zero entries
  idx <- utils::combn(non_zero, 2)

  # Convert matrix indices to vector index for dist object
  # Using the formula: k = n*(i-1) - i*(i-1)/2 + j-i  for i < j
  n_total <- length(p)
  vec_indices <- convert_to_dist_indices(idx[1, ], idx[2, ], n_total)

  sum(2 * p[idx[1, ]] * p[idx[2, ]] * diss[vec_indices])
}


#' Convert Matrix Indices to Distance Vector Position
#'
#' Converts row and column indices from a full distance matrix to the corresponding
#' position in the compact vector representation used by R's `dist` objects.
#' This function maps `(i, j)` matrix coordinates to the linear index `k` in the
#' lower-triangular distance vector.
#'
#' @param i Integer vector of row indices (should be less than `j`)
#' @param j Integer vector of column indices (should be greater than `i`)
#' @param n Integer, the size of the square distance matrix (number of objects)
#'
#' @return Integer vector of the same length as `i` and `j`, containing the
#'   positions in the distance vector corresponding to the input matrix indices.
#'
#' @details
#' R's `dist` objects store distance matrices in a compact form by only storing
#' the lower triangle (excluding the diagonal) as a vector. The storage order follows:
#' \itemize{
#'   \item Row-major order of the lower triangle: (2,1), (3,1), (3,2), (4,1), (4,2), (4,3), ...
#'   \item For an n×n matrix, the vector has length n×(n-1)/2
#'   \item The formula used is: \eqn{k = (i-1) \times n - \frac{(i-1) \times i}{2} + (j-i)} for \eqn{i < j}
#' }
#'
#' @note
#' This function assumes that `i < j` for all input pairs. If `i > j`, the result
#' will be incorrect as it corresponds to the upper triangle, which is not stored
#' in `dist` objects. The function does not check for valid indices.
#' @export
convert_to_dist_indices <- function(i, j, n) {
  # formula: k = n*(i-1) - i*(i-1)/2 + j-i  for i < j
  (i - 1) * n - (i - 1) * i / 2 + (j - i)
}


#' Compute Distance Matrix Between Clusters
#'
#' Computes a distance matrix between clusters based on pairwise distances
#' between elements in different clusters. Supports various linkage methods
#' for aggregating distances between cluster members.
#'
#' @param diss A distance matrix or `dist` object containing pairwise
#'   distances between all elements
#' @param clust Integer or factor vector of cluster assignments for each element
#' @param clust_ids_order Integer or factor vector specifying the order of
#'   cluster IDs to use in the output distance matrix
#' @param method Character string specifying the linkage method for
#'   aggregating distances between clusters. One of: "average", "min",
#'   "max", or "sigma". Sigma sets all distances between clusters to be
#'    of value `sig`. Default is "sigma".
#' @param sig Numeric value specifying the maximum distance threshold
#'   when using the "sigma" method. Distances greater than sigma are
#'   capped at this value. Default is 1.
#'
#' @return A distance matrix or `dist` object (matching the input type)
#'   containing distances between clusters. The dimensions correspond to
#'   `length(clust_ids_order)`.
#'
#' @details
#' This function computes inter-cluster distances by aggregating all pairwise
#' distances between elements belonging to different clusters. The available
#' linkage methods are:
#' \itemize{
#'   \item \code{"average"}: Mean distance between all cross-cluster pairs
#'   \item \code{"min"}: Minimum distance between any cross-cluster pair
#'   \item \code{"max"}: Maximum distance between any cross-cluster pair
#'   \item \code{"sigma"}: Returns a matrix where all inter-cluster distances
#'     are set to \code{sigma} (useful for creating complete graphs)
#' }
#'
#' The function preserves the input type: if \code{diss} is a \code{dist} object,
#' the output will also be a \code{dist} object; if it's a full matrix, the
#' output will be a full matrix.
#'
#' @note
#' The \code{clust_ids_order} parameter determines the order of clusters in the
#' output distance matrix. This is useful when you want to maintain a specific
#' cluster ordering rather than using the natural order of cluster IDs found
#' in \code{clust}.
#'
#' When using the "sigma" method, all inter-cluster distances are set to the
#' specified sigma value, creating a complete graph where all clusters are
#' equally distant from each other.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' n <- 10
#' diss <- dist(matrix(rnorm(n * 3), ncol = 3))
#' clust <- c(1, 1, 2, 2, 3, 3, 4, 4, 1, 2)
#' clust_ids <- sort(unique(clust))
#'
#' # Compute cluster distances using different methods
#' clust_dist_avg <- cluster_distance_matrix(diss, clust, clust_ids, "average")
#' clust_dist_min <- cluster_distance_matrix(diss, clust, clust_ids, "min")
#' clust_dist_max <- cluster_distance_matrix(diss, clust, clust_ids, "max")
#' clust_dist_sigma <- cluster_distance_matrix(diss, clust, clust_ids, "sigma", sig = 2)
#'
#' # Compare the results
#' print(as.matrix(clust_dist_avg))
#' print(as.matrix(clust_dist_min))
#' print(as.matrix(clust_dist_max))
#' print(clust_dist_sigma)
#'
#' # Use with full matrix input
#' diss_matrix <- as.matrix(diss)
#' clust_dist_matrix <- cluster_distance_matrix(diss_matrix, clust, clust_ids, "average")
#'
#'
#' @export
cluster_distance_matrix <- function(
  diss,
  clust,
  clust_ids_order,
  method = "sigma",
  sig = 1
) {
  n_clust <- length(clust_ids_order)
  n_elem <- length(clust)

  is_dist <- FALSE
  if (inherits(diss, "dist")) {
    is_dist <- TRUE
  }

  dist_method <- NULL

  if (method == "average") {
    dist_method <- function(x) {
      mean(x)
    }
  } else if (method == "min") {
    dist_method <- function(x) {
      min(x)
    }
  } else if (method == "max") {
    dist_method <- function(x) {
      max(x)
    }
  } else if (method == "sigma") {
    diss_clust <- matrix(sig, nrow = n_clust, ncol = n_clust)
    diag(diss_clust) <- 0

    if (is_dist) {
      return(stats::as.dist(diss_clust))
    } else {
      return(diss_clust)
    }
  } else {
    stop(paste0(
      "Cluster distance method: ",
      method,
      " not supported"
    ))
  }

  diss_clust <- c()
  for (j in 2:n_clust) {
    for (i in 1:(j - 1)) {
      idx <- convert_to_dist_indices(i, j, n_clust)

      ci <- which(clust == clust_ids_order[i])
      cj <- which(clust == clust_ids_order[j])

      all_pairs <- expand.grid(ci = ci, cj = cj)
      indices <- data.frame(
        col = pmin(all_pairs$ci, all_pairs$cj),
        row = pmax(all_pairs$ci, all_pairs$cj)
      )

      indices <- unique(indices)

      dist_indices <- convert_to_dist_indices(indices$col, indices$row, n_elem)
      diss_clust[idx] <- dist_method(diss[dist_indices])
    }
  }

  diss_clust <- structure(
    diss_clust,
    class = "dist",
    Size = n_clust,
    Diag = FALSE,
    Upper = FALSE
  )

  if (!is_dist) {
    return(as.matrix(diss_clust))
  }

  diss_clust
}
