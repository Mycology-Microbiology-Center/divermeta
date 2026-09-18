#' Rao quadratic entropy (Q) (Rao 1982)
#'
#' Computes Rao quadratic entropy \eqn{Q}{Q}, a diversity measure that accounts for
#' element abundances and their pairwise dissimilarities.
#'
#' @param ab Numeric vector of element abundances.
#' @param diss Numeric square matrix or `dist` object of pairwise dissimilarities among elements.
#'
#' @return Numeric scalar, Rao quadratic entropy \eqn{Q}{Q}.
#' @references
#' \itemize{
#' \item Rao CR (1982) Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21. \doi{10.1016/0040-5809(82)90004-1}.
#' }
#' @export
raoQuadratic <- function(ab, diss) {
  # Validate inputs
  if (any(is.na(ab)) || any(is.na(diss))) {
    stop("Input contains NA values\n")
  }
  if (any(ab < 0)) {
    stop("Abundances must be non-negative\n")
  }

  total_ab <- sum(ab, na.rm = TRUE)
  if (total_ab == 0) {
    stop("Total abundance cannot be zero")
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

  # If not already a dist object, coerce
  if (!inherits(diss, "dist")) {
    diss <- stats::as.dist(diss)
  }

  p <- as.vector(ab / total_ab)
  
  res <- dist_quadratic_form(p, diss)
  
  res
}


#' Rao quadratic entropy (Q) by blocks
#'
#' Computes Rao quadratic entropy \eqn{Q}{Q} from a compact three-column distance
#' table. Pairs listed in `diss_frame` use their (capped) distance; every pair not
#' listed is assumed to be at distance `sigma`.
#'
#' With \eqn{p_i} the relative abundances and \eqn{P} the set of unordered pairs
#' listed in `diss_frame`:
#' \deqn{Q = \sigma\left(1 - \sum_i p_i^2\right) - 2\sigma \sum_{(i,j) \in P} p_i p_j + 2 \sum_{(i,j) \in P} p_i p_j d_{ij}}{Q = sigma (1 - sum p_i^2) - 2 sigma sum_P p_i p_j + 2 sum_P p_i p_j d_ij}
#'
#' @param ids Vector of element identifiers (same order as `ab`).
#' @param ab Numeric vector of element abundances.
#' @param diss_frame Data frame with columns `ID1`, `ID2`, `Distance`, each unordered pair listed once.
#' @param sigma Numeric distance assumed for unlisted pairs; listed distances are capped at it.
#'
#' @return Numeric scalar, Rao quadratic entropy \eqn{Q}{Q}.
#' @keywords internal
#' @noRd
raoQuadratic.by_blocks <- function(ids, ab, diss_frame, sigma = 1) {
  if (
    !is.data.frame(diss_frame) ||
      !all(c("ID1", "ID2", "Distance") %in% colnames(diss_frame))
  ) {
    # Try to rename if columns exist but have different names
    if (!is.null(ncol(diss_frame)) && ncol(diss_frame) == 3) {
      diss_frame <- as.data.frame(diss_frame)
      colnames(diss_frame) <- c("ID1", "ID2", "Distance")
    } else {
      stop(
        "`diss_frame` must be a data frame with columns `ID1`, `ID2`, `Distance`"
      )
    }
  }

  # Normalize abundances
  p <- ab / sum(ab)
  names(p) <- ids

  # Cap distances at sigma
  diss_frame$Distance[diss_frame$Distance > sigma] <- sigma

  # Builds matrix multiplication by join
  diss_block <- diss_frame[diss_frame$ID1 != diss_frame$ID2, ]
  diss_block <- merge(
    diss_block,
    data.frame(ID1 = ids, Abundance1 = p),
    by = "ID1",
    all = FALSE
  )
  diss_block <- merge(
    diss_block,
    data.frame(ID2 = ids, Abundance2 = p),
    by = "ID2",
    all = FALSE
  )

  # Accounts for the listed distances and assumes unlisted distances = sigma
  sigma +
    2 *
      sum(diss_block$Abundance1 * diss_block$Distance * diss_block$Abundance2) -
    sigma * sum(p^2) -
    2 * sigma * sum(diss_block$Abundance1 * diss_block$Abundance2)
}
