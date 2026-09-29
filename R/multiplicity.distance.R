#' Distance-based multiplicity
#'
#' Computes distance-based multiplicity \eqn{\delta M_{\sigma}}{delta M_sigma} of every sample as
#' the ratio of distance-based diversity before vs after clustering, using a cutoff
#' \eqn{\sigma}{sigma} to cap pairwise distances. This index quantifies the diversity lost when
#' clustering elements into operational units, incorporating pairwise dissimilarities (e.g.,
#' genetic, functional, or phylogenetic distances) among elements.
#'
#' @inheritParams multiplicity.inventory
#' @inheritParams raoQuadratic
#' @param diss Data frame with three columns (taken by position): the identifiers of two subunits
#'   and the dissimilarity between them, in either orientation, or the path to a file with them
#'   (see [distance-files]). With `method = "sigma"` only the pairs of subunits of the same
#'   cluster are needed (every one of them must be listed); with any other method every pair of
#'   subunits must be listed. Only the subunits present in at least one sample count: pairs with
#'   a subunit that has zero abundance in every sample may be missing, with a warning (see
#'   Details for the linkage methods). Rows with identifiers not in `abund` and self pairs are
#'   ignored, and so are, with `method = "sigma"`, the rows of subunits of different clusters,
#'   also when their distance is invalid (NA, infinite or negative).
#' @param method Character string specifying the linkage method for
#'   aggregating distances between clusters. One of: "average", "min",
#'   "max", "sigma" or "custom". Sigma sets all distances between clusters to be
#'    of value `sig`. Default is "sigma".
#' @param sig Numeric cutoff \eqn{\sigma}{sigma} (default `1`) defining the threshold distance
#'   at which two units are considered maximally different. All distances greater than `sig`
#'   are capped at `sig`. This parameter should be chosen based on the biological meaning
#'   of distances in your dataset (e.g., maximum expected genetic distance, functional
#'   dissimilarity threshold).
#' @param diss_clust Data frame with three columns (taken by position): the identifiers of two
#'   clusters and the dissimilarity between them. Every pair of clusters must be listed. Only used
#'   when `method` is "custom" (see [unit_distances()]). Default is `NULL`.
#'
#' @return Named numeric vector with the distance-based multiplicity
#'   \eqn{\delta M_{\sigma}}{delta M_sigma} of every sample. This value represents the effective
#'   number of functionally distinct elements per cluster. A value of 1 indicates minimal functional
#'   diversity within clusters, while higher values indicate greater functional diversity lost
#'   through clustering. Samples with zero total abundance get `NA`.
#'
#' @details
#' Distance-based multiplicity is the ratio of the distance-based functional diversities
#' \eqn{\delta D_{\sigma}}{delta D_sigma} (see [diversity.functional()]) of every sample before and
#' after clustering.
#'
#' The function supports several methods for computing inter-cluster distances:
#' \itemize{
#'   \item \code{"sigma"}: All inter-cluster distances are set to \code{sig}, creating maximally
#'     different clusters (default). Subunits of different clusters are then also at distance
#'     \code{sig}: their listed distances are ignored.
#'   \item \code{"average"}: Mean distance between all cross-cluster element pairs
#'   \item \code{"min"}: Minimum distance between any cross-cluster element pair
#'   \item \code{"max"}: Maximum distance between any cross-cluster element pair
#'   \item \code{"custom"}: Use pre-computed cluster distances provided via \code{diss_clust}
#' }
#'
#' Cluster distances are computed from the uncapped subunit distances (over every subunit in
#' `abund`) and then capped at `sig`.
#'
#' Subunits with zero abundance in every sample do not need distances. With the linkage methods
#' ("average", "min" and "max"), their distances are still used for the cluster distances when
#' they are listed, so the result can change depending on whether they are listed. With a
#' distance file, duplicated rows of their pairs are not detected, also with
#' `check_large_distance_file = TRUE`. A cluster whose subunits all have zero abundance in every
#' sample needs no distance to the other clusters.
#'
#' @seealso [unit_distances()] for computing cluster distances,
#'   [raoQuadratic()] for Rao's quadratic entropy,
#'   [diversity.functional()] for distance-based functional diversity
#'
#' @export
#'
#' @examples
#' # Two clusters of two subunits
#' clust <- c(e1 = 1, e2 = 1, e3 = 2, e4 = 2)
#' abund <- matrix(
#'   c(10, 10, 10, 10,
#'     10, 0, 5, 20),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(c("S1", "S2"), names(clust))
#' )
#'
#' # Example 1: Using sigma method (default): only pairs inside clusters are needed
#' diss_within <- data.frame(ID1 = c("e1", "e3"), ID2 = c("e2", "e4"), Distance = c(0.3, 0.4))
#' multiplicity.distance(abund, diss_within, clust, method = "sigma", sig = 1)
#'
#' # Example 2: Using average linkage: every pair is needed
#' diss <- data.frame(
#'   ID1 = c("e1", "e1", "e1", "e2", "e2", "e3"),
#'   ID2 = c("e2", "e3", "e4", "e3", "e4", "e4"),
#'   Distance = c(0.3, 0.9, 1.0, 0.8, 1.0, 0.4)
#' )
#' multiplicity.distance(abund, diss, clust, method = "average", sig = 1)
#'
#' # Example 3: Using custom cluster distances
#' diss_clust <- data.frame(ID1 = 1, ID2 = 2, Distance = 0.8)
#' multiplicity.distance(abund, diss, clust, method = "custom", diss_clust = diss_clust)
#'
multiplicity.distance <- function(
  abund,
  diss,
  clust,
  method = "sigma",
  sig = 1,
  diss_clust = NULL,
  chunk_size = 1e6,
  check_large_distance_file = FALSE,
  header = TRUE
) {
  .check_multiplicity_distance_args(method, sig, diss_clust)
  ab <- .parse_abund(abund)
  cl <- .parse_clust(clust, ab$subunits)
  acc <- .multiplicity.distance_acc(ab, cl, method, sig, diss_clust)
  src <- .diss_source(
    diss, ab$subunits, chunk_size, check_large_distance_file, acc$scope, cl,
    header = header, need = .present_subunits(ab)
  )
  .consume_pairs(src, list(acc), ab$subunits, cl)[[1]]
}


.check_multiplicity_distance_args <- function(method, sig, diss_clust) {
  .check_sig(sig)
  if (!is.character(method) || length(method) != 1 ||
    !(method %in% c("sigma", "min", "max", "average", "custom"))) {
    stop(paste0("Cluster distance method: ", method, " not supported"))
  }
  if (method == "custom" && is.null(diss_clust)) {
    stop("diss_clust cannot be NULL if method is 'custom'")
  }
}


.multiplicity.distance_core <- function(ab, src, cl, method, sig, diss_clust) {
  acc <- .multiplicity.distance_acc(ab, cl, method, sig, diss_clust)
  .consume_pairs(src, list(acc), ab$subunits, cl)[[1]]
}


# Accumulates Rao's quadratic entropy before clustering over the pairs read (and
# the distances between clusters for the linkage methods)
.multiplicity.distance_acc <- function(ab, cl, method, sig, diss_clust) {
  n_units <- length(cl$unit_ids)
  P <- .relative(ab)
  Q_before <- numeric(nrow(P))

  if (method == "sigma") {
    # Subunits of different clusters are at distance sigma, so only the pairs
    # inside clusters differ from sigma
    update <- function(pairs) {
      within <- cl$unit_of[pairs$i] == cl$unit_of[pairs$j]
      w <- pmin(pairs$d[within], sig) - sig
      Q_before <<- Q_before + 2 * .pair_sums(P, pairs$i[within], pairs$j[within], w)[, 1]
    }
  } else {
    if (method == "custom") {
      pairs_clust <- .parse_diss(diss_clust, cl$unit_ids, name = "diss_clust")
      .require_all_pairs(pairs_clust, cl$unit_ids, name = "diss_clust", what = "clusters")
    } else {
      counter <- .unit_distance_counter(cl$unit_of, n_units, method)
    }
    update <- function(pairs) {
      # Cluster distances from the uncapped distances
      if (method != "custom") {
        counter$update(pairs)
      }
      pairs$d <- pmin(pairs$d, sig)
      Q_before <<- Q_before + .rao_q(P, pairs)
    }
  }

  finish <- function() {
    P_clust <- P %*% .unit_indicator(cl$unit_of, n_units)
    if (method == "sigma") {
      Q_before <- sig * (1 - as.numeric(Matrix::rowSums(P^2))) + Q_before
      Q_after <- sig * (1 - as.numeric(Matrix::rowSums(P_clust^2)))
    } else {
      if (method != "custom") {
        pairs_clust <- counter$finish()
        # Units with no subunit present in any sample have zero weight, and
        # their distances may be missing
        present <- as.numeric(Matrix::colSums(ab$A)) > 0
        .require_unit_pairs(pairs_clust, cl$unit_ids, sort(unique(cl$unit_of[present])))
      }
      # Cap distances at sigma
      pairs_clust$d <- pmin(pairs_clust$d, sig)
      Q_after <- .rao_q(P_clust, pairs_clust)
    }

    # Samples where the denominator is effectively zero
    den <- sig - Q_before
    degenerate <- ab$N > 0 & abs(den) < .Machine$double.eps
    if (any(degenerate)) {
      warning(paste0(
        "Cannot compute multiplicity: denominator (sigma - Q_before) is effectively zero. ",
        "Returning NA for samples: ",
        paste(ab$samples[degenerate], collapse = ", ")
      ))
    }

    res <- (sig - Q_after) / den
    res[degenerate] <- NA_real_
    .by_sample(res, ab)
  }

  list(scope = if (method == "sigma") "within" else "all", update = update, finish = finish)
}
