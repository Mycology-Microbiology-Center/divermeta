
#' Visualize distances as a sized-tile plot
#'
#' Creates a compact visualization of the pairwise distances or dissimilarities between subunits,
#' given as a three-column table as for the indices. They are drawn as a symmetric matrix where
#' each cell is a square whose area is proportional to the distance.
#'
#' @param diss Data frame with three columns (taken by position): the identifiers of two subunits
#'   and the distance between them, in either orientation. Only in-memory tables are accepted.
#' @param ids Optional vector with the identifiers of the subunits to show, in the order they are
#'   drawn. If `NULL` (default), every identifier in the first two columns of `diss`, in order of
#'   first appearance. Rows with other identifiers are ignored.
#' @param labels Optional character vector with one label per subunit (in the order of `ids`).
#'   If `NULL` (default), the identifiers are used.
#' @param clabel.row Numeric multiplier for the relative size of y-axis text.
#' @param clabel.col Numeric multiplier for the relative size of x-axis text.
#' @param csize Numeric multiplier controlling the maximum symbol (tile) size.
#' @param clegend If greater than 0, a size legend is shown; otherwise the
#'   legend is hidden.
#' @param grid Logical; whether to draw a light grid behind the tiles.
#'
#' @return A `ggplot2` object representing the visualization.
#'
#' @details The distance of a subunit to itself is 0, and pairs not listed in `diss` are left
#' blank. Values are visualized by magnitude; larger values produce larger squares. The first
#' subunit is displayed at the top. `diss` is validated as for the indices (numeric, finite,
#' non-negative distances, and a pair listed again only with the same distance). Requires the
#' `ggplot2` package.
#'
#' @seealso \code{\link[ade4]{table.value}}

#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' diss <- data.frame(
#'   ID1 = c("a", "a", "a", "b", "b", "c"),
#'   ID2 = c("b", "c", "d", "c", "d", "d"),
#'   Distance = c(0.2, 0.7, 1.0, 0.6, 0.9, 0.4)
#' )
#' visualize_dist(diss)
#'
#' @export
#'
visualize_dist <- function(
  diss,
  ids        = NULL,
  labels     = NULL,
  clabel.row = 1,     # relative y-axis text size
  clabel.col = 1,     # relative x-axis text size
  csize      = 1,     # overall symbol size multiplier
  clegend    = 1,     # show legends if > 0
  grid       = TRUE   # show panel grid
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for `visualize_dist()`. Please install it with install.packages('ggplot2').", call. = FALSE)
  }
  if (.is_diss_file(diss)) {
    stop("`visualize_dist()` only accepts a distance table in memory, not a file")
  }
  if ((!is.data.frame(diss) && !is.matrix(diss)) || ncol(diss) != 3) {
    stop("`diss` must be a data frame with three columns: id_1, id_2, distance")
  }

  if (is.null(ids)) {
    ids <- unique(c(rbind(.as_ids(diss[, 1]), .as_ids(diss[, 2]))))
  }
  if (length(ids) == 0 || anyNA(ids)) {
    stop("`ids` must be a non-empty vector of subunit identifiers")
  }
  ids <- .as_ids(ids)
  if (anyDuplicated(ids) > 0) {
    stop("`ids` (subunit identifiers) must be unique")
  }
  labels <- .plot_labels(labels, ids, "labels", "subunits")
  pairs <- .parse_diss(diss, ids)

  # Symmetric matrix: zero on the diagonal, NA for the pairs not listed
  n <- length(ids)
  M <- matrix(NA_real_, n, n)
  diag(M) <- 0
  M[cbind(pairs$i, pairs$j)] <- pairs$d
  M[cbind(pairs$j, pairs$i)] <- pairs$d

  # Long format, in the (column-major) order of as.vector()
  dat <- data.frame(
    row = rep(labels, times = n),
    col = rep(labels, each = n),
    mag = as.vector(M),
    stringsAsFactors = FALSE
  )
  dat <- dat[!is.na(dat$mag), , drop = FALSE]

  # Factor levels to put the first subunit on top
  dat$row <- factor(dat$row, levels = rev(unique(labels)))
  dat$col <- factor(dat$col, levels = unique(labels))

  # Breaks for legend similar to pretty() midpoints in ade4 legend
  # scale_size_area makes point area ~ value.
  breaks_mag <- pretty(dat$mag, n = 4)
  breaks_mag <- breaks_mag[breaks_mag > 0]
  if (length(breaks_mag) == 0) breaks_mag <- unique(dat$mag)

  # Plot
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = col, y = row)) +
    ggplot2::geom_point(
      ggplot2::aes(size = mag),
      shape = 22,        # filled square
      fill = "black",
      color = "black",
      stroke = 0.3) +
    # Area-proportional sizing, max size scaled by csize
    ggplot2::scale_size_area(
      max_size = 10 * csize,
      breaks   = breaks_mag,
      name     = "Distance",
      guide    = if (clegend > 0) "legend" else "none") +
    ggplot2::coord_fixed() +
    # Keep every subunit on the axes, also those with no listed pair
    ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
    ggplot2::scale_y_discrete(position = "right", drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL)

  .tile_theme(p, clabel.row, clabel.col, grid)
}
