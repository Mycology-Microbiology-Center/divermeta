#' Visualize an abundance table as a sized-tile plot
#'
#' Creates a compact visualization of an abundance table where each cell is
#' drawn as a square whose area is proportional to the abundance value.
#' The plot follows the layout of the table: samples in rows and subunits
#' (OTUs, species, genes, etc.) in columns.
#'
#' @param abund Numeric matrix, data.frame or sparse `Matrix` of abundances with samples in rows
#'   and subunits in columns, as for the indices. Column names are the subunit identifiers.
#' @param sample.labels Optional character vector with one label per sample (row). If `NULL`
#'   (default), the row names of `abund` are used, or the row numbers when it has none.
#' @param subunit.labels Optional character vector with one label per subunit (column). If
#'   `NULL` (default), the column names of `abund` are used.
#' @param clabel.row Numeric multiplier for the relative size of y-axis text.
#' @param clabel.col Numeric multiplier for the relative size of x-axis text.
#' @param csize Numeric multiplier controlling the maximum symbol (tile) size.
#' @param clegend If greater than 0, a size legend is shown; otherwise the
#'   legend is hidden.
#' @param grid Logical; whether to draw a light grid behind the tiles.
#' @param transform One of `"identity"`, `"log1p"`, or `"sqrt"` indicating
#'   a simple non-negative transform applied to the abundances for sizing.
#'
#' @return A `ggplot2` object representing the visualization.
#'
#' @details Values are visualized by magnitude; larger values produce larger
#' squares. The first sample is displayed at the top. Legend tick labels are
#' shown on the original (back-transformed) abundance scale for readability.
#' Requires the `ggplot2` package.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' abund <- matrix(
#'   data = floor(rlnorm(n = 4 * 5, meanlog = 1, sdlog = 2.5)),
#'   nrow = 4, ncol = 5,
#'   dimnames = list(paste0("S", 1:4), paste0("f", 1:5))
#' )
#' visualize_abundance(abund)
#'
#' @export
#'
visualize_abundance <- function(
  abund,
  sample.labels  = NULL,
  subunit.labels = NULL,
  clabel.row = 1,
  clabel.col = 1,
  csize      = 1,
  clegend    = 1,
  grid       = TRUE,
  transform  = c("identity", "log1p", "sqrt")
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for `visualize_abundance()`. Please install it with install.packages('ggplot2').", call. = FALSE)
  }

  ab <- .parse_abund(abund)
  A <- as.matrix(ab$A)
  sample.labels <- .plot_labels(sample.labels, ab$samples, "sample.labels", "samples (rows of `abund`)")
  subunit.labels <- .plot_labels(subunit.labels, ab$subunits, "subunit.labels", "subunits (columns of `abund`)")

  # Transform values for sizing
  tf <- match.arg(transform)
  raw_vals <- as.vector(A)
  vals <- switch(tf,
    identity = raw_vals,
    log1p    = log1p(raw_vals),
    sqrt     = sqrt(raw_vals)
  )

  # Long format, in the (column-major) order of as.vector()
  dat <- data.frame(
    row = rep(sample.labels, times = ncol(A)),
    col = rep(subunit.labels, each = nrow(A)),
    value = vals,
    stringsAsFactors = FALSE
  )

  # Factor levels to put the first sample on top
  dat$row <- factor(dat$row, levels = rev(unique(sample.labels)))
  dat$col <- factor(dat$col, levels = unique(subunit.labels))

  # Legend breaks: compute on original scale, display back-transformed labels
  positive_orig <- raw_vals[raw_vals > 0]
  legend_orig_breaks <- pretty(positive_orig, n = 4)
  # If no positive values, keep empty to suppress legend (unless user forces it)
  if (length(legend_orig_breaks) == 0) {
    legend_trans_breaks <- numeric(0)
    legend_labels <- character(0)
  } else {
    forward_transform <- switch(tf,
      identity = function(x) x,
      log1p    = function(x) log1p(x),
      sqrt     = function(x) sqrt(x)
    )
    legend_trans_breaks <- forward_transform(legend_orig_breaks)
    legend_labels <- prettyNum(legend_orig_breaks, digits = 3, drop0trailing = TRUE)
  }

  # Plot
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = col, y = row)) +
    ggplot2::geom_point(
      ggplot2::aes(size = value),
      shape = 22,
      fill = "black",
      color = "black",
      stroke = 0.3
    ) +
    ggplot2::scale_size_area(
      max_size = 10 * csize,
      breaks   = legend_trans_breaks,
      labels   = legend_labels,
      name     = "Abundance",
      guide    = if (clegend > 0 && length(legend_trans_breaks) > 0) "legend" else "none"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::labs(x = NULL, y = NULL)

  .tile_theme(p, clabel.row, clabel.col, grid)
}


# Labels of the rows or columns of a plot: `labels` if given (one per element),
# `default` otherwise
.plot_labels <- function(labels, default, name, what) {
  if (is.null(labels)) {
    return(as.character(default))
  }
  if (length(labels) != length(default) || anyNA(labels)) {
    stop(paste0(
      "`", name, "` must have one label (not NA) per element: ", length(default), " ", what,
      ", ", length(labels), " labels"
    ))
  }
  as.character(labels)
}


# Theme and grid shared by the sized-tile plots
.tile_theme <- function(p, clabel.row, clabel.col, grid) {
  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(clabel.col), angle = 90, vjust = 0.5, hjust = 0),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(clabel.row)),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(grid)) {
    p + base_theme +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey85"))
  } else {
    p + base_theme +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }
}
