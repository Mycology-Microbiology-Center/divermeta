#' Visualize an abundance matrix as a sized-tile plot
#'
#' Creates a compact visualization of an abundance matrix where each cell is
#' drawn as a square whose area is proportional to the abundance value.
#' Rows are features (OTUs, species, genes, etc.) and columns are samples.
#'
#' @param abund A numeric matrix or `data.frame` with features in rows and
#'   samples in columns.
#' @param feature.labels Optional character vector of row labels (features).
#'   If `NA` (default), labels are taken from `rownames(abund)` when available,
#'   otherwise sequential indices are used.
#' @param sample.labels Optional character vector of column labels (samples).
#'   If `NA` (default), labels are taken from `colnames(abund)` when available,
#'   otherwise sequential indices are used.
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
#' squares. The first row (first feature) is displayed at the top. When
#' `transform = "log1p"` or `"sqrt"`, non-negative guards are applied via
#' `pmax(0, x)` prior to the transform. Legend tick labels are shown on the
#' original (back-transformed) abundance scale for readability.
#'
#' @examples
#' set.seed(1)
#' abund <- matrix(  
#'   data = floor(rlnorm(n = 5 * 4, meanlog = 1, sdlog = 2.5)),
#'   nrow = 5, ncol = 4)
#' rownames(abund) <- paste0("feat", 1:5)
#' colnames(abund) <- paste0("S", 1:4)
#' visualize_abundance(abund)
#'
#' @export
#'
visualize_abundance <- function(
  abund,
  feature.labels = NA,
  sample.labels  = NA,
  clabel.row = 1,
  clabel.col = 1,
  csize      = 1,
  clegend    = 1,
  grid       = TRUE,
  transform  = c("identity", "log1p", "sqrt")
) {

  stopifnot(is.matrix(abund) || is.data.frame(abund))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for `visualize_abundance()`. Please install it with install.packages('ggplot2').", call. = FALSE)
  }

  # Coerce to data.frame for consistent handling
  if (is.matrix(abund)) {
    abund <- as.data.frame(abund, stringsAsFactors = FALSE)
  }

  # Default labels taken from dimnames (if NA)
  if (is.na(feature.labels) || is.na(sample.labels)) {
    feature.labels <- rownames(abund)
    sample.labels  <- colnames(abund)
  }

  # Default sequential labels (if NULL)
  if (is.null(feature.labels)) { feature.labels <- as.character(seq_len(nrow(abund))) }
  if (is.null(sample.labels))  { sample.labels  <- as.character(seq_len(ncol(abund))) }

  # Ensure numeric matrix
  if (!all(vapply(abund, is.numeric, logical(1)))) {
    stop("All columns of abund must be numeric.")
  }

  # Transform values for sizing
  tf <- match.arg(transform)
  raw_vals <- as.numeric(unlist(abund))
  if (tf == "log1p") {
    vals <- log1p(pmax(0, raw_vals))
  } else if (tf == "sqrt") {
    vals <- sqrt(pmax(0, raw_vals))
  } else {
    vals <- pmax(0, raw_vals)
  }

  # Long format
  dat <- data.frame(
    row = rep(feature.labels, each  = ncol(abund)),
    col = rep(sample.labels,  times = nrow(abund)),
    value = vals,
    stringsAsFactors = FALSE
  )

  # Factor levels to put first row on top
  dat$row <- factor(dat$row, levels = rev(feature.labels))
  dat$col <- factor(dat$col, levels = sample.labels)

  # Legend breaks: compute on original scale, display back-transformed labels
  orig_vals <- pmax(0, raw_vals)
  positive_orig <- orig_vals[orig_vals > 0]
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

  # Theme / grid
  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(clabel.col), angle = 90, vjust = 0.5, hjust = 0),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(clabel.row)),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(grid)) {
    p <- p + base_theme +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey85"))
  } else {
    p <- p + base_theme +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }

  return(p)
}


