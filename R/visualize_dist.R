
#' Visualize a distance matrix as a sized-tile plot
#'
#' Creates a compact visualization of a distance or dissimilarity matrix where
#' each cell is drawn as a square whose area is proportional to the magnitude
#' of the corresponding value. Labels are taken from the input when available
#' and can be customized.
#'
#' @param diss A distance object (class `dist`), numeric matrix, or
#'   `data.frame` containing pairwise distances/dissimilarities.
#' @param row.labels Optional character vector of row labels. If `NA` (default),
#'   labels are taken from `diss` when available, otherwise sequential indices
#'   are used.
#' @param col.labels Optional character vector of column labels. 
#'   If `NA` (default), labels are taken from `diss` when available, otherwise
#'   sequential indices are used.
#' @param clabel.row Numeric multiplier for the relative size of y-axis text.
#' @param clabel.col Numeric multiplier for the relative size of x-axis text.
#' @param csize Numeric multiplier controlling the maximum symbol (tile) size.
#' @param clegend If greater than 0, a size legend is shown; otherwise the
#'   legend is hidden.
#' @param grid Logical; whether to draw a light grid behind the tiles.
#'
#' @return A `ggplot2` object representing the visualization.
#'
#' @details If `diss` is a `dist` object, labels are taken from the `Labels`
#' attribute when present. Values are visualized by absolute magnitude; larger
#' values produce larger squares. The first row is displayed at the top.
#'
#' @seealso \code{\link[ade4]{table.value}}

#' @examples
#' x <- 1:10
#' names(x) <- LETTERS[1:10]
#' diss <- stats::dist(x)
#' visualize_dist(diss)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_size_area coord_fixed
#' @importFrom ggplot2 scale_x_discrete scale_y_discrete labs theme_minimal
#' @importFrom ggplot2 theme element_text rel element_blank element_line
#' 
visualize_dist <- function(
  diss,
  row.labels = NA,
  col.labels = NA,
  clabel.row = 1,     # relative y-axis text size
  clabel.col = 1,     # relative x-axis text size
  csize      = 1,     # overall symbol size multiplier
  clegend    = 1,     # show legends if > 0
  grid       = TRUE   # show panel grid
) {


  stopifnot("dist" %in% class(diss) || is.matrix(diss) || is.data.frame(diss))

  # Default labels taken from diss (if NA)
  if (is.na(row.labels) || is.na(col.labels)){
    if ("dist" %in% class(diss)) {
      if (! is.null(attr(diss, "Labels")) ) {
          row.labels <- attr(diss, "Labels")
          col.labels <- attr(diss, "Labels")
        } else {
          row.labels <- NULL   # will be set later
          col.labels <- NULL
        }
    } else {
        row.labels <- rownames(diss)
        col.labels <- colnames(diss)
    }
  }

  # Convert to data.frame
  if (is.matrix(diss)) {
    diss <- as.data.frame(diss, stringsAsFactors = FALSE)
  } else if ("dist" %in% class(diss)) {
    diss <- as.data.frame( as.matrix(diss) )
  }

  # TO DO - check the order of rows/cols

  # Default sequential labels (if NULL)
  if (is.null(row.labels)) { row.labels <- as.character(seq_len(nrow(diss))) }
  if (is.null(col.labels)) { col.labels <- as.character(seq_len(ncol(diss))) }

  # Ensure numeric matrix
  if (!all(vapply(diss, is.numeric, logical(1)))) {
    stop("All columns of diss must be numeric.")
  }

  # Reshape to long format
  dat <- data.frame(
    row = rep(row.labels, each  = ncol(diss)),
    col = rep(col.labels, times = nrow(diss)),
    z   = as.numeric(unlist(diss)),
    stringsAsFactors = FALSE)

  # Factor levels to put first row on top
  dat$row <- factor(dat$row, levels = rev(row.labels))
  dat$col <- factor(dat$col, levels = col.labels)

  # Encodings
  dat$mag  <- pmax(0, abs(dat$z))      # magnitude for area

  # Breaks for legend similar to pretty() midpoints in ade4 legend
  # We show |value| directly; scale_size_area makes point area ~ value.
  breaks_mag <- pretty(dat$mag, n = 4)
  breaks_mag <- breaks_mag[breaks_mag > 0]
  if (length(breaks_mag) == 0) breaks_mag <- unique(dat$mag)

  # Plot
  p <- ggplot(dat, aes(x = col, y = row)) +
    geom_point(
      aes(size = mag),
      shape = 22,        # filled square
      fill = "black",
      color = "black",
      stroke = 0.3) +
    # Area-proportional sizing, max size scaled by csize
    scale_size_area(
      max_size = 10 * csize,
      breaks   = breaks_mag,
      name     = "Distance",
      guide    = if (clegend > 0) "legend" else "none") +
    coord_fixed() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL)

  # Theme / grid
  base_theme <- theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(size = rel(clabel.col), angle = 90, vjust = 0.5, hjust = 0),
      axis.text.y = element_text(size = rel(clabel.row)),
      panel.grid.minor = element_blank() )

  if (isTRUE(grid)) {
    p <- p + base_theme +
      theme(panel.grid.major = element_line(color = "grey85"))
  } else {
    p <- p + base_theme +
      theme(panel.grid.major = element_blank())
  }

  return(p)
}

