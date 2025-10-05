
# x <- 1:10
# names(x) <- LETTERS[1:10]
# diss <- dist(x)
# visualize_dist(diss)

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
  dat$sign <- ifelse(dat$z >= 0, "pos", "neg")

  # Breaks for legend similar to pretty() midpoints in ade4 legend
  # We show |value| directly; scale_size_area makes point area ~ value.
  breaks_mag <- pretty(dat$mag, n = 4)
  breaks_mag <- breaks_mag[breaks_mag > 0]
  if (length(breaks_mag) == 0) breaks_mag <- unique(dat$mag)

  # Plot
  p <- ggplot(dat, aes(x = col, y = row)) +
    geom_point(
      aes(size = mag, fill = sign),
      shape = 22,        # filled square
      color = "black",
      stroke = 0.3) +
    # Area-proportional sizing, max size scaled by csize
    scale_size_area(
      max_size = 10 * csize,
      breaks   = breaks_mag,
      name     = "|value|",
      guide    = if (clegend > 0) "legend" else "none") +
    scale_fill_manual(
      values = c(pos = "black", neg = "white"),
      name   = "sign",
      guide  = if (clegend > 0){ guide_legend(override.aes = list(shape = 22, size = 5)) } else { "none" } ) +
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

