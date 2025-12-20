# Visualize an abundance matrix as a sized-tile plot

Creates a compact visualization of an abundance matrix where each cell
is drawn as a square whose area is proportional to the abundance value.
Rows are features (OTUs, species, genes, etc.) and columns are samples.

## Usage

``` r
visualize_abundance(
  abund,
  feature.labels = NA,
  sample.labels = NA,
  clabel.row = 1,
  clabel.col = 1,
  csize = 1,
  clegend = 1,
  grid = TRUE,
  transform = c("identity", "log1p", "sqrt")
)
```

## Arguments

- abund:

  A numeric matrix or `data.frame` with features in rows and samples in
  columns.

- feature.labels:

  Optional character vector of row labels (features). If `NA` (default),
  labels are taken from `rownames(abund)` when available, otherwise
  sequential indices are used.

- sample.labels:

  Optional character vector of column labels (samples). If `NA`
  (default), labels are taken from `colnames(abund)` when available,
  otherwise sequential indices are used.

- clabel.row:

  Numeric multiplier for the relative size of y-axis text.

- clabel.col:

  Numeric multiplier for the relative size of x-axis text.

- csize:

  Numeric multiplier controlling the maximum symbol (tile) size.

- clegend:

  If greater than 0, a size legend is shown; otherwise the legend is
  hidden.

- grid:

  Logical; whether to draw a light grid behind the tiles.

- transform:

  One of `"identity"`, `"log1p"`, or `"sqrt"` indicating a simple
  non-negative transform applied to the abundances for sizing.

## Value

A `ggplot2` object representing the visualization.

## Details

Values are visualized by magnitude; larger values produce larger
squares. The first row (first feature) is displayed at the top. When
`transform = "log1p"` or `"sqrt"`, non-negative guards are applied via
`pmax(0, x)` prior to the transform. Legend tick labels are shown on the
original (back-transformed) abundance scale for readability.

## Examples

``` r
set.seed(1)
abund <- matrix(  
  data = floor(rlnorm(n = 5 * 4, meanlog = 1, sdlog = 2.5)),
  nrow = 5, ncol = 4)
rownames(abund) <- paste0("feat", 1:5)
colnames(abund) <- paste0("S", 1:4)
visualize_abundance(abund)

```
