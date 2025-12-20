# Visualize a distance matrix as a sized-tile plot

Creates a compact visualization of a distance or dissimilarity matrix
where each cell is drawn as a square whose area is proportional to the
magnitude of the corresponding value. Labels are taken from the input
when available and can be customized.

## Usage

``` r
visualize_dist(
  diss,
  row.labels = NA,
  col.labels = NA,
  clabel.row = 1,
  clabel.col = 1,
  csize = 1,
  clegend = 1,
  grid = TRUE
)
```

## Arguments

- diss:

  A distance object (class `dist`), numeric matrix, or `data.frame`
  containing pairwise distances/dissimilarities.

- row.labels:

  Optional character vector of row labels. If `NA` (default), labels are
  taken from `diss` when available, otherwise sequential indices are
  used.

- col.labels:

  Optional character vector of column labels. If `NA` (default), labels
  are taken from `diss` when available, otherwise sequential indices are
  used.

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

## Value

A `ggplot2` object representing the visualization.

## Details

If `diss` is a `dist` object, labels are taken from the `Labels`
attribute when present. Values are visualized by absolute magnitude;
larger values produce larger squares. The first row is displayed at the
top.

## See also

[`table.value`](https://adeverse.github.io/ade4/reference/table.value.html)

## Examples

``` r
x <- 1:10
names(x) <- LETTERS[1:10]
diss <- stats::dist(x)
visualize_dist(diss)

```
