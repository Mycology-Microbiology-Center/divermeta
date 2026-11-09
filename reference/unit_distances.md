# Distances between units

Computes the distance between every pair of units (clusters) from the
pairwise distances between their subunits, given as a three-column
table. The result is itself a three-column table, which can be passed as
`diss_clust` to
[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
with `method = "custom"`.

## Usage

``` r
unit_distances(
  diss,
  clust,
  method = "average",
  sig = 1,
  chunk_size = 1e+06,
  check_large_distance_file = FALSE,
  header = TRUE
)
```

## Arguments

- diss:

  Data frame with three columns (taken by position): the identifiers of
  two subunits and the distance between them, or the path to a file with
  them (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  Every pair of subunits of different units must be listed (in either
  orientation). Ignored for `method = "sigma"`, where it can be `NULL`.

- clust:

  Named vector of cluster memberships: names are subunit identifiers and
  values are the units they belong to. A data frame is also accepted,
  with two columns (subunit identifiers, units) or one column of units
  with the subunit identifiers as row names.

- method:

  Character string specifying the linkage method for aggregating the
  distances between the subunits of two units. One of: "average", "min",
  "max", or "sigma". Sigma sets all distances between units to `sig`.
  Default is "average".

- sig:

  Numeric distance between units used by the "sigma" method. Default is
  1.

- chunk_size:

  Number of rows of a distance file read at a time (default `1e6`).
  Lower values use less memory. Ignored when `diss` is in memory.

- check_large_distance_file:

  Logical (default `FALSE`). If `TRUE` and `diss` is a file, the file is
  fully checked with
  [`check_distance_file()`](https://mycology-microbiology-center.github.io/divermeta/reference/check_distance_file.md)
  before computing the index: every needed pair must be listed exactly
  once, with a valid distance. This reads the file once more and uses
  temporary disk space. If `FALSE`, only the number of rows is checked
  and a warning is given. Ignored when `diss` is in memory, which is
  always fully checked. See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md).

- header:

  Logical (default `TRUE`). Whether a distance file starts with a header
  line (with several files, every one of them). With `TRUE`, a first
  line that looks like data (a number as third field) is an error; use
  `FALSE` for files without a header. Ignored when `diss` is in memory.
  See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md).

## Value

Data frame with columns `ID1`, `ID2` and `Distance`, one row per pair of
units.

## Details

The distances are read sequentially: every row of `diss` updates a
counter (minimum, maximum or sum) of the pair of units of its subunits,
so memory depends on the number of pairs of units, not on the number of
pairs of subunits. Distances are not capped.

Only the pairs of subunits of different units are used: rows of pairs
inside a unit are ignored, also when their distance is invalid (NA,
infinite or negative). With a distance file, only the pairs of subunits
of different units are checked by `check_large_distance_file`.

## See also

[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md),
[distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)

## Examples

``` r
clust <- c(a1 = "A", a2 = "A", b1 = "B", c1 = "C")
diss <- data.frame(
  ID1 = c("a1", "a1", "a1", "a2", "a2", "b1"),
  ID2 = c("a2", "b1", "c1", "b1", "c1", "c1"),
  Distance = c(0.1, 0.5, 0.9, 0.7, 0.8, 0.6)
)
unit_distances(diss, clust, method = "average")
#>   ID1 ID2 Distance
#> 1   A   B     0.60
#> 2   A   C     0.85
#> 3   B   C     0.60
unit_distances(diss, clust, method = "min")
#>   ID1 ID2 Distance
#> 1   A   B      0.5
#> 2   A   C      0.8
#> 3   B   C      0.6
```
