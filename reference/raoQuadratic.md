# Rao quadratic entropy (Q) (Rao 1982)

Computes Rao quadratic entropy \\Q\\, a diversity measure that accounts
for element abundances and their pairwise dissimilarities, for every
sample.

## Usage

``` r
raoQuadratic(
  abund,
  diss,
  chunk_size = 1e+06,
  check_large_distance_file = FALSE,
  header = TRUE
)
```

## Arguments

- abund:

  Numeric matrix, data.frame or sparse `Matrix` of abundances with
  samples in rows and subunits in columns. Column names are the subunit
  identifiers.

- diss:

  Data frame with three columns (taken by position): the identifiers of
  two subunits and the dissimilarity between them, or the path to a file
  with them, read in chunks (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  Every pair of subunits present in at least one sample must be listed,
  in either orientation. Pairs with a subunit that has zero abundance in
  every sample are not needed: if some of them are missing, a warning is
  given. Rows with identifiers not in `abund` and self pairs are
  ignored.

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

Named numeric vector with Rao quadratic entropy \\Q\\ of every sample.
Samples with zero total abundance get 0.

## References

- Rao CR (1982) Diversity and dissimilarity coefficients: A unified
  approach. Theoretical Population Biology, 21.
  [doi:10.1016/0040-5809(82)90004-1](https://doi.org/10.1016/0040-5809%2882%2990004-1)
  .

## Examples

``` r
abund <- matrix(
  c(2, 1, 0,
    1, 1, 1),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("S1", "S2"), c("a", "b", "c"))
)
diss <- data.frame(
  ID1 = c("a", "a", "b"),
  ID2 = c("b", "c", "c"),
  Distance = c(0.4, 0.9, 0.6)
)
raoQuadratic(abund, diss)
#>        S1        S2 
#> 0.1777778 0.4222222 
```
