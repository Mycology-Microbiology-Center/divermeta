# Average redundancy

Computes the average redundancy \\\overline{Re}\\ of every sample: the
average, across the units present in the sample, of the functional
redundancy \\Re\\ of each unit computed as in
[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)
over its own subunits only. Only the distances between subunits of the
same unit are used, so units that are not comparable with each other
(e.g. different enzyme families) can be summarized together, and no unit
dominates the index because of its size or abundance.

## Usage

``` r
average.redundancy(
  abund,
  diss,
  clust,
  normalize = FALSE,
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
  two subunits and the distance between them, or the path to a file with
  them (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  Distances between every pair of subunits of the same unit must be
  listed, in either orientation; rows for subunits of different units
  are ignored, also when their distance is invalid (NA, infinite or
  negative). Pairs with a subunit that has zero abundance in every
  sample may be missing, with a warning, except with
  `assume_homogeneous_abundance = TRUE` and
  `assume_max_reference_distance = FALSE`, where every subunit counts in
  the reference.

- clust:

  Vector or factor of cluster memberships: the unit of every subunit. If
  named, names are subunit identifiers; every column of `abund` must be
  listed, and listed subunits that are not in `abund` have zero
  abundance in every sample. Otherwise it must follow the columns of
  `abund`. A data frame is also accepted, with two columns (subunit
  identifiers, units) or one column of units (named by its row names, if
  set).

- normalize:

  Logical (default `FALSE`). If `TRUE`, the redundancy of every unit is
  divided by its Simpson's diversity \\S\\, giving \\1 - Q / S\\, before
  averaging (see
  [`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)).

- chunk_size:

  Number of rows of a distance file read at a time (default `1e6`).
  Lower values use less memory. Ignored when `diss` is in memory.

- check_large_distance_file:

  Logical (default `FALSE`). If `TRUE` and `diss` is a file, the file is
  fully checked with
  [`check_distance_file()`](https://mycology-microbiology-center.github.io/divermeta/reference/check_distance_file.md)
  before computing the index: every pair of subunits of the same unit
  must be listed exactly once, with a valid distance. This reads the
  file once more and uses temporary disk space. If `FALSE`, only the
  number of rows is checked and a warning is given. Ignored when `diss`
  is in memory, which is always fully checked. See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md).

- header:

  Logical (default `TRUE`). Whether a distance file starts with a header
  line (with several files, every one of them). With `TRUE`, a first
  line that looks like data (a number as third field) is an error; use
  `FALSE` for files without a header. Ignored when `diss` is in memory.
  See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md).

## Value

Named numeric vector with the average redundancy of every sample.
Samples with zero total abundance get 0.

## Details

Let \\Re(f)\\ be the functional redundancy of the composition \\f\\ of a
unit in the sample, computed from the relative abundances of its
subunits within the unit, and \\N_f\\ its total abundance. Then
\$\$\overline{Re} = \frac{1}{\left\|\\f : N_f \> 0\\\right\|} \sum\_{f
\\:\\ N_f \> 0} Re(f),\$\$ with \\Re(f) = S(f) - Q(f)\\, or \\Re(f) =
1 - Q(f) / S(f)\\ with `normalize = TRUE`.

A unit with a single present subunit has no redundancy: it contributes 0
to the average, also when normalized. Only units present in the sample
are averaged. As in
[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md),
distances greater than 1 are capped at 1, so the index is never
negative.

## See also

[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md),
[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md)

## Examples

``` r
# Two units: A with subunits a1, a2, a3 and B with subunits b1, b2
clust <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B")
diss <- data.frame(
  ID1 = c("a1", "a1", "a2", "b1"),
  ID2 = c("a2", "a3", "a3", "b2"),
  Distance = c(0.3, 0.8, 0.6, 0.4)
)
abund <- matrix(
  c(5, 3, 0, 2, 0,
    0, 1, 4, 2, 6,
    2, 0, 0, 0, 1),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("S1", "S2", "S3"), names(clust))
)

average.redundancy(abund, diss, clust)
#>        S1        S2        S3 
#> 0.1640625 0.1765000 0.0000000 
average.redundancy(abund, diss, clust, normalize = TRUE)
#>   S1   S2   S3 
#> 0.35 0.50 0.00 
```
