# Relative multiplicity from reference diversities

Computes relative multiplicity \\RM\\ of every sample, as
[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md)
does, when the reference diversity \\\delta D\_{\sigma}(\tilde{f})\\ of
every unit is already known (e.g. from an external database).

## Usage

``` r
relative.multiplicity.ref_div(
  abund,
  diss,
  clust,
  ref_div,
  sigma = 1,
  cap_at = NULL,
  average_over = "sample",
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

- ref_div:

  Numeric vector of reference diversities, one per unit. They must be
  strictly positive, except with `cap_at`, where they may also be 0: the
  ratio of such a unit is `cap_at` in the samples where it is present,
  and 0 where it is absent. Such units only count with `average_over`
  `"sample"` or `"union"`, and only where they are present. If named,
  values are matched to units by name (extra names are ignored, and
  numeric names are matched as numbers, as for `sigma`); otherwise they
  follow the order of `unique(clust)`. They should be computed with the
  same `sigma` as the one given here, e.g. with
  [`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md).

- sigma:

  Numeric cutoff \\\sigma\\ (default `1`) at which two subunits are
  considered maximally different. Either a single value or a vector of
  strictly positive values with one value per unit. If named, values are
  matched to units by name (numeric names are matched as numbers, so
  `"1e+05"` matches the unit `100000`); otherwise they follow the order
  of `unique(clust)`.

- cap_at:

  `NULL` (default) or a single positive number. If `NULL`, the ratios
  are not capped. Otherwise, every unit's ratio higher than `cap_at` is
  set to `cap_at`, e.g. `cap_at = 1` for the cases where the unit is
  more diverse in the sample than in its reference.

- average_over:

  The units every sample averages over, one of `"sample"` (default, the
  units present in the sample), `"reference"` (the units with a strictly
  positive reference diversity), `"intersection"` (the units present in
  the sample with a strictly positive reference diversity) or `"union"`
  (the units present in the sample or with a strictly positive reference
  diversity). Absent units contribute 0, and a sample without units to
  average over gets 0.

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

Named numeric vector with the relative multiplicity \\RM\\ of every
sample. Samples with zero total abundance get 0.

## Details

`sigma` is only used for the diversities of the samples. As in
[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md),
the units every sample averages over are set by `average_over`; units
with zero abundance in every sample count as 0 with `"reference"` and
`"union"` if their reference diversity is strictly positive. With zero
reference diversities, the options differ: e.g. a sample where only
units with zero reference diversity are present gets 0 with
`"reference"` and `"intersection"`.

## See also

[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md),
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)

## Examples

``` r
clust <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B")
diss <- data.frame(
  ID1 = c("a1", "a1", "a2", "b1"),
  ID2 = c("a2", "a3", "a3", "b2"),
  Distance = c(0.3, 0.8, 0.6, 0.4)
)
abund <- matrix(
  c(5, 3, 0, 2, 0,
    0, 1, 4, 2, 6),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("S1", "S2"), names(clust))
)

# Reference diversities from an external database
relative.multiplicity.ref_div(abund, diss, clust, ref_div = c(A = 2.5, B = 1.8))
#>        S1        S2 
#> 0.5105051 0.5743221 

# Unit B without reference diversity: its ratio is cap_at where present
relative.multiplicity.ref_div(abund, diss, clust, ref_div = c(A = 2.5, B = 0), cap_at = 1)
#>        S1        S2 
#> 0.7327273 0.7475248 

# Only the units with a reference diversity count
relative.multiplicity.ref_div(
  abund, diss, clust, ref_div = c(A = 2.5, B = 0), cap_at = 1, average_over = "reference"
)
#>        S1        S2 
#> 0.4654545 0.4950495 
```
