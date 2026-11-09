# Relative multiplicity

Computes relative multiplicity \\RM\\ of every sample: the average,
across the units present in the sample, of each unit's distance-based
functional diversity \\\delta D\_{\sigma}\\ relative to the diversity of
a reference composition of that same unit. Because every unit is
compared only against itself, units that are not comparable with each
other (e.g. different enzyme families) can be summarized together, and
no unit dominates the index because of its size or abundance. The
reference composition of every unit is built from the whole study, as
set by `assume_max_reference_distance` and
`assume_homogeneous_abundance`, and its diversity is computed only once.

## Usage

``` r
relative.multiplicity(
  abund,
  diss,
  clust,
  sigma = 1,
  cap_at = NULL,
  average_over = "sample",
  assume_max_reference_distance = FALSE,
  assume_homogeneous_abundance = FALSE,
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

- assume_max_reference_distance:

  Logical (default `FALSE`). If `FALSE`, the reference distances of a
  unit are the distances among all its subunits. If `TRUE`, every pair
  of subunits in the reference is at distance \\\sigma\\.

- assume_homogeneous_abundance:

  Logical (default `FALSE`). If `FALSE`, the reference abundances of a
  unit are those of its subunits pooled across all samples. If `TRUE`,
  all subunits of the unit have the same abundance in the reference.

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

Let \\f\\ be the composition of a unit in the sample, \\N_f\\ its total
abundance and \\\tilde{f}\\ its reference composition. Then \$\$RM =
\frac{1}{\left\|\\f : N_f \> 0\\\right\|} \sum\_{f \\:\\ N_f \> 0}
\frac{\delta D\_{\sigma}(f)}{\delta D\_{\sigma}(\tilde{f})},\$\$ with
\\\delta D\_{\sigma}\\ computed as in
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md).

By default (`average_over = "sample"`), only units present in the sample
are averaged: a unit with zero total abundance is ignored. In general,
the average is taken over a set of units \\U\\ set by `average_over`,
\$\$RM = \frac{1}{\|U\|} \sum\_{f \in U} \frac{\delta
D\_{\sigma}(f)}{\delta D\_{\sigma}(\tilde{f})},\$\$ where absent units
contribute 0, and \\RM = 0\\ if \\U\\ is empty.

Every unit that occurs in the study has a strictly positive reference
diversity, so here `"intersection"` is the same as `"sample"` and
`"union"` the same as `"reference"`. With the pooled reference
abundances (`assume_homogeneous_abundance = FALSE`), a unit with zero
abundance across the whole study has no reference and never counts; with
`assume_homogeneous_abundance = TRUE`, it has a reference and counts as
0 with `"reference"` and `"union"`.

With `assume_max_reference_distance = TRUE`, the reference diversity of
a unit is the inverse Simpson index of its reference abundances; if
`assume_homogeneous_abundance` is also `TRUE`, it is simply the number
of subunits of the unit.

Every unit's diversity depends on the distances only through the sum
over its pairs of subunits \\\sum\_{i\<j} a_i a_j (\min(d\_{ij},
\sigma) - \sigma)\\, accumulated while reading the distance table. The
samples and the reference are accumulated in the same read.

## See also

[`relative.multiplicity.ref_div()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.ref_div.md)
when the reference diversities are already known,
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md),
[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)

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

relative.multiplicity(abund, diss, clust)
#>        S1        S2        S3 
#> 0.7803188 0.8759285 0.7278825 
relative.multiplicity(
  abund, diss, clust,
  assume_max_reference_distance = TRUE,
  assume_homogeneous_abundance = TRUE
)
#>        S1        S2        S3 
#> 0.4439394 0.5003883 0.4166667 
```
