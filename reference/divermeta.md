# Estimate diversity and multiplicity indices across samples

Estimate several diversity and multiplicity indices for many samples at
once. Abundances are provided with samples in rows and subunits
(species/OTUs/genes) in columns. Dissimilarities are provided as a
three-column table of pairs of subunits.

## Usage

``` r
divermeta(
  abund,
  diss = NULL,
  indices = c("multiplicity_inventory"),
  clust = NULL,
  q = 1,
  sig = 1,
  method = "sigma",
  diss_clust = NULL,
  average_over = "sample",
  cap_at = NULL,
  percentile = 100,
  assume_max_reference_distance = FALSE,
  assume_homogeneous_abundance = FALSE,
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

  Optional data frame with three columns (taken by position): the
  identifiers of two subunits and the dissimilarity between them, or the
  path to a file with them (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  Every pair of subunits present in at least one sample must be listed
  (see
  [`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
  for the pairs needed with `method = "sigma"`); pairs with a subunit
  that has zero abundance in every sample may be missing, with a
  warning, except when `relative_multiplicity` is requested with
  `assume_homogeneous_abundance = TRUE` and
  `assume_max_reference_distance = FALSE`, where every subunit counts in
  the reference. A file is read only once for all the requested indices.

- indices:

  Character vector of index names to compute (see list above).

- clust:

  Optional vector/factor of cluster memberships for each subunit. If
  named, names are subunit identifiers; otherwise it follows the columns
  of `abund`. A data frame is also accepted, with two columns (subunit
  identifiers, units) or one column of units (named by its row names, if
  set).

- q:

  Numeric order for Hill-number based indices (used by `FD_q` and
  `multiplicity_inventory`). Default `1`.

- sig:

  Numeric cutoff `σ` for distance-based measures (used by `FD_sigma`,
  `multiplicity_distance` and `relative_multiplicity`). Default `1`.

- method:

  Character string specifying the linkage method for aggregating
  distances between clusters. One of: "average", "min", "max", "sigma"
  or "custom". Sigma sets all distances between clusters to be of value
  `sig`. Default is "sigma".

- diss_clust:

  Optional three-column table of distances between clusters, used by
  `multiplicity_distance` when `method` is "custom".

- average_over:

  The units every sample averages over, used by `relative_multiplicity`:
  `"sample"` (default, the units present in the sample), `"reference"`
  (the units with a strictly positive reference diversity),
  `"intersection"` or `"union"` of both. Absent units contribute 0. See
  [`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md).
  `average_redundancy` always averages only the units present in each
  sample.

- cap_at:

  `NULL` (default) or a single positive number, used by
  `relative_multiplicity` and `proportional_abundance`. If set, each
  unit's ratio higher than `cap_at` is set to `cap_at`.

- percentile:

  Numeric in (0, 100\] (default `100`), used by
  `proportional_abundance`. Percentile of the nonzero abundances of
  every unit across samples used as its reference level.

- assume_max_reference_distance:

  Logical (default `FALSE`), used by `relative_multiplicity`. If `TRUE`,
  every pair of subunits in the reference of a unit is at distance
  `sig`, instead of their listed distance.

- assume_homogeneous_abundance:

  Logical (default `FALSE`), used by `relative_multiplicity`. If `TRUE`,
  all subunits of a unit have the same abundance in its reference,
  instead of their abundances pooled across all samples.

- normalize:

  Logical indicating whether index values should be normalized to \[0,
  1\] by dividing each column by its maximum value. Useful for comparing
  indices with different scales. Default `FALSE`.

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

data.frame with one row per sample: a first column `Sample` with the
sample identifiers (the row names of `abund`, or the row numbers when it
has none), then one column per requested index. Empty samples (all zero
abundances) return 0 for all indices.

## Details

Indices available (use these labels in `indices`):

- "multiplicity_inventory": inventory multiplicity \\^{q}M\\ (order
  \\q\\)

- "multiplicity_distance": distance-based multiplicity \\\delta
  M\_{\sigma}\\ (cutoff `sig`)

- "raoQ": Rao quadratic entropy \\Q\\

- "FD_sigma": distance-based functional diversity \\\delta D\_{\sigma}\\
  (cutoff `sig`)

- "FD_q": distance-based functional diversity \\^{q}FD\\ (order \\q\\)

- "redundancy": functional redundancy \\Re\\

- "relative_multiplicity": relative multiplicity \\RM\\ (cutoff `sig`)

- "average_redundancy": average normalized redundancy \\\overline{Re}\\

- "proportional_abundance": proportional abundance of the units
  (`percentile`, `cap_at`)

Notes:

- Indices that use dissimilarities (`multiplicity_distance`, `raoQ`,
  `FD_sigma`, `FD_q`, `redundancy`, `relative_multiplicity`,
  `average_redundancy`) require `diss`.

- Indices that use clustering (`multiplicity_inventory`,
  `multiplicity_distance`, `relative_multiplicity`,
  `average_redundancy`, `proportional_abundance`) require `clust`.

- `proportional_abundance` is computed on the abundances of the units
  (the sum of their subunits) as in
  [`proportional.abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.md);
  it does not require `diss`. Known reference abundances are available
  in
  [`proportional.abundance.ref_abu()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.ref_abu.md).

- `relative_multiplicity` and `average_redundancy` are computed as in
  [`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md)
  and
  [`average.redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/average.redundancy.md)
  (with `normalize = TRUE`), which also offer a cutoff per unit, known
  reference diversities
  ([`relative.multiplicity.ref_div()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.ref_div.md))
  and the unnormalized redundancy. Subunits listed in a named `clust`
  that are not in `abund` have zero abundance in every sample: they
  count in the units of these two indices and do not change the others.

- Every entry of `indices` gives one column, named after the index
  (aliases such as `"M_inventory"` are renamed to
  `"multiplicity_inventory"`). An index requested more than once, by the
  same label or through an alias, gives one identical column per
  request, and the repeated columns get unique names from
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html) (e.g. `raoQ`
  and `raoQ.1`).

## See also

[`multiplicity.inventory()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.inventory.md),
[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md),
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md),
[`diversity.functional.traditional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.traditional.md),
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md),
[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md),
[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md),
[`average.redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/average.redundancy.md),
[`proportional.abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.md)

## Examples

``` r

# Abundance table: rows = samples, columns = subunits
abund <- matrix(
  c(
    10, 0, 4, 6, # Sample1
    5,  8, 3, 0, # Sample2
    0, 12, 7, 9  # Sample3
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("Sample1", "Sample2", "Sample3"), c("s1", "s2", "s3", "s4"))
)

# Dissimilarities among subunits (0-1), one row per pair
diss <- data.frame(
  ID1 = c("s1", "s1", "s1", "s2", "s2", "s3"),
  ID2 = c("s2", "s3", "s4", "s3", "s4", "s4"),
  Distance = c(0.2, 0.7, 1.0, 0.6, 0.9, 0.4)
)

# Clusters per subunit (used by multiplicity_* indices)
clust <- c(s1 = "A", s2 = "A", s3 = "B", s4 = "B")

## Run divermeta
divermeta(abund,
  clust = clust,
  diss = diss,
  q = 1,
  sig = 0.8,
  indices = c("multiplicity_inventory", "multiplicity_distance", "raoQ", "redundancy")
)
#>    Sample multiplicity_inventory multiplicity_distance      raoQ redundancy
#> 1 Sample1               1.400047              1.136364 0.4880000  0.1320000
#> 2 Sample2               1.718327              1.126582 0.2570313  0.3601562
#> 3 Sample3               1.479358              1.186944 0.4408163  0.2096939

## Indices that average over the clusters
divermeta(abund,
  clust = clust,
  diss = diss,
  sig = 0.8,
  indices = c("relative_multiplicity", "average_redundancy")
)
#>    Sample relative_multiplicity average_redundancy
#> 1 Sample1             0.9323921                0.3
#> 2 Sample2             0.8728202                0.4
#> 3 Sample3             0.9363820                0.3
```
