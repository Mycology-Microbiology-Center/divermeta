# Check a large distance file

Checks that a distance file lists every needed pair of subunits exactly
once, with a valid distance, without loading the file into memory. The
indices run this check when called with
`check_large_distance_file = TRUE`; it can also be run on its own, e.g.
once, before computing several indices from the same file.

## Usage

``` r
check_distance_file(
  file,
  ids,
  clust = NULL,
  scope = "all",
  chunk_size = 1e+06,
  max_pairs_in_memory = 1e+07,
  tmp_dir = tempdir(),
  header = TRUE
)
```

## Arguments

- file:

  Path to the distance file, or several paths read one after the other.
  See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)
  for the format.

- ids:

  Character vector with the identifiers of the subunits, e.g.
  `colnames(abund)`. Rows with other identifiers are ignored.

- clust:

  Vector of cluster memberships of the subunits: named by subunit
  identifier, or following the order of `ids`. A data frame is also
  accepted, with two columns (subunit identifiers, units) or one column
  of units. Only needed for `scope` "within" or "between".

- scope:

  Which pairs of subunits must be listed: "all" (default) every pair, as
  needed by
  [`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md),
  [`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md),
  [`diversity.functional.traditional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.traditional.md),
  [`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)
  and
  [`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
  with a linkage method; "within" the pairs of subunits of the same
  unit, as needed by
  [`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
  with `method = "sigma"`,
  [`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md),
  [`relative.multiplicity.ref_div()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.ref_div.md)
  and
  [`average.redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/average.redundancy.md);
  "between" the pairs of subunits of different units, as needed by
  [`unit_distances()`](https://mycology-microbiology-center.github.io/divermeta/reference/unit_distances.md).
  [`divermeta()`](https://mycology-microbiology-center.github.io/divermeta/reference/divermeta.md)
  needs "all" when any of the requested indices does, and "within"
  otherwise. Rows of pairs outside the scope are ignored.

- chunk_size:

  Number of rows read at a time (default `1e6`).

- max_pairs_in_memory:

  Maximum number of pairs expected in memory at a time (default `1e7`,
  about 160 MB). Lower values use less memory and more temporary files.

- tmp_dir:

  Directory for the temporary files (default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- header:

  Logical (default `TRUE`). Whether every file starts with a header
  line. See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md).

## Value

A list of class `divermeta_distance_check`, printed as a summary, with:

- ok:

  `TRUE` if every needed pair is listed exactly once with a valid
  distance.

- scope:

  The scope checked.

- n_rows:

  Number of rows read.

- n_expected:

  Number of pairs needed.

- n_used:

  Number of valid rows of pairs in the scope.

- n_ignored:

  Number of rows ignored: unknown identifiers, self pairs and pairs
  outside the scope.

- n_invalid:

  Number of rows with a non-numeric distance, and of rows in the scope
  with an NA, infinite or negative distance. NA, infinite and negative
  distances outside the scope are ignored, as the indices do.

- n_missing:

  Number of needed pairs with no row at all (a row with an invalid
  distance lists its pair, and is counted in `n_invalid`).

- n_duplicated:

  Number of rows with a valid distance listing a pair already listed
  with a valid distance (in either orientation).

- n_conflicting:

  Number of those rows with a different distance.

- examples:

  An example of every problem found.

## Details

The file is read twice, never as a whole:

1.  It is read in chunks of `chunk_size` rows. Every row of a needed
    pair is written to a temporary binary file on disk (a bucket),
    chosen by its first subunit so that every bucket holds about
    `max_pairs_in_memory` pairs.

2.  Every bucket is loaded in turn to find the duplicated and
    conflicting pairs and to count the pairs of every subunit, which
    gives the missing pairs.

The temporary files take about 16 bytes of disk per pair (e.g. about 80
GB for \\5 \cdot 10^9\\ pairs) and are deleted at the end.

## See also

[distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)

## Examples

``` r
diss <- data.frame(
  ID1 = c("a", "a", "b", "b"),
  ID2 = c("b", "c", "c", "a"),
  Distance = c(0.4, 0.9, 0.6, 0.5)
)
diss_file <- tempfile(fileext = ".tsv")
write.table(diss, diss_file, sep = "\t", quote = FALSE, row.names = FALSE)

# The pair a-b is listed twice, with different distances
check_distance_file(diss_file, ids = c("a", "b", "c"))
#> Distance file check (every pair of subunits): FAILED
#>   Rows read: 4
#>   Pairs needed: 3
#>   Rows used: 4
#>   Rows ignored: 0
#>   Invalid distances: 0
#>   Missing pairs: 0
#>   Duplicated pairs: 1 (e.g. a-b)
#>   Conflicting duplicates: 1 (e.g. a-b (0.4 and 0.5))
```
