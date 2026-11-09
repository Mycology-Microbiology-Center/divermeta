# Metagenomic Alpha-Diversity Index (MAD) (Finn 2024)

Computes the Metagenomic Alpha-Diversity Index (MAD) of every sample, a
metric that measures the average dissimilarity of elements (e.g.,
protein-encoding genes) within clusters relative to cluster
representatives. Unlike multiplicity, MAD does not account for element
abundances and decreases as the number of elements per cluster
increases.

## Usage

``` r
metagenomic.alpha.index(
  abund,
  diss,
  clust,
  representatives = NULL,
  chunk_size = 1e+06,
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
  two subunits and the dissimilarity between them, scaled to the range
  \[0, 1\], or the path to a file with them (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  The distance between every representative and every other subunit of
  its cluster present in a sample must be listed. Only these pairs are
  kept while reading, so missing pairs and pairs listed again with a
  different distance are always detected, also in files (there is no
  need for `check_large_distance_file`).

- clust:

  Vector or factor of cluster memberships for each subunit. If named,
  names are subunit identifiers; otherwise it follows the columns of
  `abund`. A data frame is also accepted, with two columns (subunit
  identifiers, units) or one column of units (named by its row names, if
  set).

- representatives:

  Optional named vector mapping cluster names to the identifier of the
  representative subunit of each cluster. Numeric identifiers are
  matched as numbers, so a name `"1e+05"` matches the cluster `100000`.
  If `NULL` (default), the first subunit (in the column order of
  `abund`) of each cluster present in the sample is used as the
  representative.

- chunk_size:

  Number of rows of a distance file read at a time (default `1e6`).
  Lower values use less memory. Ignored when `diss` is in memory.

- header:

  Logical (default `TRUE`). Whether a distance file starts with a header
  line (with several files, every one of them). With `TRUE`, a first
  line that looks like data (a number as third field) is an error; use
  `FALSE` for files without a header. Ignored when `diss` is in memory.
  See
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md).

## Value

Named numeric vector with the Metagenomic Alpha-Diversity Index (MAD) of
every sample. Higher values indicate greater average dissimilarity
within clusters. Samples with zero total abundance get 0.

## Details

Abundances are only used to decide which subunits are present (nonzero
abundance) in every sample: the index of a sample is computed over its
present subunits and clusters. A given representative is used as the
reference of its cluster even when it is absent from a sample.

Note: Unlike multiplicity indices, MAD does not incorporate element
abundances and decreases as cluster size increases, which may not
reflect biological complexity.

Note on duplicated pairs: an in-memory `diss` is validated as a whole,
so a pair listed again with a different distance is an error even when
MAD does not use that pair. When `diss` is a file, only the pairs MAD
uses (representative and subunit of its cluster) are checked for
conflicting duplicates, and the others are ignored. In both cases every
distance must be a valid (non-negative, not NA) number.

## References

- Finn DR (2024) A metagenomic alpha-diversity index for microbial
  functional biodiversity. FEMS Microbiology Ecology 100(3):fiae019.
  [doi:10.1093/femsec/fiae019](https://doi.org/10.1093/femsec/fiae019)

## See also

[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
for abundance-weighted distance-based multiplicity,
[`multiplicity.inventory()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.inventory.md)
for inventory multiplicity

## Examples

``` r
# Example: Compute MAD for gene clusters
clust <- c(g1 = 1, g2 = 1, g3 = 2, g4 = 2, g5 = 2, g6 = 3)
diss <- data.frame(
  ID1 = c("g1", "g3", "g3", "g4"),
  ID2 = c("g2", "g4", "g5", "g5"),
  Distance = c(0.2, 0.1, 0.15, 0.12)
)
abund <- matrix(
  c(1, 1, 1, 1, 1, 1,
    1, 0, 0, 1, 1, 1),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("S1", "S2"), names(clust))
)

# Use default (first present element as representative)
metagenomic.alpha.index(abund, diss, clust)
#>        S1        S2 
#> 0.5305556 0.7650000 

# Specify representatives explicitly
reps <- c("1" = "g1", "2" = "g3", "3" = "g6")
metagenomic.alpha.index(abund, diss, clust, representatives = reps)
#>        S1        S2 
#> 0.5305556 0.7812500 
```
