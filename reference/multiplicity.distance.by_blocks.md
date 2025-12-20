# Distance-based multiplicity by blocks

Computes distance-based multiplicity \\\delta M\_{\sigma}\\ from a
compact three-column distance table. This is an efficient implementation
for large datasets where storing the full distance matrix would be
memory-intensive. Distances between elements from different clusters are
assumed to be equal to \\\sigma\\ (maximally different).

## Usage

``` r
multiplicity.distance.by_blocks(ids, ab, diss_frame, clust, sigma = 1)
```

## Arguments

- ids:

  Character or integer vector of element identifiers (length `n`). Must
  match the identifiers used in `diss_frame`.

- ab:

  Numeric vector of pre-clustering abundances (length `n`, same order as
  `ids`).

- diss_frame:

  Data frame with columns `ID1`, `ID2`, `Distance` containing pairwise
  distances for unique pairs of elements. Should only include
  within-cluster pairs or pairs where distances are less than `sigma`;
  cross-cluster distances are automatically set to `sigma`.

- clust:

  Vector or factor of cluster memberships for each element (length `n`,
  same order as `ids`).

- sigma:

  Numeric cutoff \\\sigma\\ (default `1`) at which two units are
  considered maximally different. All distances in `diss_frame` greater
  than `sigma` are capped, and distances between elements from different
  clusters are set to `sigma`.

## Value

Numeric scalar, the distance-based multiplicity \\\delta M\_{\sigma}\\.

## See also

[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
for the standard implementation using full matrices,
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
for Rao's quadratic entropy,
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)
for distance-based functional diversity

## Examples

``` r
# Example: Compute multiplicity from a distance table
ids <- c("elem1", "elem2", "elem3", "elem4")
ab <- c(10, 15, 20, 25)
clust <- c(1, 1, 2, 2)

# Distance table: only within-cluster pairs (cross-cluster assumed = sigma)
diss_frame <- data.frame(
  ID1 = c("elem1", "elem3"),
  ID2 = c("elem2", "elem4"),
  Distance = c(0.3, 0.4),
  stringsAsFactors = FALSE
)

multiplicity.distance.by_blocks(ids, ab, diss_frame, clust, sigma = 1)
#> [1] 1.226852
```
