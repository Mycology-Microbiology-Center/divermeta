# Functional redundancy (Re) by blocks

Computes functional redundancy \\Re\\ from a compact three-column
distance table. This is an efficient implementation of
[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)
for large datasets where storing the full distance matrix would be
memory-intensive. Only the distances between elements inside the same
unit (block) need to be provided; every pair of elements not listed in
the table is assumed to be maximally different (distance equal to 1).

## Usage

``` r
redundancy.by_blocks(ids, ab, diss_frame)
```

## Arguments

- ids:

  Character or integer vector of element identifiers (length `n`). Must
  match the identifiers used in `diss_frame`.

- ab:

  Numeric vector of element abundances (length `n`, same order as
  `ids`).

- diss_frame:

  Data frame with columns `ID1`, `ID2`, `Distance` containing pairwise
  dissimilarities, scaled to the range \[0, 1\], for unique pairs of
  elements. Each unordered pair must be listed only once. Should only
  include within-unit pairs or pairs where distances are less than 1;
  unlisted pairs are automatically set to 1 and distances greater than 1
  are capped at 1.

## Value

Numeric scalar, functional redundancy `Re`.

## Details

Let \\p_i\\ be the relative abundance of element \\i\\ and \\P\\ the set
of unordered pairs \\(i, j)\\, \\i \neq j\\, listed in `diss_frame`.
Assuming every pair not in \\P\\ is at distance 1, Rao's quadratic
entropy is \$\$Q = \left(1 - \sum_i p_i^2\right) - 2 \sum\_{(i,j) \in P}
p_i p_j + 2 \sum\_{(i,j) \in P} p_i p_j d\_{ij},\$\$ so that the
redundancy \\Re = \left(1 - \sum_i p_i^2\right) - Q\\ reduces to \$\$Re
= 2 \sum\_{(i,j) \in P} p_i p_j \left(1 - d\_{ij}\right).\$\$ Only the
listed (within-unit) pairs contribute, so the computation scales with
the number of rows of `diss_frame` instead of \\n^2\\.

## References

- Ricotta C, Pavoine S (2025) What do functional diversity, redundancy,
  rarity, and originality actually measure? A theoretical guide for
  ecologists and conservationists. Ecological Complexity 61.
  [doi:10.1016/j.ecocom.2025.101116](https://doi.org/10.1016/j.ecocom.2025.101116)
  .
  <https://www.sciencedirect.com/science/article/pii/S1476945X25000017>

- Rao CR (1982) Diversity and dissimilarity coefficients: A unified
  approach. Theoretical Population Biology 21.
  [doi:10.1016/0040-5809(82)90004-1](https://doi.org/10.1016/0040-5809%2882%2990004-1)
  .

## See also

[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)
for the standard implementation using full matrices,
[`multiplicity.distance.by_blocks()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.by_blocks.md)
for distance-based multiplicity by blocks,
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
for Rao's quadratic entropy

## Examples

``` r
# Example: Compute redundancy from a distance table
ids <- c("elem1", "elem2", "elem3", "elem4")
ab <- c(10, 15, 20, 25)

# Distance table: only within-unit pairs (other pairs assumed = 1)
diss_frame <- data.frame(
  ID1 = c("elem1", "elem3"),
  ID2 = c("elem2", "elem4"),
  Distance = c(0.3, 0.4),
  stringsAsFactors = FALSE
)

redundancy.by_blocks(ids, ab, diss_frame)
#> [1] 0.1653061
```
