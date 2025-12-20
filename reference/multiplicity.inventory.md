# Inventory multiplicity

Computes inventory multiplicity \\^{q}M\\, the within-cluster diversity
component under Hill-number partitioning. It summarizes the average
diversity inside clusters given element abundances and their cluster
memberships. Use `q` to control abundance weighting (e.g., `q = 0`
richness-like, `q = 1` Shannon-type).

## Usage

``` r
multiplicity.inventory(ab, clust, q = 1)
```

## Arguments

- ab:

  Numeric vector of element (subunit) abundances. Elements with zero
  abundance are automatically removed before computation.

- clust:

  Vector or factor of cluster memberships for each element of `ab`. Must
  have the same length as `ab`

- q:

  Numeric order of the Hill number (default `1`). Controls abundance
  weighting:

  `q = 0`

  :   Richness-like weighting (all elements weighted equally)

  `q = 1`

  :   Shannon-type weighting (default, proportional to abundance)

  `q = 2`

  :   Simpson-type weighting (emphasizes abundant elements)

  Must be non-negative.

## Value

Numeric scalar, the inventory multiplicity \\^{q}M\\. This value can be
interpreted as the effective number of equally abundant elements per
cluster. A value of 1 indicates each cluster contains essentially one
element (no diversity lost), while higher values indicate greater
intra-cluster diversity.

## References

- Hill MO (1973) Diversity and evenness: a unifying notation and its
  consequences. Ecology 54(2):427-432.
  [doi:10.2307/1934352](https://doi.org/10.2307/1934352)

- Jost L (2007) Partitioning diversity into independent alpha and beta
  components. Ecology 88(10):2427-2439.
  [doi:10.1890/06-1736.1](https://doi.org/10.1890/06-1736.1)

## See also

[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
for distance-based multiplicity,
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)
for functional diversity indices

## Examples

``` r
# Example 1: High multiplicity (many diverse elements per cluster)
# Three clusters, each with 3 equally abundant elements
ab <- rep(10, 9)
clust <- c(rep(1, 3), rep(2, 3), rep(3, 3))
multiplicity.inventory(ab, clust)
#> [1] 3
# Result: 3 (each cluster has effective diversity of 3)

# Example 2: Low multiplicity (one element per cluster)
# Three clusters, each with one element
ab <- c(10, 20, 30)
clust <- c(1, 2, 3)
multiplicity.inventory(ab, clust)
#> [1] 1
# Result: 1 (no diversity lost, each cluster is a single element)

# Example 3: Unequal abundances within clusters
ab <- c(10, 5, 2,  # cluster 1: unequal abundances
        8, 8, 8,    # cluster 2: equal abundances
        20, 1)      # cluster 3: very unequal abundances
clust <- c(rep(1, 3), rep(2, 3), rep(3, 2))
multiplicity.inventory(ab, clust, q = 1)
#> [1] 2.103125

# Example 4: Using different q values
ab <- rep(1:12)
clust <- c(rep(1, 4), rep(2, 4), rep(3, 4))
multiplicity.inventory(ab, clust, q = 0)  # Richness-like
#> [1] 4
multiplicity.inventory(ab, clust, q = 1)  # Shannon-type (default)
#> [1] 3.914214
multiplicity.inventory(ab, clust, q = 2)  # Simpson-type
#> [1] 3.907692
```
