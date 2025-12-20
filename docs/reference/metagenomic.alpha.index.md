# Metagenomic Alpha-Diversity Index (MAD) (Finn 2024)

Computes the Metagenomic Alpha-Diversity Index (MAD), a metric that
measures the average dissimilarity of elements (e.g., protein-encoding
genes) within clusters relative to cluster representatives. Unlike
multiplicity, MAD does not account for element abundances and decreases
as the number of elements per cluster increases.

## Usage

``` r
metagenomic.alpha.index(clust, diss, representatives = NULL)
```

## Arguments

- clust:

  Vector or factor of cluster memberships for each element (e.g., gene).
  Must have the same length as the number of rows/columns in `diss`.

- diss:

  Numeric square matrix of pairwise dissimilarities among elements.
  Should be scaled to the range \[0, 1\], where 0 indicates identical
  elements and 1 indicates maximally different elements.

- representatives:

  Optional named vector mapping cluster names to the index (position) of
  the representative element for each cluster. If `NULL` (default), the
  first element in each cluster is used as the representative.

## Value

Numeric scalar, the Metagenomic Alpha-Diversity Index (MAD). Higher
values indicate greater average dissimilarity within clusters.

## Details

Note: Unlike multiplicity indices, MAD does not incorporate element
abundances and decreases as cluster size increases, which may not
reflect biological complexity.

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
clust <- c(1, 1, 2, 2, 2, 3)
diss <- matrix(c(
  0.0, 0.2, 0.8, 0.9, 0.9, 0.9,
  0.2, 0.0, 0.8, 0.9, 0.9, 0.9,
  0.8, 0.8, 0.0, 0.1, 0.15, 0.9,
  0.9, 0.9, 0.1, 0.0, 0.12, 0.9,
  0.9, 0.9, 0.15, 0.12, 0.0, 0.9,
  0.9, 0.9, 0.9, 0.9, 0.9, 0.0
), nrow = 6, byrow = TRUE)

# Use default (first element as representative)
metagenomic.alpha.index(clust, diss)
#> [1] 0.5305556

# Specify representatives explicitly
reps <- c("1" = 1, "2" = 3, "3" = 6)
metagenomic.alpha.index(clust, diss, representatives = reps)
#> [1] 0.5305556
```
