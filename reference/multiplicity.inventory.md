# Inventory multiplicity

Computes inventory multiplicity \\^{q}M\\, the within-cluster diversity
component under Hill-number partitioning, for every sample. It
summarizes the average diversity inside clusters given element
abundances and their cluster memberships. Use `q` to control abundance
weighting (e.g., `q = 0` richness-like, `q = 1` Shannon-type).

## Usage

``` r
multiplicity.inventory(abund, clust, q = 1)
```

## Arguments

- abund:

  Numeric matrix, data.frame or sparse `Matrix` of abundances with
  samples in rows and subunits in columns. Column names are the subunit
  identifiers.

- clust:

  Vector or factor of cluster memberships for each subunit. If named,
  names are subunit identifiers and every column of `abund` must be
  listed; otherwise it must follow the columns of `abund`. A data frame
  is also accepted, with two columns (subunit identifiers, units) or one
  column of units (named by its row names, if set).

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

Named numeric vector with the inventory multiplicity \\^{q}M\\ of every
sample. This value can be interpreted as the effective number of equally
abundant elements per cluster. A value of 1 indicates each cluster
contains essentially one element (no diversity lost), while higher
values indicate greater intra-cluster diversity. Samples with zero total
abundance get 0. Since the index is otherwise at least 1, 0 flags an
empty sample rather than a diversity value.

## Details

Subunits with zero abundance in a sample are ignored for that sample.

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
# Three clusters with three subunits each (a, b, c)
clust <- c(a1 = "a", a2 = "a", a3 = "a", b1 = "b", b2 = "b", b3 = "b",
           c1 = "c", c2 = "c", c3 = "c")
abund <- rbind(
  High = rep(10, 9),                     # every subunit equally abundant
  Low = c(10, 0, 0, 20, 0, 0, 30, 0, 0), # one subunit per cluster
  Uneven = c(10, 5, 2, 8, 8, 8, 20, 1, 0)
)
colnames(abund) <- names(clust)

multiplicity.inventory(abund, clust)          # High: 3, Low: 1
#>     High      Low   Uneven 
#> 3.000000 1.000000 2.103125 
multiplicity.inventory(abund, clust, q = 0)  # Richness-like
#>     High      Low   Uneven 
#> 3.000000 1.000000 2.666667 
multiplicity.inventory(abund, clust, q = 2)  # Simpson-type
#>     High      Low   Uneven 
#> 3.000000 1.000000 1.808864 
```
