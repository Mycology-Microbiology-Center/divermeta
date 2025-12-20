# Distance-based multiplicity

Computes distance-based multiplicity \\\delta M\_{\sigma}\\ as the ratio
of distance-based diversity before vs after clustering, using a cutoff
\\\sigma\\ to cap pairwise distances. This index quantifies the
diversity lost when clustering elements into operational units,
incorporating pairwise dissimilarities (e.g., genetic, functional, or
phylogenetic distances) among elements.

## Usage

``` r
multiplicity.distance(
  ab,
  diss,
  clust,
  method = "sigma",
  sig = 1,
  clust_ids_order = NULL,
  diss_clust = NULL
)
```

## Arguments

- ab:

  Numeric vector of element abundances before clustering.

- diss:

  Numeric matrix or `dist` object of pairwise dissimilarities among
  elements before clustering. Must be square and match the length of
  `ab`.

- clust:

  Vector or factor of cluster memberships for each element of `ab`. Must
  have the same length as `ab`.

- method:

  Character string specifying the linkage method for aggregating
  distances between clusters. One of: "average", "min", "max", or
  "sigma". Sigma sets all distances between clusters to be of value
  `sig`. Default is "sigma".

- sig:

  Numeric cutoff \\\sigma\\ (default `1`) defining the threshold
  distance at which two units are considered maximally different. All
  distances greater than `sig` are capped at `sig`. This parameter
  should be chosen based on the biological meaning of distances in your
  dataset (e.g., maximum expected genetic distance, functional
  dissimilarity threshold).

- clust_ids_order:

  Integer or factor vector specifying the order of cluster IDs as they
  appear in the output. If `NULL` (default), uses the natural order of
  unique cluster IDs found in `clust`.

- diss_clust:

  Numeric matrix or `dist` object of pairwise dissimilarities among
  clusters after clustering. Only used when `method` is "custom".
  Default is `NULL`.

## Value

Numeric scalar, the distance-based multiplicity \\\delta M\_{\sigma}\\.
This value represents the effective number of functionally distinct
elements per cluster. A value of 1 indicates minimal functional
diversity within clusters, while higher values indicate greater
functional diversity lost through clustering.

## Details

Unlike inventory multiplicity, distance-based multiplicity is only
defined for \\q = 1\\ (Shannon-type weighting), meaning all elements are
weighted proportionally to their abundance.

The function supports several methods for computing inter-cluster
distances:

- `"sigma"`: All inter-cluster distances are set to `sig`, creating
  maximally different clusters (default)

- `"average"`: Mean distance between all cross-cluster element pairs

- `"min"`: Minimum distance between any cross-cluster element pair

- `"max"`: Maximum distance between any cross-cluster element pair

- `"custom"`: Use a pre-computed cluster distance matrix provided via
  `diss_clust`

When using the "custom" method, you must provide `diss_clust` and ensure
its dimensions match the number of clusters. The `clust_ids_order`
parameter determines how cluster abundances are aligned with the
distance matrix.

## See also

[`multiplicity.distance.by_blocks()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.by_blocks.md)
for computing from a compact distance table,
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
for Rao's quadratic entropy,
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)
for distance-based functional diversity,
[`cluster_distance_matrix()`](https://mycology-microbiology-center.github.io/divermeta/reference/cluster_distance_matrix.md)
for computing cluster distance matrices

## Examples

``` r
# Example 1: Using sigma method (default)
# Two clusters with maximally different clusters
ab <- c(10, 10, 10, 10)  # 4 elements
clust <- c(1, 1, 2, 2)    # 2 clusters
diss <- matrix(c(
  0.0, 0.3, 1.0, 1.0,
  0.3, 0.0, 1.0, 1.0,
  1.0, 1.0, 0.0, 0.4,
  1.0, 1.0, 0.4, 0.0
), nrow = 4, byrow = TRUE)

multiplicity.distance(ab, diss, clust, method = "sigma", sig = 1)
#> [1] 1.212121

# Example 2: Using average linkage
# Compute cluster distances as mean of cross-cluster element distances
multiplicity.distance(ab, diss, clust, method = "average", sig = 1)
#> [1] 1.636364

# Example 3: Using custom cluster distances
# Provide pre-computed cluster distance matrix
ab <- c(10, 10, 10, 10)
clust <- c("A", "A", "B", "B")
clust_ids <- c("A", "B")
diss <- dist(matrix(rnorm(4 * 3), ncol = 3))
diss_clust <- matrix(c(0, 0.8, 0.8, 0), nrow = 2)

multiplicity.distance(
  ab, diss, clust,
  method = "custom",
  diss_clust = diss_clust,
  clust_ids_order = clust_ids
)
#> [1] 2.4

# Example 4: Low diversity within clusters
# Elements within clusters are functionally similar (distance ~0)
diss_low <- matrix(c(
  0.0, 0.05, 1.0, 1.0,
  0.05, 0.0, 1.0, 1.0,
  1.0, 1.0, 0.0, 0.05,
  1.0, 1.0, 0.05, 0.0
), nrow = 4, byrow = TRUE)

multiplicity.distance(ab, diss_low, clust, method = "sigma", sig = 1)
#> [1] 1.025641
```
