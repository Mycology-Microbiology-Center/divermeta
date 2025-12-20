# Compute Distance Matrix Between Clusters

Computes a distance matrix between clusters based on pairwise distances
between elements in different clusters. Supports various linkage methods
for aggregating distances between cluster members.

## Usage

``` r
cluster_distance_matrix(
  diss,
  clust,
  clust_ids_order,
  method = "sigma",
  sig = 1
)
```

## Arguments

- diss:

  A distance matrix or `dist` object containing pairwise distances
  between all elements

- clust:

  Integer or factor vector of cluster assignments for each element

- clust_ids_order:

  Integer or factor vector specifying the order of cluster IDs to use in
  the output distance matrix

- method:

  Character string specifying the linkage method for aggregating
  distances between clusters. One of: "average", "min", "max", or
  "sigma". Sigma sets all distances between clusters to be of value
  `sig`. Default is "sigma".

- sig:

  Numeric value specifying the maximum distance threshold when using the
  "sigma" method. Distances greater than sigma are capped at this value.
  Default is 1.

## Value

A distance matrix or `dist` object (matching the input type) containing
distances between clusters. The dimensions correspond to
`length(clust_ids_order)`.

## Details

This function computes inter-cluster distances by aggregating all
pairwise distances between elements belonging to different clusters. The
available linkage methods are:

- `"average"`: Mean distance between all cross-cluster pairs

- `"min"`: Minimum distance between any cross-cluster pair

- `"max"`: Maximum distance between any cross-cluster pair

- `"sigma"`: Returns a matrix where all inter-cluster distances are set
  to `sigma` (useful for creating complete graphs)

The function preserves the input type: if `diss` is a `dist` object, the
output will also be a `dist` object; if it's a full matrix, the output
will be a full matrix.

## Note

The `clust_ids_order` parameter determines the order of clusters in the
output distance matrix. This is useful when you want to maintain a
specific cluster ordering rather than using the natural order of cluster
IDs found in `clust`.

When using the "sigma" method, all inter-cluster distances are set to
the specified sigma value, creating a complete graph where all clusters
are equally distant from each other.

## Examples

``` r
# Create sample data
set.seed(123)
n <- 10
diss <- dist(matrix(rnorm(n * 3), ncol = 3))
clust <- c(1, 1, 2, 2, 3, 3, 4, 4, 1, 2)
clust_ids <- sort(unique(clust))

# Compute cluster distances using different methods
clust_dist_avg <- cluster_distance_matrix(diss, clust, clust_ids, "average")
clust_dist_min <- cluster_distance_matrix(diss, clust, clust_ids, "min")
clust_dist_max <- cluster_distance_matrix(diss, clust, clust_ids, "max")
clust_dist_sigma <- cluster_distance_matrix(diss, clust, clust_ids, "sigma", sig = 2)

# Compare the results
print(as.matrix(clust_dist_avg))
#>          1        2        3        4
#> 1 0.000000 1.863846 2.091610 2.487629
#> 2 1.863846 0.000000 2.130240 2.289524
#> 3 2.091610 2.130240 0.000000 3.049326
#> 4 2.487629 2.289524 3.049326 0.000000
print(as.matrix(clust_dist_min))
#>           1         2         3        4
#> 1 0.0000000 0.6430504 1.0645847 1.269369
#> 2 0.6430504 0.0000000 0.6771221 1.391801
#> 3 1.0645847 0.6771221 0.0000000 1.833059
#> 4 1.2693691 1.3918006 1.8330587 0.000000
print(as.matrix(clust_dist_max))
#>          1        2        3        4
#> 1 0.000000 2.877945 2.824485 3.488312
#> 2 2.877945 0.000000 4.292038 3.868997
#> 3 2.824485 4.292038 0.000000 5.133808
#> 4 3.488312 3.868997 5.133808 0.000000
print(clust_dist_sigma)
#>   1 2 3
#> 2 2    
#> 3 2 2  
#> 4 2 2 2

# Use with full matrix input
diss_matrix <- as.matrix(diss)
clust_dist_matrix <- cluster_distance_matrix(diss_matrix, clust, clust_ids, "average")

```
