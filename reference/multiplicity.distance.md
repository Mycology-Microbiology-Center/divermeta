# Distance-based multiplicity

Computes distance-based multiplicity \\\delta M\_{\sigma}\\ of every
sample as the ratio of distance-based diversity before vs after
clustering, using a cutoff \\\sigma\\ to cap pairwise distances. This
index quantifies the diversity lost when clustering elements into
operational units, incorporating pairwise dissimilarities (e.g.,
genetic, functional, or phylogenetic distances) among elements.

## Usage

``` r
multiplicity.distance(
  abund,
  diss,
  clust,
  method = "sigma",
  sig = 1,
  diss_clust = NULL,
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

  Data frame with three columns (taken by position): the identifiers of
  two subunits and the dissimilarity between them, in either
  orientation, or the path to a file with them (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  With `method = "sigma"` only the pairs of subunits of the same cluster
  are needed (every one of them must be listed); with any other method
  every pair of subunits must be listed. Only the subunits present in at
  least one sample count: pairs with a subunit that has zero abundance
  in every sample may be missing, with a warning (see Details for the
  linkage methods). Rows with identifiers not in `abund` and self pairs
  are ignored, and so are, with `method = "sigma"`, the rows of subunits
  of different clusters, also when their distance is invalid (NA,
  infinite or negative).

- clust:

  Vector or factor of cluster memberships for each subunit. If named,
  names are subunit identifiers and every column of `abund` must be
  listed; otherwise it must follow the columns of `abund`. A data frame
  is also accepted, with two columns (subunit identifiers, units) or one
  column of units (named by its row names, if set).

- method:

  Character string specifying the linkage method for aggregating
  distances between clusters. One of: "average", "min", "max", "sigma"
  or "custom". Sigma sets all distances between clusters to be of value
  `sig`. Default is "sigma".

- sig:

  Numeric cutoff \\\sigma\\ (default `1`) defining the threshold
  distance at which two units are considered maximally different. All
  distances greater than `sig` are capped at `sig`. This parameter
  should be chosen based on the biological meaning of distances in your
  dataset (e.g., maximum expected genetic distance, functional
  dissimilarity threshold).

- diss_clust:

  Data frame with three columns (taken by position): the identifiers of
  two clusters and the dissimilarity between them. Every pair of
  clusters must be listed. Only used when `method` is "custom" (see
  [`unit_distances()`](https://mycology-microbiology-center.github.io/divermeta/reference/unit_distances.md)).
  Default is `NULL`.

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

Named numeric vector with the distance-based multiplicity \\\delta
M\_{\sigma}\\ of every sample. This value represents the effective
number of functionally distinct elements per cluster. A value of 1
indicates minimal functional diversity within clusters, while higher
values indicate greater functional diversity lost through clustering.
Samples with zero total abundance get 0. Since the index is otherwise
always positive, 0 flags an empty sample rather than a multiplicity
value.

## Details

Distance-based multiplicity is the ratio of the distance-based
functional diversities \\\delta D\_{\sigma}\\ (see
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md))
of every sample before and after clustering.

The function supports several methods for computing inter-cluster
distances:

- `"sigma"`: All inter-cluster distances are set to `sig`, creating
  maximally different clusters (default). Subunits of different clusters
  are then also at distance `sig`: their listed distances are ignored.

- `"average"`: Mean distance between all cross-cluster element pairs

- `"min"`: Minimum distance between any cross-cluster element pair

- `"max"`: Maximum distance between any cross-cluster element pair

- `"custom"`: Use pre-computed cluster distances provided via
  `diss_clust`

Cluster distances are computed from the uncapped subunit distances (over
every subunit in `abund`) and then capped at `sig`.

Subunits with zero abundance in every sample do not need distances. With
the linkage methods ("average", "min" and "max"), their distances are
still used for the cluster distances when they are listed, so the result
can change depending on whether they are listed. With a distance file,
duplicated rows of their pairs are not detected, also with
`check_large_distance_file = TRUE`. A cluster whose subunits all have
zero abundance in every sample needs no distance to the other clusters.

## See also

[`unit_distances()`](https://mycology-microbiology-center.github.io/divermeta/reference/unit_distances.md)
for computing cluster distances,
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
for Rao's quadratic entropy,
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)
for distance-based functional diversity

## Examples

``` r
# Two clusters of two subunits
clust <- c(e1 = 1, e2 = 1, e3 = 2, e4 = 2)
abund <- matrix(
  c(10, 10, 10, 10,
    10, 0, 5, 20),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("S1", "S2"), names(clust))
)

# Example 1: Using sigma method (default): only pairs inside clusters are needed
diss_within <- data.frame(ID1 = c("e1", "e3"), ID2 = c("e2", "e4"), Distance = c(0.3, 0.4))
multiplicity.distance(abund, diss_within, clust, method = "sigma", sig = 1)
#>       S1       S2 
#> 1.212121 1.124031 

# Example 2: Using average linkage: every pair is needed
diss <- data.frame(
  ID1 = c("e1", "e1", "e1", "e2", "e2", "e3"),
  ID2 = c("e2", "e3", "e4", "e3", "e4", "e4"),
  Distance = c(0.3, 0.9, 1.0, 0.8, 1.0, 0.4)
)
multiplicity.distance(abund, diss, clust, method = "average", sig = 1)
#>       S1       S2 
#> 1.194444 1.164122 

# Example 3: Using custom cluster distances
diss_clust <- data.frame(ID1 = 1, ID2 = 2, Distance = 0.8)
multiplicity.distance(abund, diss, clust, method = "custom", diss_clust = diss_clust)
#>       S1       S2 
#> 1.333333 1.259542 
```
