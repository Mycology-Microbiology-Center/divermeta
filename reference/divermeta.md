# Estimate diversity and multiplicity indices across samples

Estimate several diversity and multiplicity indices for many samples at
once. Abundances are provided with samples in columns and features
(species/OTUs/genes) in rows. If supplied, the dissimilarity matrix must
correspond to feature names.

## Usage

``` r
divermeta(
  abund,
  diss = NULL,
  indices = c("multiplicity_inventory"),
  clusters = NULL,
  q = 1,
  sig = 1,
  method = "sigma",
  normalize = FALSE
)
```

## Arguments

- abund:

  Numeric matrix or data.frame of abundances with samples in columns and
  features in rows.

- diss:

  Optional numeric dissimilarity matrix or `dist` object among features;
  square with row/column names matching `rownames(abund)` when present.

- indices:

  Character vector of index names to compute (see list above).

- clusters:

  Optional vector/factor of cluster memberships for each feature (row of
  `abund`).

- q:

  Numeric order for Hill-number based indices (used by `FD_q` and
  `multiplicity_inventory`). Default `1`.

- sig:

  Numeric cutoff `σ` for distance-based measures (used by `FD_sigma` and
  `multiplicity_distance`). Default `1`.

- method:

  Character string specifying the linkage method for aggregating
  distances between clusters. One of: "average", "min", "max", or
  "sigma". Sigma sets all distances between clusters to be of value
  `sig`. Default is "sigma".

- normalize:

  Logical indicating whether index values should be normalized to \[0,
  1\] by dividing each column by its maximum value. Useful for comparing
  indices with different scales. Default `FALSE`.

## Value

data.frame with one row per sample and one column per requested index.
Empty samples (all zero abundances) return `NA` for all indices.

## Details

Indices available (use these labels in `indices`):

- "multiplicity_inventory": inventory multiplicity \\^{q}M\\ (order
  \\q\\)

- "multiplicity_distance": distance-based multiplicity \\\delta
  M\_{\sigma}\\ (cutoff `sig`)

- "raoQ": Rao quadratic entropy \\Q\\

- "FD_sigma": distance-based functional diversity \\\delta D\_{\sigma}\\
  (cutoff `sig`)

- "FD_q": distance-based functional diversity \\^{q}FD\\ (order \\q\\)

- "redundancy": functional redundancy \\Re\\

Notes:

- Indices that use dissimilarities (`multiplicity_distance`, `raoQ`,
  `FD_sigma`, `FD_q`, `redundancy`) require `diss`.

- Indices that use clustering (`multiplicity_inventory`,
  `multiplicity_distance`) require `clusters` (same length/order as
  feature rows).

## See also

[`multiplicity.inventory()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.inventory.md),
[`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md),
[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md),
[`diversity.functional.traditional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.traditional.md),
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md),
[`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)

## Examples

``` r
# Abundance matrix: rows = features, columns = samples
species <- c("s1", "s2", "s3", "s4")
samples <- c("Sample1", "Sample2", "Sample3")
abund <- matrix(
  c(
    10, 5,  0, # s1
    0,  8, 12, # s2
    4,  3,  7, # s3
    6,  0,  9 # s4
  ),
  nrow = length(species),
  ncol = length(samples),
  byrow = TRUE,
  dimnames = list(species, samples)
)

# Dissimilarity matrix among features (0-1), names must match abund rownames
diss <- matrix(
  c(
    0.0, 0.2, 0.7, 1.0,
    0.2, 0.0, 0.6, 0.9,
    0.7, 0.6, 0.0, 0.4,
    1.0, 0.9, 0.4, 0.0
  ),
  nrow = length(species),
  ncol = length(species),
  byrow = TRUE,
  dimnames = list(species, species)
)

# Clusters per feature (used by multiplicity_* indices)
clusters <- c(s1 = "A", s2 = "A", s3 = "B", s4 = "B")

## Run divermeta
divermeta(abund,
  clusters = clusters,
  diss = diss,
  q = 1,
  sig = 0.8,
  indices = c("multiplicity_inventory", "multiplicity_distance", "raoQ", "redundancy")
)
#>    Sample multiplicity_inventory multiplicity_distance      raoQ redundancy
#> 1 Sample1               1.400047              1.075269 0.4880000  0.1320000
#> 2 Sample2               1.718327              1.024460 0.2570312  0.3601563
#> 3 Sample3               1.479358              1.055409 0.4408163  0.2096939
```
