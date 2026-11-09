# Functional redundancy (Re) (Ricotta & Pavoine 2025)

Computes functional redundancy \\Re\\, a measure of the degree to which
distinct elements are functionally similar given their abundances and
pairwise dissimilarities, for every sample. This implementation follows
the Simpson–Rao family and corresponds to the `q = 2` case.

## Usage

``` r
redundancy(
  abund,
  diss,
  normalize = FALSE,
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
  two subunits and the dissimilarity between them, scaled to the range
  \[0, 1\] (distances above 1 are capped at 1), or the path to a file
  with them (see
  [distance-files](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)).
  Every pair of subunits present in at least one sample must be listed,
  in either orientation; pairs with a subunit that has zero abundance in
  every sample may be missing, with a warning.

- normalize:

  Logical (default `FALSE`). If `TRUE`, the redundancy is divided by
  Simpson's diversity \\S\\, giving \\(S - Q) / S = 1 - Q / S\\.

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

Named numeric vector with the functional redundancy `Re` of every
sample. Samples with zero total abundance get 0.

## Details

With relative abundances \\p_i\\, Simpson's diversity is \\S = 1 -
\sum_i p_i^2\\ and Rao's quadratic entropy is \\Q = \sum\_{i \neq j} p_i
p_j d\_{ij}\\. Then \$\$Re = S - Q,\$\$ and, with `normalize = TRUE`,
\$\$Re = \frac{S - Q}{S} = 1 - \frac{Q}{S}.\$\$ A sample with a single
present subunit has no redundancy: its value is 0, also when normalized.

Distances greater than 1 are capped at 1 (\\d\_{ij} = \min(d\_{ij},
1)\\) before computing \\Q\\, so that \\Q \le S\\ and the redundancy is
never negative: it lies between 0 and \\S\\, or between 0 and 1 when
normalized.

## References

- Ricotta C, Pavoine S (2025) What do functional diversity, redundancy,
  rarity, and originality actually measure? A theoretical guide for
  ecologists and conservationists. Ecological Complexity 61.
  [doi:10.1016/j.ecocom.2025.101116](https://doi.org/10.1016/j.ecocom.2025.101116)
  .
  <https://www.sciencedirect.com/science/article/pii/S1476945X25000017>

- Ricotta, C., de Bello, F., Moretti, M., Caccianiga, M., Cerabolini,
  B.E.L., & Pavoine, S. (2016). Measuring the functional redundancy of
  biological communities: a quantitative guide. Methods in Ecology and
  Evolution, 7(11), 1386–1395. <https://doi.org/10.1111/2041-210X.12604>

- Rao CR (1982) Diversity and dissimilarity coefficients: A unified
  approach. Theoretical Population Biology 21.
  [doi:10.1016/0040-5809(82)90004-1](https://doi.org/10.1016/0040-5809%2882%2990004-1)
  .

## See also

[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md),
[`average.redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/average.redundancy.md)

## Examples

``` r
abund <- matrix(c(2, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
                dimnames = list(c("S1", "S2"), c("a", "b", "c")))
diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))
redundancy(abund, diss)
#>        S1        S2 
#> 0.2666667 0.2444444 
redundancy(abund, diss, normalize = TRUE)
#>        S1        S2 
#> 0.6000000 0.3666667 
```
