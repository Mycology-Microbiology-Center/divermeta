# Proportional abundance with known reference abundances

Computes the proportional abundance of every sample, as
[`proportional.abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.md)
does, when the reference abundance of every element is already known
(e.g. from an external database or a larger collection of samples)
instead of being a percentile of the given samples.

## Usage

``` r
proportional.abundance.ref_abu(abund, ref_abu, cap_at = NULL)
```

## Arguments

- abund:

  Numeric matrix, data.frame or sparse `Matrix` of abundances with
  samples in rows and elements (units or subunits) in columns. Column
  names are the element identifiers.

- ref_abu:

  Numeric vector of strictly positive reference abundances, one per
  element (column of `abund`). If named, values are matched to the
  columns by name (extra names are ignored); otherwise they follow the
  order of the columns. With `cap_at`, reference abundances may also be
  0: the ratio of such an element is `cap_at` in the samples where it is
  present, and 0 where it is absent.

- cap_at:

  `NULL` (default) or a single positive number. If `NULL`, the ratios
  are not capped. Otherwise, every element's ratio higher than `cap_at`
  is set to `cap_at`. With `percentile = 100` every ratio is at most 1,
  so only a `cap_at` below 1 changes the result.

## Value

Named numeric vector with the proportional abundance of every sample.
Samples with zero total abundance get 0.

## Details

With \\R_j\\ the reference abundance of element \\j\\, \$\$PA_i =
\sum\_{j \\:\\ x\_{ij} \> 0} \frac{x\_{ij}}{R_j},\$\$ or, with a numeric
`cap_at` \\c\\, \$\$PA_i = \sum\_{j \\:\\ x\_{ij} \> 0}
\min\left(\frac{x\_{ij}}{R_j}, c\right).\$\$

Ratios are not bounded by 1: abundances above their reference give
ratios above 1, unless capped. Elements absent from a sample contribute
0 to it.

## See also

[`proportional.abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.md)

## Examples

``` r
abund <- matrix(
  c(10, 0, 4, 6,
    5,  8, 3, 0,
    0, 12, 7, 9),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("S1", "S2", "S3"), c("e1", "e2", "e3", "e4"))
)

# Reference abundances from an external database
proportional.abundance.ref_abu(abund, ref_abu = c(e1 = 8, e2 = 10, e3 = 5, e4 = 12))
#>    S1    S2    S3 
#> 2.550 2.025 3.350 

# Element e4 without reference abundance: its ratio is cap_at where present
proportional.abundance.ref_abu(abund, ref_abu = c(e1 = 8, e2 = 10, e3 = 5, e4 = 0), cap_at = 1)
#>    S1    S2    S3 
#> 2.800 2.025 3.000 
```
