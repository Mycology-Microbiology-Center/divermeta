# Proportional abundance

Computes the proportional abundance of every sample: the abundance of
every element in the sample, divided by a reference level of that
element across all samples, added up over the elements. The reference
level of an element is a percentile of its nonzero abundances in all
samples, so every element is measured on its own scale and abundant
elements do not dominate the index.

## Usage

``` r
proportional.abundance(abund, percentile = 100, cap_at = NULL)
```

## Arguments

- abund:

  Numeric matrix, data.frame or sparse `Matrix` of abundances with
  samples in rows and elements (units or subunits) in columns. Column
  names are the element identifiers.

- percentile:

  Single number in (0, 100\] (default `100`, the maximum). Percentile of
  the nonzero abundances of every element across the samples used as its
  reference level.

- cap_at:

  `NULL` (default) or a single positive number. If `NULL`, the ratios
  are not capped. Otherwise, every element's ratio higher than `cap_at`
  is set to `cap_at`. With `percentile = 100` every ratio is at most 1,
  so only a `cap_at` below 1 changes the result.

## Value

Named numeric vector with the proportional abundance of every sample.
Samples with zero total abundance get 0.

## Details

The function does not know whether the columns of `abund` are units or
subunits: it treats every column as one element.

Let \\x\_{ij}\\ be the abundance of element \\j\\ in sample \\i\\, and
\\P_j\\ the `percentile` (type 7 quantile, see
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html)) of its
nonzero abundances across samples. Then \$\$PA_i = \sum\_{j \\:\\
x\_{ij} \> 0} \frac{x\_{ij}}{P_j},\$\$ or, with a numeric `cap_at`
\\c\\, \$\$PA_i = \sum\_{j \\:\\ x\_{ij} \> 0}
\min\left(\frac{x\_{ij}}{P_j}, c\right).\$\$

With `percentile = 100`, \\P_j\\ is the maximum abundance of the
element, every ratio is in (0, 1\] and the index is at most the number
of elements present. Elements absent from all samples have no reference
level and contribute 0. The index does not change if all abundances are
multiplied by the same constant.

## See also

[`proportional.abundance.ref_abu()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.ref_abu.md)
when the reference abundances are known,
[`multiplicity.inventory()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.inventory.md),
[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md)

## Examples

``` r
abund <- matrix(
  c(10, 0, 4, 6,
    5,  8, 3, 0,
    0, 12, 7, 9),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("S1", "S2", "S3"), c("e1", "e2", "e3", "e4"))
)

proportional.abundance(abund)
#>       S1       S2       S3 
#> 2.238095 1.595238 3.000000 
proportional.abundance(abund, percentile = 50)
#>       S1       S2       S3 
#> 3.133333 2.216667 4.150000 
proportional.abundance(abund, percentile = 50, cap_at = 1)
#>       S1       S2       S3 
#> 2.800000 2.216667 3.000000 
```
