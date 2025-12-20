# Quadratic Form for a Distance Object

Computes the quadratic form \\\sum\_{i,j} 2 p_i p_j m\_{ij}\\ for a
vector `p` and a `dist` object `m`. This is equivalent to \\P^T M P\\,
where `M` is the full distance matrix, without explicitly constructing
the full matrix.

## Usage

``` r
dist_quadratic_form(p, diss)
```

## Arguments

- p:

  Numeric vector of length \\ n \\, representing the weights or
  abundances.

- diss:

  A `dist` object of size \\ n \\ (as returned by `vegdist` or `dist`),
  containing the pairwise distances between observations.

## Value

A numeric scalar, the value of the quadratic form \\P^T M P\\.

## Details

The function uses the fact that a `dist` object stores only the lower
triangle of the distance matrix in column-major order. It removes zero
entries and then constructs all pairs of indices corresponding to this
storage order and computes the sum efficiently.

## Examples

``` r
set.seed(123)
p <- runif(5)  # abundances
diss <- stats::dist(matrix(runif(5 * 3), ncol = 3)) # distance matrix
dist_quadratic_form(p, diss)
#> [1] 5.89332
```
