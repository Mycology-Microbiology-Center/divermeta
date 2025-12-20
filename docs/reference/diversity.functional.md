# Distance-based functional diversity (q = 1) (Chiu & Chao 2014)

Computes distance-based functional diversity \\\delta D\_{\sigma}\\
following Chiu & Chao (2014) with Shannon-type weighting (order \\q =
1\\). Pairwise distances are capped at the cutoff \\\sigma\\.

## Usage

``` r
diversity.functional(ab, diss, sig = 1)
```

## Arguments

- ab:

  Numeric vector of element abundances.

- diss:

  Numeric matrix or `dist` object of pairwise dissimilarities among
  elements.

- sig:

  Numeric cutoff \\\sigma\\ at which two units are considered different
  (default `1`).

## Value

Numeric scalar, the distance-based functional diversity \\\delta
D\_{\sigma}\\.

## References

- Chiu CH, Chao A (2014) Distance-based functional diversity measures
  and their decomposition: A framework based on Hill numbers. PLOS ONE
  9(7).
  [doi:10.1371/journal.pone.0100014](https://doi.org/10.1371/journal.pone.0100014)
  .
  <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014>

## See also

[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md),
[`diversity.functional.traditional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.traditional.md)
