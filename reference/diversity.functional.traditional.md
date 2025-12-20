# Distance-based functional diversity (order q) (Chiu & Chao 2014)

Computes distance-based functional diversity \\^{q}FD\\ following Chiu &
Chao (2014). This corresponds to D(Q) in their paper and \\^{q}FD\\ in
the divermeta manuscript. For \\q = 1\\, an approximation is used
internally.

## Usage

``` r
diversity.functional.traditional(ab, diss, q = 1)
```

## Arguments

- ab:

  Numeric vector of element abundances.

- diss:

  Numeric matrix of pairwise dissimilarities among elements.

- q:

  Numeric order of the Hill number (non-negative).

## Value

Numeric scalar, the distance-based functional diversity \\^{q}FD\\.

## References

- Chiu CH, Chao A (2014) Distance-based functional diversity measures
  and their decomposition: A framework based on Hill numbers. PLOS ONE
  9(7).
  [doi:10.1371/journal.pone.0100014](https://doi.org/10.1371/journal.pone.0100014)
  .
  <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014>

## See also

[`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md),
[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
