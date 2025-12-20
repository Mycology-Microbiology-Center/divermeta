# Rao quadratic entropy (Q) (Rao 1982)

Computes Rao quadratic entropy \\Q\\, a diversity measure that accounts
for element abundances and their pairwise dissimilarities.

## Usage

``` r
raoQuadratic(ab, diss)
```

## Arguments

- ab:

  Numeric vector of element abundances.

- diss:

  Numeric square matrix or `dist` object of pairwise dissimilarities
  among elements.

## Value

Numeric scalar, Rao quadratic entropy \\Q\\.

## References

- Rao CR (1982) Diversity and dissimilarity coefficients: A unified
  approach. Theoretical Population Biology, 21.
  [doi:10.1016/0040-5809(82)90004-1](https://doi.org/10.1016/0040-5809%2882%2990004-1)
  .
