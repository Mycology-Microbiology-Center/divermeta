# Functional redundancy (Re) (Ricotta & Pavoine 2025)

Computes functional redundancy \\Re\\, a measure of the degree to which
distinct elements are functionally similar given their abundances and
pairwise dissimilarities. This implementation follows the Simpson–Rao
family and corresponds to the `q = 2` case.

## Usage

``` r
redundancy(ab, diss)
```

## Arguments

- ab:

  Numeric vector of element abundances.

- diss:

  Numeric square matrix or `dist` object of pairwise dissimilarities
  among elements scaled to the range \[0, 1\].

## Value

Numeric scalar, functional redundancy `Re`.

## References

- Ricotta C, Pavoine S (2025) What do functional diversity, redundancy,
  rarity, and originality actually measure? A theoretical guide for
  ecologists and conservationists. Ecological Complexity 61.
  [doi:10.1016/j.ecocom.2025.101116](https://doi.org/10.1016/j.ecocom.2025.101116)
  .
  <https://www.sciencedirect.com/science/article/pii/S1476945X25000017>

- Rao CR (1982) Diversity and dissimilarity coefficients: A unified
  approach. Theoretical Population Biology 21.
  [doi:10.1016/0040-5809(82)90004-1](https://doi.org/10.1016/0040-5809%2882%2990004-1)
  .

## See also

[`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
