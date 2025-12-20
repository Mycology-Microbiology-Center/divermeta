# Convert Matrix Indices to Distance Vector Position

Converts row and column indices from a full distance matrix to the
corresponding position in the compact vector representation used by R's
`dist` objects. This function maps `(i, j)` matrix coordinates to the
linear index `k` in the lower-triangular distance vector.

## Usage

``` r
convert_to_dist_indices(i, j, n)
```

## Arguments

- i:

  Integer vector of row indices (should be less than `j`)

- j:

  Integer vector of column indices (should be greater than `i`)

- n:

  Integer, the size of the square distance matrix (number of objects)

## Value

Integer vector of the same length as `i` and `j`, containing the
positions in the distance vector corresponding to the input matrix
indices.

## Details

R's `dist` objects store distance matrices in a compact form by only
storing the lower triangle (excluding the diagonal) as a vector. The
storage order follows:

- Row-major order of the lower triangle: (2,1), (3,1), (3,2), (4,1),
  (4,2), (4,3), ...

- For an n×n matrix, the vector has length n×(n-1)/2

- The formula used is: \\k = (i-1) \times n - \frac{(i-1) \times i}{2} +
  (j-i)\\ for \\i \< j\\

## Note

This function assumes that `i < j` for all input pairs. If `i > j`, the
result will be incorrect as it corresponds to the upper triangle, which
is not stored in `dist` objects. The function does not check for valid
indices.
