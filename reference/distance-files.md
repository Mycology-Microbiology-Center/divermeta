# Distances from files

Every index that uses distances accepts `diss` either as a three-column
table in memory or as the path to a file (or several files, read one
after the other) with that table. A file is read in chunks of
`chunk_size` rows, so the table never has to fit in memory: memory
depends on the number of samples, subunits and units, and on
`chunk_size`, not on the number of pairs.

## File format

- Delimited text with at least three columns, taken by position: the
  identifiers of two subunits and the distance between them. Further
  columns are ignored.

- Files ending in `.csv` are comma separated and files ending in `.tsv`
  are tab separated. For other extensions the separator (tab, comma or
  white space) is detected from the first line.

- With `header = TRUE` (default) the first line of every file is a
  header and is skipped; if its third field is a number, the line looks
  like data and the index stops, so a row is never dropped silently.
  Files without a header need `header = FALSE`, and then a first row
  with a non-numeric distance is an error like any other row.

- Identifiers are read as text and matched as text against the column
  names of `abund` (so `100000` in a file does not match a column named
  `1e+05`).

- Files may be plain or compressed with gzip (`.gz`), bzip2 (`.bz2`) or
  xz (`.xz`). Files compressed with zstd (`.zst`) are read through the
  `zstd` command line tool, which must be installed.

## Every pair exactly once

Pairs are not kept in memory, so they cannot be de-duplicated while
reading: a file must list every pair of subunits the index needs
**exactly once** (in either orientation). Rows with identifiers not in
`abund`, self pairs and pairs the index does not use (e.g. pairs of
subunits of different units for indices that only use the pairs inside
units) are ignored, also when their distance is NA, infinite or
negative. The exception is
[`metagenomic.alpha.index()`](https://mycology-microbiology-center.github.io/divermeta/reference/metagenomic.alpha.index.md),
which checks the distance of every row with two known identifiers:
there, an NA, infinite or negative distance is an error even in a pair
it does not use.

Pairs with a subunit that has zero abundance in every sample are not
needed (except for the homogeneous reference of
[`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md),
also when requested through
[`divermeta()`](https://mycology-microbiology-center.github.io/divermeta/reference/divermeta.md)):
if some of them are missing, a warning is given and the index is
computed without them. When they are listed, they count as any other
row. When `diss` is a file:

- The file is always checked to exist, to be readable and to have three
  columns, and every chunk is checked for NA, infinite or negative
  distances in the pairs the index uses. Every distance must still be
  readable as a number.

- The number of rows of the needed pairs is compared with the number of
  pairs: too few or too many rows is an error.

- With `check_large_distance_file = FALSE` (default), nothing else is
  checked and a warning says so. A missing pair would count as distance
  0 (as if the two subunits were identical) and a duplicated pair twice,
  and a missing pair and a duplicated one cancel out in the row count.

- With `check_large_distance_file = TRUE`,
  [`check_distance_file()`](https://mycology-microbiology-center.github.io/divermeta/reference/check_distance_file.md)
  is run first, which finds every missing, duplicated, conflicting or
  invalid row, and the index stops if there is any. The check reads the
  file once more and needs about 16 bytes of temporary disk space per
  pair.

In-memory tables are always fully checked, and there a pair listed twice
with the same distance is accepted.

## See also

[`check_distance_file()`](https://mycology-microbiology-center.github.io/divermeta/reference/check_distance_file.md)

## Examples

``` r
abund <- matrix(c(2, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
                dimnames = list(c("S1", "S2"), c("a", "b", "c")))
diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))

diss_file <- tempfile(fileext = ".csv")
write.csv(diss, diss_file, row.names = FALSE)

raoQuadratic(abund, diss_file, check_large_distance_file = TRUE)
#>        S1        S2 
#> 0.1777778 0.4222222 

# A file without a header
write.table(diss, diss_file, sep = ",", row.names = FALSE, col.names = FALSE)
raoQuadratic(abund, diss_file, check_large_distance_file = TRUE, header = FALSE)
#>        S1        S2 
#> 0.1777778 0.4222222 
```
