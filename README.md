
<!-- README.md is generated from README.Rmd. Please edit that file -->

# divermeta

<!-- badges: start -->
<!-- badges: end -->

The goal of divermeta is to provide the users / researchers with
functions to compute both taxonomic and functional multiplicities.

## Installation

You can install the development version of divermeta from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Mycology-Microbiology-Center/divermeta")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(divermeta)
## basic example code

# Taxonomic
# --------

# High Multiplicity
# Three clusters of each ten different units
# Multiplicity should be 10
ab <- rep(10,30)
clust <- c(rep(1,10), rep(2,10), rep(3,10))
print(paste("Multiplicity: ",multiplicity.taxonomic(ab, clust)))
#> [1] "Multiplicity:  10"

# Low Multiplicity
# Three clusters of each one different units
# Multiplicity should be 1
ab <- rep(10,3)
clust <- c(1,2,3)
print(paste("Multiplicity: ",multiplicity.taxonomic(ab, clust)))
#> [1] "Multiplicity:  1"
```
