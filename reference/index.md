# Package index

## Core function

Main function of the package for computing the multiplicity index:

- [`divermeta()`](https://mycology-microbiology-center.github.io/divermeta/reference/divermeta.md)
  : Estimate diversity and multiplicity indices across samples

## Multiplicity indices

Functions to compute individual multiplicity indices and related
quantities:

- [`multiplicity.inventory()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.inventory.md)
  : Inventory multiplicity
- [`multiplicity.distance()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.md)
  : Distance-based multiplicity
- [`relative.multiplicity()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.md)
  : Relative multiplicity
- [`relative.multiplicity.ref_div()`](https://mycology-microbiology-center.github.io/divermeta/reference/relative.multiplicity.ref_div.md)
  : Relative multiplicity from reference diversities

## Diversity and redundancy

Diversity, redundancy, and related indices:

- [`metagenomic.alpha.index()`](https://mycology-microbiology-center.github.io/divermeta/reference/metagenomic.alpha.index.md)
  : Metagenomic Alpha-Diversity Index (MAD) (Finn 2024)
- [`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)
  : Distance-based functional diversity (Chiu & Chao 2014)
- [`diversity.functional.traditional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.traditional.md)
  : Distance-based functional diversity (order q) (Chiu & Chao 2014)
- [`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)
  : Functional redundancy (Re) (Ricotta & Pavoine 2025)
- [`average.redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/average.redundancy.md)
  : Average redundancy
- [`proportional.abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.md)
  : Proportional abundance
- [`proportional.abundance.ref_abu()`](https://mycology-microbiology-center.github.io/divermeta/reference/proportional.abundance.ref_abu.md)
  : Proportional abundance with known reference abundances
- [`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
  : Rao quadratic entropy (Q) (Rao 1982)

## Distance and support utilities

Helper functions:

- [`unit_distances()`](https://mycology-microbiology-center.github.io/divermeta/reference/unit_distances.md)
  : Distances between units

## Large distance files

Reading the distances from files in chunks, and checking those files:

- [`distance-files`](https://mycology-microbiology-center.github.io/divermeta/reference/distance-files.md)
  : Distances from files
- [`check_distance_file()`](https://mycology-microbiology-center.github.io/divermeta/reference/check_distance_file.md)
  : Check a large distance file

## Visualisation

Function for visualizing abundances and distance structures:

- [`visualize_abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/visualize_abundance.md)
  : Visualize an abundance table as a sized-tile plot
- [`visualize_dist()`](https://mycology-microbiology-center.github.io/divermeta/reference/visualize_dist.md)
  : Visualize distances as a sized-tile plot
