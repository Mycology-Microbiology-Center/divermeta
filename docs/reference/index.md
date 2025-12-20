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
- [`multiplicity.distance.by_blocks()`](https://mycology-microbiology-center.github.io/divermeta/reference/multiplicity.distance.by_blocks.md)
  : Distance-based multiplicity by blocks

## Diversity and redundancy

Diversity, redundancy, and related indices:

- [`metagenomic.alpha.index()`](https://mycology-microbiology-center.github.io/divermeta/reference/metagenomic.alpha.index.md)
  : Metagenomic Alpha-Diversity Index (MAD) (Finn 2024)
- [`diversity.functional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.md)
  : Distance-based functional diversity (q = 1) (Chiu & Chao 2014)
- [`diversity.functional.traditional()`](https://mycology-microbiology-center.github.io/divermeta/reference/diversity.functional.traditional.md)
  : Distance-based functional diversity (order q) (Chiu & Chao 2014)
- [`redundancy()`](https://mycology-microbiology-center.github.io/divermeta/reference/redundancy.md)
  : Functional redundancy (Re) (Ricotta & Pavoine 2025)
- [`raoQuadratic()`](https://mycology-microbiology-center.github.io/divermeta/reference/raoQuadratic.md)
  : Rao quadratic entropy (Q) (Rao 1982)

## Distance and support utilities

Helper functions:

- [`dist_quadratic_form()`](https://mycology-microbiology-center.github.io/divermeta/reference/dist_quadratic_form.md)
  : Quadratic Form for a Distance Object
- [`cluster_distance_matrix()`](https://mycology-microbiology-center.github.io/divermeta/reference/cluster_distance_matrix.md)
  : Compute Distance Matrix Between Clusters
- [`convert_to_dist_indices()`](https://mycology-microbiology-center.github.io/divermeta/reference/convert_to_dist_indices.md)
  : Convert Matrix Indices to Distance Vector Position

## Visualisation

Function for visualizing abundances and distance structures:

- [`visualize_abundance()`](https://mycology-microbiology-center.github.io/divermeta/reference/visualize_abundance.md)
  : Visualize an abundance matrix as a sized-tile plot
- [`visualize_dist()`](https://mycology-microbiology-center.github.io/divermeta/reference/visualize_dist.md)
  : Visualize a distance matrix as a sized-tile plot
