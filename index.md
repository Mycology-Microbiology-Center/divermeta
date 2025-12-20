# divermeta

**divermeta** is an R package for quantifying the diversity hidden
within clustered units in metagenomic and ecological analyses. It
implements the *multiplicity* index, a novel metric that measures the
average diversity within clusters and quantifies how much diversity is
lost when sequences or features are grouped into operational units
(e.g., OTUs, gene clusters, functional groups).

## What is multiplicity?

In metagenomic and ecological studies, sequences are typically clustered
into operational units to simplify analysis. While this clustering is
necessary and useful, it inevitably obscures the diversity that exists
*within* each cluster.

**Multiplicity** quantifies this hidden diversity. It answers the
question: “On average, how many distinct elements are contained within
each cluster?”

- **High multiplicity** means clusters contain many diverse elements
  (e.g., multiple gene variants or haplotypes)  
- **Low multiplicity** means clusters are dominated by similar elements
  (e.g., nearly identical sequences within an OTU)

Multiplicity complements traditional diversity indices by partitioning
total diversity into:

- Between-cluster diversity (what traditional indices measure)  
- Within-cluster diversity (what multiplicity measures)

This partitioning helps researchers understand how clustering choices
affect diversity estimates and can reveal important biological patterns
that would otherwise remain hidden.

## Installation

You can install the development version of divermeta from
[GitHub](https://github.com/Mycology-Microbiology-Center/divermeta)
with:

``` r
# install.packages("pak")
pak::pak("Mycology-Microbiology-Center/divermeta")
```

Or using `remotes`:

``` r
# install.packages("remotes")
remotes::install_github("Mycology-Microbiology-Center/divermeta")
```

## Quick start

The main function
[`divermeta()`](https://mycology-microbiology-center.github.io/divermeta/reference/divermeta.md)
computes multiplicity and other diversity indices across all samples.
Here’s a minimal example:

``` r
library(divermeta)

# Abundance matrix: rows = features, columns = samples
abund <- matrix(
  c(10, 5,  8, 12,  6, 3,  4, 9),
  nrow = 4,
  ncol = 2,
  dimnames = list(
    c("gene1", "gene2", "gene3", "gene4"),
    c("Sample_A", "Sample_B")
  )
)

# Assign features to clusters
clusters <- c(gene1 = "Group_A", gene2 = "Group_A", 
              gene3 = "Group_B", gene4 = "Group_B")

# Compute inventory multiplicity
divermeta(abund, clusters = clusters, indices = "multiplicity_inventory")
#>     Sample multiplicity_inventory
#> 1 Sample_A               1.929710
#> 2 Sample_B               1.868481
```

## Learn more

For detailed examples, use cases, and technical documentation:

- **[Vignette](https://mycology-microbiology-center.github.io/divermeta/articles/divermeta.html)**:
  Comprehensive guide with examples for both inventory and
  distance-based multiplicity
- **[Function
  reference](https://mycology-microbiology-center.github.io/divermeta/reference/index.html)**:
  Complete documentation of all functions and parameters
- **[Manuscript](https://github.com/Mycology-Microbiology-Center/divermeta)**:
  Theoretical background, mathematical foundations, and real-world
  applications

## Citation

If you use **divermeta** in your research, please cite:

> González-Casabianca F, Dulya O, Flores-Almaraz V, Tedersoo L,
> Mikryukov V. Multiplicity quantifies intra-cluster diversity in
> metagenomic analysis // in prep.
