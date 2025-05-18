library(divermeta)

# Experiments
# ----------------------------

# Inventory
# --------

# High Multiplicity
# Three clusters of each ten different units
# Multiplicity should be 10
ab <- rep(10,30)
clust <- c(rep(1,10), rep(2,10), rep(3,10))
print(paste("Multiplicity: ",multiplicity.inventory(ab, clust)))

# Low Multiplicity
# Three clusters of each one different units
# Multiplicity should be 1
ab <- rep(10,3)
clust <- c(1,2,3)
print(paste("Multiplicity: ",multiplicity.inventory(ab, clust)))


# Functional
# ----------
# The idea of this experiments is to see why we needed a tweaked definition of
# functional diversity if we want to use ratio of the unclustered / clustersed
# as the definition of Multiplicity.

# The traditional functional diversity by Chao fails in some cases to capture
# what we want

# Basic
# ----------------------------
# n units with max distance and equal abundance
# functional diversity should be n
# Both work

sig <- 0.5
n <- 10
ab <- matrix(1, nrow = 1, ncol = n)
diss <- matrix(sig, nrow = n, ncol = n)
diag(diss) <- 0

print("Basic")
print(paste("Functional diversity (traditional):", diversity.functional.traditional(ab, diss)))
print(paste("Functional diversity", diversity.functional(ab, diss, sig)))
print("")
print("")
print("")


# Doubling property
# ----------------------------
# Builds an assemblage of three identical groups, inter Clustered distance are maximum.
# The diversity of the assemblage should be three times the individual
# diversity.
# Note: We use ratio to check for doubeling property, but this is not
# multiplicity since there is there is no actual clustering.
# Traditional is pretty close, the tweaked one we use works fine.
# Ratio should be 1
sig <- 0.9
n <- 10
ab_unit <- rbind(runif(n, 100, 1000))
diss_unit <- random_matrix <- matrix(runif(n * n, min = 0.3, max = 0.6), nrow = n, ncol = n)
diag(diss_unit) <- 0

div_trad <- diversity.functional.traditional(ab_unit, diss_unit)
div <- diversity.functional(ab_unit, diss_unit, sig)

# Assemblage
diss <- matrix(sig, nrow = 3 * n, ncol = 3 * n)
diss[1:n, 1:n] <- diss_unit
diss[(n + 1):(2 * n), (n + 1):(2 * n)] <- diss_unit
diss[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- diss_unit

ab <- cbind(ab_unit, ab_unit, ab_unit)


print("Doubleling Property")
print(paste("Unit - Functional diversity (traditional):", diversity.functional.traditional(ab_unit, diss_unit)))
print(paste("Assemblage  - Functional diversity (traditional):", diversity.functional.traditional(ab, diss)))
print(paste("Ratio:", diversity.functional.traditional(ab, diss) / diversity.functional.traditional(ab_unit, diss_unit)))
print("")
print(paste("Unit - Functional diversity", diversity.functional(ab_unit, diss_unit, sig)))
print(paste("Assemblage  - Functional diversity", diversity.functional(ab, diss, sig)))
print(paste("Ratio:", diversity.functional(ab, diss, sig) / diversity.functional(ab_unit, diss_unit, sig)))
print("")
print("")
print("")



# Assemblage  of units
# ----------------------------
# Diversity of an assemblage  of groups of identical units (distance zero
# between them) should be the same as an assemblage of representative units (a
# matrix of ones with zero on the diagonal)
# Only the tweaked one we use works, the traditional one fails.
# ratio (or multiplicity) should be 1
sig <- 1
n <- 10
ab_unit <- rbind(runif(n, 100, 1000))
diss_unit <- matrix(0, nrow = n, ncol = n)

div_trad <- diversity.functional.traditional(ab_unit, diss_unit)
div <- diversity.functional(ab_unit, diss_unit, sig)

# Assemblage
diss <- matrix(sig, nrow = 3 * n, ncol = 3 * n)
diss[1:n, 1:n] <- diss_unit
diss[(n + 1):(2 * n), (n + 1):(2 * n)] <- diss_unit
diss[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- diss_unit

ab <- cbind(ab_unit, ab_unit, ab_unit)

# Equivalent (Clusteres)
ab_clust <- cbind(sum(ab_unit), sum(ab_unit), sum(ab_unit))
diss_clust <- matrix(sig, ncol = 3, nrow = 3)
diag(diss_clust) <- 0

print("Assemblage of Identical units")
print(paste("Assemblage  - Functional diversity (traditional):", diversity.functional.traditional(ab, diss)))
print(paste("Representatives - Functional diversity (traditional):", diversity.functional.traditional(ab_clust, diss_clust)))
print(paste("Ratio:", diversity.functional.traditional(ab, diss) / diversity.functional.traditional(ab_clust, diss_clust)))
print("")
print(paste("Assemblage  - Functional diversity", diversity.functional(ab, diss, sig)))
print(paste("Representatives - Functional diversity", diversity.functional(ab_clust, diss_clust, sig)))
print(paste("Ratio (Multiplicity):",  diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig)))
print("")
print("")
print("")



# Assemblage  of units with some multiplicity
# Same as the previous experiment but instead of identical, they will have
# some variability, increasing the multiplicity.
# Move around some of the parameters to see build intuition

# Parameters
sig <- 0.5
n <- 10
min_abundance = 100
max_abundance = 1000
min_intra_distance = 0.1
max_intra_distance = 0.5



ab_unit <- rbind(runif(n, min_abundance, max_abundance))
diss_unit <- matrix(runif(n * n, min = min_intra_distance, max = max_intra_distance), nrow = n, ncol = n)
diag(diss_unit) <- 0


div_trad <- diversity.functional.traditional(ab_unit, diss_unit)
div <- diversity.functional(ab_unit, diss_unit, sig)

# Assemblage
diss <- matrix(sig, nrow = 3 * n, ncol = 3 * n)
diss[1:n, 1:n] <- diss_unit
diss[(n + 1):(2 * n), (n + 1):(2 * n)] <- diss_unit
diss[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- diss_unit

ab <- cbind(ab_unit, ab_unit, ab_unit)

# Equivalent (Clusteres)
ab_clust <- cbind(sum(ab_unit), sum(ab_unit), sum(ab_unit))
diss_clust <- matrix(sig, ncol = 3, nrow = 3)
diag(diss_clust) <- 0

print("Assemblage of Identical units")
print(paste("Assemblage  - Functional diversity (traditional):", diversity.functional.traditional(ab, diss)))
print(paste("Representatives - Functional diversity (traditional):", diversity.functional.traditional(ab_clust, diss_clust)))
print(paste("Ratio:", diversity.functional.traditional(ab, diss) / diversity.functional.traditional(ab_clust, diss_clust)))
print("")
print(paste("Assemblage  - Functional diversity", diversity.functional(ab, diss, sig)))
print(paste("Representatives - Functional diversity", diversity.functional(ab_clust, diss_clust, sig)))
print(paste("Ratio (Multiplicity):",  diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig)))
print("")
print("")
print("")

multiplicity.distance(ab,diss, ab_clust, diss_clust, sig)
#

