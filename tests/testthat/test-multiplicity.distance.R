

# Tests for distance-based multiplicity

# Basic
# ----------------------------
# n units with max distance and equal abundance
# functional diversity should be n
# Both work

test_that("Basic", {

  sig <- 0.5
  n <- 10
  ab <- rep(1, n)
  diss <- matrix(sig, nrow = n, ncol = n)
  diag(diss) <- 0

  expect_equal(round(diversity.functional.traditional(ab, diss), 3), n)
  expect_equal(diversity.functional(ab, diss, sig), n)

})


# Deterministic numeric checks (two-species and sigma capping)
# -----------------------------------------------------------
test_that("Two-species numeric and sigma capping", {

  # Two species with unequal abundances, distance above sigma
  ab <- c(2, 1)
  d <- 0.7
  sig <- 0.5

  diss <- matrix(c(0, d, d, 0), nrow = 2, ncol = 2)

  # diversity.functional must cap to sigma; equivalent to using min(d, sig)
  diss_capped <- matrix(c(0, sig, sig, 0), nrow = 2, ncol = 2)
  expect_equal(diversity.functional(ab, diss, sig), diversity.functional(ab, diss_capped, sig))

  # multiplicity.distance closed form for 2 species
  # raoQ_before = 2 * p1 * p2 * min(d, sig); raoQ_after with two clusters at distance sig: 2 * p1 * p2 * sig
  p <- ab / sum(ab)
  Q_before <- 2 * p[1] * p[2] * sig
  Q_after <- 2 * p[1] * p[2] * sig

  # After clustering matrix uses inter-cluster distance = sig
  ab_clust <- ab
  diss_clust <- diss_capped

  expect_equal(
    multiplicity.distance(ab, diss, ab_clust, diss_clust, sig),
    (sig - Q_after) / (sig - Q_before)
  )
})



# Doubling property
# ----------------------------
# Builds an assemblage of three identical groups, inter Clustered distance are maximum.
# The diversity of the assemblage should be three times the individual
# diversity.
# Note: We use ratio to check for doubling property, but this is not
# multiplicity since there is there is no actual clustering.
# Traditional is pretty close, the tweaked one we use works fine.
# Ratio should be 3
test_that("Doubling Property", {

  set.seed(42)
  sig <- 0.9
  n <- 10
  ab_unit <- runif(n, 100, 1000)
  diss_unit <- random_matrix <- matrix(runif(n * n, min = 0.3, max = 0.6), nrow = n, ncol = n)
  diag(diss_unit) <- 0

  div_trad <- diversity.functional.traditional(ab_unit, diss_unit)
  div <- diversity.functional(ab_unit, diss_unit, sig)

  # Assemblage
  diss <- matrix(sig, nrow = 3 * n, ncol = 3 * n)
  diss[1:n, 1:n] <- diss_unit
  diss[(n + 1):(2 * n), (n + 1):(2 * n)] <- diss_unit
  diss[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- diss_unit

  ab <- c(ab_unit, ab_unit, ab_unit)

  expect_equal(round(diversity.functional.traditional(ab, diss) / diversity.functional.traditional(ab_unit, diss_unit)), 3)
  expect_equal(diversity.functional(ab, diss, sig) / diversity.functional(ab_unit, diss_unit, sig), 3)

})





# Assemblage  of units
# ----------------------------
# Diversity of an assemblage  of groups of identical units (distance zero
# between them) should be the same as an assemblage of representative units (a
# matrix of ones with zero on the diagonal)
# Only the tweaked one we use works, the traditional one fails.
# ratio (or multiplicity) should be 1
test_that("Assemblage  of units", {

    sig <- 1
    n <- 10
    ab_unit <- runif(n, 100, 1000)
    diss_unit <- matrix(0, nrow = n, ncol = n)

    div_trad <- diversity.functional.traditional(ab_unit, diss_unit)
    div <- diversity.functional(ab_unit, diss_unit, sig)

    # Assemblage
    diss <- matrix(sig, nrow = 3 * n, ncol = 3 * n)
    diss[1:n, 1:n] <- diss_unit
    diss[(n + 1):(2 * n), (n + 1):(2 * n)] <- diss_unit
    diss[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- diss_unit

    ab <- c(ab_unit, ab_unit, ab_unit)

    # Equivalent (Clusteres)
    ab_clust <- c(sum(ab_unit), sum(ab_unit), sum(ab_unit))
    diss_clust <- matrix(sig, ncol = 3, nrow = 3)
    diag(diss_clust) <- 0

    expect_false((diversity.functional.traditional(ab, diss) / diversity.functional.traditional(ab_clust, diss_clust)) == 1)
    expect_equal(diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig), 1)

})


# Ratio Equivalence. Distance based multiplicity should be equivalent to the
# ratio of the functional diversities
test_that("Ratio equivalence", {

  set.seed(42)
  for(i in 1:10)
  {
    # Parameters
    sig <- runif(1)
    n <- 10
    min_abundance = 100
    max_abundance = 1000
    min_intra_distance = 0.1
    max_intra_distance = 0.5



    ab_unit <- runif(n, min_abundance, max_abundance)
    diss_unit <- matrix(runif(n * n, min = min_intra_distance, max = max_intra_distance), nrow = n, ncol = n)
    diag(diss_unit) <- 0


    div_trad <- diversity.functional.traditional(ab_unit, diss_unit)
    div <- diversity.functional(ab_unit, diss_unit, sig)

    # Assemblage
    diss <- matrix(sig, nrow = 3 * n, ncol = 3 * n)
    diss[1:n, 1:n] <- diss_unit
    diss[(n + 1):(2 * n), (n + 1):(2 * n)] <- diss_unit
    diss[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- diss_unit

    ab <- c(ab_unit, ab_unit, ab_unit)

    # Equivalent (Clusteres)
    ab_clust <- c(sum(ab_unit), sum(ab_unit), sum(ab_unit))
    diss_clust <- matrix(sig, ncol = 3, nrow = 3)
    diag(diss_clust) <- 0

    ratio <- diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig)
    m <- multiplicity.distance(ab,diss, ab_clust, diss_clust, sig)

    expect_equal(round(ratio,5), round(m,5))

  }

})

# Deterministic by-blocks equals manual multiplicity
# --------------------------------------------------
test_that("by_blocks equals manual for simple case", {

  # Four elements in two clusters (1-2, 3-4)
  ids <- c("a", "b", "c", "d")
  ab <- c(2, 3, 5, 7)
  clust <- c(1, 1, 2, 2)
  sig <- 0.8

  # Within-cluster distances (symmetric); between clusters implicitly = sig for by_blocks
  # Build diss_frame with unique pairs only
  # Only include within-cluster unique pairs; cross-cluster assumed to be sigma internally
  df <- data.frame(
    ID1 = c("b", "d"),
    ID2 = c("a", "c"),
    Distance = c(0.2, 0.15),
    stringsAsFactors = FALSE
  )

  # by_blocks computation
  mb <- multiplicity.distance.by_blocks(ids = ids, ab = ab, diss_frame = df, clust = clust, sigma = sig)

  # Manual full matrices: within clusters as above, between clusters = sig
  diss <- matrix(sig, nrow = 4, ncol = 4, dimnames = list(ids, ids))
  diag(diss) <- 0
  diss["a", "b"] <- diss["b", "a"] <- 0.2
  diss["c", "d"] <- diss["d", "c"] <- 0.15

  ab_clust <- tapply(ab, clust, sum)
  diss_clust <- matrix(sig, nrow = 2, ncol = 2)
  diag(diss_clust) <- 0

  mm <- multiplicity.distance(ab, diss, ab_clust, diss_clust, sig)

  expect_equal(mb, mm, tolerance = 1e-12)
})

# Block equivalence
# Checks that the by blocks implementation gives the same results as the
# normal implementation
test_that("Implementation equivalence",{

  set.seed(42)

  for(w_ in 1:10)
  {
      total <- 10
      sigma  <- 0.3

      # Example matrices
      A <- matrix(runif(total**2, 0.1, sigma), ncol = total, nrow = total)
      B <- matrix(runif(total**2, 0.1, sigma), ncol = total, nrow = total)
      C <- matrix(runif(total**2, 0.1, sigma), ncol = total, nrow = total)
      D <- matrix(runif(total**2, 0.1, sigma), ncol = total, nrow = total)


      A <- (A+t(A))/2
      B <- (B+t(B))/2
      C <- (C+t(C))/2
      D <- (D+t(D))/2


      # Assign row and column names
      dfs <- c()
      i <- 0
      clust <- c()
      blocks <- list(A,B,C,D)
      for(M in blocks)
      {

        i <- i + 1
        clust <- c(clust,rep(i, total))
        rownames(M) <- colnames(M) <- (1+(i-1)*total):(i*total)
        df_M <- as.data.frame(as.table(M))
        colnames(df_M) <- c("ID1","ID2","Distance")
        df_M[["ID1"]] <- as.numeric(df_M[["ID1"]]) + (i-1)*total
        df_M[["ID2"]] <- as.numeric(df_M[["ID2"]]) + (i-1)*total
        df_M <- df_M[df_M$ID1 > df_M$ID2,]
        dfs[[i]] <- df_M

      }

      ids <- 1:(total*length(blocks))
      ab <- runif(total*length(blocks), 1, 100)
      ab_clust <- tapply(ab, clust, sum)
      diss_frame <- do.call(rbind, dfs)



      # Create block diagonal matrix
      diss <- matrix(sigma, nrow = total * length(blocks), ncol = total * length(blocks))

      offset <- 0
      for (M in blocks) {
        diss[(1:total) + offset, (1:total) + offset] <- M
        offset <- offset + total
      }

      diag(diss) <- 0


      diss_clust <- matrix(sigma, ncol = length(blocks), nrow = length(blocks))
      diag(diss_clust) <- 0

      byBlocks <- multiplicity.distance.by_blocks(ids, ab, diss_frame, clust, sigma)
      classic <- multiplicity.distance(ab, diss, ab_clust, diss_clust, sigma)

      expect_equal(byBlocks - classic, 0)
  }
})

