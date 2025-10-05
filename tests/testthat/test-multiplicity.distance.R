

# Distance Based

# Basic
# ----------------------------
# n units with max distance and equal abundance
# functional diversity should be n
# Both work

test_that("Basic", {

  sig <- 0.5
  n <- 10
  ab <- matrix(1, nrow = 1, ncol = n)
  diss <- matrix(sig, nrow = n, ncol = n)
  diag(diss) <- 0

  expect_equal(round(diversity.functional.traditional(ab, diss), 3), n)
  expect_equal(diversity.functional(ab, diss, sig), n)

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

    expect_false((diversity.functional.traditional(ab, diss) / diversity.functional.traditional(ab_clust, diss_clust)) == 1)
    expect_equal(diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig), 1)

})


# Ratio Equivalence. Distance based multiplicity should be equivalent to the
# ratio of the functional diversities
test_that("Ratio equivalence", {
#
  for(i in 1:10)
  {
    # Parameters
    sig <- runif(1)
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

    ratio <- diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig)
    m <- multiplicity.distance(ab,diss, ab_clust, diss_clust, sig)

    expect_equal(round(ratio,5), round(m,5))

  }

})

# Block equivalence
# Checks that the by blocks implementation gives the same results as the
# normal implementation
test_that("Implementation equivalence",{

  set.seed(123)

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

