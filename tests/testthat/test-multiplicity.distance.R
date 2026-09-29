# Tests for distance-based multiplicity

# Deterministic numeric checks (two-species and sigma capping)
# -----------------------------------------------------------
test_that("Two-species case: sigma capping works correctly", {
  # Two species with unequal abundances, distance above sigma
  ab <- one_row(c(2, 1))
  d <- 0.7
  sig <- 0.5

  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = d)

  # diversity.functional must cap to sigma; equivalent to using min(d, sig)
  diss_capped <- data.frame(ID1 = "1", ID2 = "2", Distance = sig)
  expect_equal(diversity.functional(ab, diss, sig), diversity.functional(ab, diss_capped, sig))

  # multiplicity.distance closed form for 2 species
  # raoQ_before = 2 * p1 * p2 * min(d, sig); raoQ_after with two clusters at distance sig: 2 * p1 * p2 * sig
  p <- c(2, 1) / 3
  Q_before <- 2 * p[1] * p[2] * sig
  Q_after <- 2 * p[1] * p[2] * sig

  for (method in c("sigma", "min", "max", "average")) {
    expect_equal(
      unname(multiplicity.distance(ab, diss, clust = c(1, 2), method = method, sig = sig)),
      (sig - Q_after) / (sig - Q_before)
    )
  }
})



# Assemblage of units
# ----------------------------
# Tests that clustering copies of identical units (distance zero within clusters)
# does not change functional diversity. The diversity of an assemblage of groups
# of identical units should equal the diversity of representative units.
# Multiplicity should be 1 (no diversity lost).
test_that("Clustering identical units: multiplicity equals 1", {
  set.seed(42)
  sig <- 1
  n <- 10
  ab_unit <- runif(n, 100, 1000)
  diss_unit <- matrix(0, nrow = n, ncol = n)

  # Assemblage
  diss <- mat_to_table(block_matrix(list(diss_unit, diss_unit, diss_unit), sig))
  ab <- one_row(c(ab_unit, ab_unit, ab_unit))
  clust <- rep(1:3, each = n)

  # Equivalent (Clusteres)
  ab_clust <- one_row(c(sum(ab_unit), sum(ab_unit), sum(ab_unit)))
  diss_clust <- mat_to_table(max_diss(3, sig))


  expect_equal(unname(diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig)), 1)
  expect_equal(unname(multiplicity.distance(ab, diss, clust, sig = sig)), 1)
})


# Ratio equivalence
# ----------------------------
# Tests that distance-based multiplicity equals the ratio of functional diversities
# before and after clustering (delta D_sigma before / delta D_sigma after).
test_that("Distance-based multiplicity equals ratio of functional diversities", {
  set.seed(42)
  for (i in 1:10)
  {
    # Parameters
    sig <- runif(1)
    n <- 10
    min_abundance <- 100
    max_abundance <- 1000
    min_intra_distance <- 0.1
    max_intra_distance <- 0.5

    clust <- c(rep(1, n), rep(2, n), rep(3, n))

    ab_unit <- runif(n, min_abundance, max_abundance)
    diss_unit <- rand_diss(n, min = min_intra_distance, max = max_intra_distance)

    # Assemblage
    diss <- mat_to_table(block_matrix(list(diss_unit, diss_unit, diss_unit), sig))
    ab <- one_row(c(ab_unit, ab_unit, ab_unit))

    # Equivalent (Clusteres)
    ab_clust <- one_row(c(sum(ab_unit), sum(ab_unit), sum(ab_unit)))
    diss_clust <- mat_to_table(max_diss(3, sig))

    ratio <- diversity.functional(ab, diss, sig) / diversity.functional(ab_clust, diss_clust, sig)
    m <- multiplicity.distance(ab, diss, clust = clust, method = "sigma", sig = sig)

    expect_equal(round(unname(ratio), 5), round(unname(m), 5))
  }
})

# Implementation equivalence
# --------------------------------------------------
# Tests that the sigma method, which only needs the distances inside clusters,
# equals the legacy implementation with full distance matrices.
test_that("sigma method with pairs inside clusters equals legacy implementation", {
  # Four elements in two clusters (1-2, 3-4)
  ids <- c("a", "b", "c", "d")
  ab <- c(2, 3, 5, 7)
  clust <- c(1, 1, 2, 2)
  sig <- 0.8

  # Within-cluster distances only
  df <- data.frame(
    ID1 = c("b", "d"),
    ID2 = c("a", "c"),
    Distance = c(0.2, 0.15),
    stringsAsFactors = FALSE
  )

  mb <- multiplicity.distance(one_row(ab, ids), df, clust, sig = sig)

  # Manual full matrices: within clusters as above, between clusters = sig
  diss <- matrix(sig, nrow = 4, ncol = 4, dimnames = list(ids, ids))
  diag(diss) <- 0
  diss["a", "b"] <- diss["b", "a"] <- 0.2
  diss["c", "d"] <- diss["d", "c"] <- 0.15

  mm <- legacy$multiplicity.distance(ab = ab, diss = diss, clust = clust, method = "sigma", sig = sig)

  expect_equal(unname(mb), mm, tolerance = 1e-12)

  # Listing the pairs between clusters too does not change it
  expect_equal(multiplicity.distance(one_row(ab, ids), mat_to_table(diss), clust, sig = sig), mb)

  # A missing pair inside a cluster is an error
  expect_error(
    multiplicity.distance(one_row(ab, ids), df[1, ], clust, sig = sig),
    "Missing distances.*unit: 2"
  )
})

# Distances between clusters
# --------------------------------------------------
# Tests that the sigma method sets the distances between elements of different
# clusters to sigma, even when they are listed in the table.
test_that("sigma method sets distances between clusters to sigma", {
  ids <- c("a", "b", "c", "d")
  ab <- one_row(c(2, 3, 5, 7), ids)
  clust <- c(1, 1, 2, 2)
  sig <- 0.8

  df_within <- data.frame(
    ID1 = c("a", "c"),
    ID2 = c("b", "d"),
    Distance = c(0.2, 0.15),
    stringsAsFactors = FALSE
  )
  # Same pairs plus two pairs across clusters closer than sigma
  df_cross <- rbind(
    df_within,
    data.frame(ID1 = c("a", "b"), ID2 = c("c", "d"), Distance = c(0.1, 0.3))
  )

  mb_within <- multiplicity.distance(ab, df_within, clust, sig = sig)
  mb_cross <- multiplicity.distance(ab, df_cross, clust, sig = sig)
  expect_equal(mb_cross, mb_within)

  # Full matrix with every distance between clusters equal to sigma
  diss <- matrix(sig, nrow = 4, ncol = 4, dimnames = list(ids, ids))
  diag(diss) <- 0
  diss["a", "b"] <- diss["b", "a"] <- 0.2
  diss["c", "d"] <- diss["d", "c"] <- 0.15
  expect_equal(
    unname(mb_cross),
    legacy$multiplicity.distance(c(2, 3, 5, 7), diss, clust, method = "sigma", sig = sig),
    tolerance = 1e-12
  )

  # The legacy implementation uses the distances between clusters: different value
  diss["a", "c"] <- diss["c", "a"] <- 0.1
  diss["b", "d"] <- diss["d", "b"] <- 0.3
  expect_false(isTRUE(all.equal(
    unname(mb_cross),
    legacy$multiplicity.distance(c(2, 3, 5, 7), diss, clust, method = "sigma", sig = sig)
  )))

  # A pair across clusters listed in both orientations is ignored, not rejected
  df_twice <- rbind(df_cross, data.frame(ID1 = "c", ID2 = "a", Distance = 0.1))
  expect_equal(
    multiplicity.distance(ab, df_twice, clust, sig = sig),
    mb_within
  )
})

# Block equivalence for multiple clusters
# ----------------------------------------
# Tests that the sigma method produces identical results to the legacy
# implementation for complex multi-cluster scenarios.
test_that("sigma method equals legacy implementation for multiple clusters", {
  set.seed(42)

  for (w_ in 1:10)
  {
    total <- 10
    sigma <- 0.3

    # Example matrices
    blocks <- lapply(1:4, function(i) rand_diss(total, 0.1, sigma))

    dfs <- list()
    clust <- c()
    for (i in seq_along(blocks))
    {
      clust <- c(clust, rep(i, total))
      dfs[[i]] <- mat_to_table(blocks[[i]], (1 + (i - 1) * total):(i * total))
    }

    ids <- 1:(total * length(blocks))
    abund <- matrix(
      runif(3 * total * length(blocks), 1, 100),
      nrow = 3,
      dimnames = list(NULL, ids)
    )
    diss_frame <- do.call(rbind, dfs)

    # Create block diagonal matrix
    diss <- block_matrix(blocks, sigma)

    diss_clust <- max_diss(length(blocks), sigma)

    new <- multiplicity.distance(abund, diss_frame, clust, sig = sigma)
    classic <- apply(abund, 1, function(ab) {
      legacy$multiplicity.distance(ab = ab, diss = diss, clust = clust, method = "custom", sig = sigma, clust_ids_order = c(1, 2, 3, 4), diss_clust = diss_clust)
    })

    expect_equal(unname(new), unname(classic))
  }
})


test_that("Distance-based multiplicity is the same using a specific link method than invoking it through custom",{
  set.seed(42)
  for (i in seq_len(10))
  {
    n <- 25
    ids <- paste0("f", seq_len(n))
    abund <- matrix(runif(4 * n, min = 1, max = 100), nrow = 4, dimnames = list(NULL, ids))
    clust <- stats::setNames(sample(c("e", "d", "c", "b", "a"), n, replace = TRUE), ids)
    diss <- rand_diss(n, min = 0.3, max = 1)
    tab <- mat_to_table(diss, ids)
    sig <- 0.7

    for(method in c("min", "average", "max"))
    {

      diss_clust <- unit_distances(tab, clust, method = method, sig = sig)
      diss_clust <- diss_clust[sample(nrow(diss_clust)), ]

      m1 <- multiplicity.distance(abund, tab, clust, method = method, sig = sig)
      m2 <- multiplicity.distance(abund, tab, clust, method = "custom", sig = sig, diss_clust = diss_clust)

      expect_equal(m1, m2)

    }

    # The sigma method also sets the distances between subunits of different
    # clusters to sigma, while custom uses them: they agree only when those
    # distances are at least sigma
    diss_clust <- unit_distances(tab, clust, method = "sigma", sig = sig)
    m_sigma <- multiplicity.distance(abund, tab, clust, method = "sigma", sig = sig)
    expect_false(isTRUE(all.equal(
      m_sigma,
      multiplicity.distance(abund, tab, clust, method = "custom", sig = sig, diss_clust = diss_clust)
    )))

    far <- diss
    far[outer(clust, clust, "!=")] <- sig + 0.1
    tab_far <- mat_to_table(far, ids)
    expect_equal(
      multiplicity.distance(abund, tab_far, clust, method = "sigma", sig = sig),
      multiplicity.distance(abund, tab_far, clust, method = "custom", sig = sig, diss_clust = diss_clust)
    )
    expect_equal(multiplicity.distance(abund, tab_far, clust, method = "sigma", sig = sig), m_sigma)
}

})


test_that("min, max and average methods equal the legacy implementation", {
  set.seed(42)
  for (i in seq_len(5)) {
    n <- 30
    ids <- paste0("f", seq_len(n))
    # Character cluster ids (legacy fails with numeric ids that are not 1..k)
    clust <- sample(c("e", "d", "c", "b", "a"), n, replace = TRUE)
    diss <- rand_diss(n, min = 0.1, max = 1.2)
    abund <- matrix(
      runif(6 * n, min = 1, max = 100) * (runif(6 * n) > 0.3),
      nrow = 6,
      dimnames = list(NULL, ids)
    )
    tab <- mat_to_table(diss, ids)

    for (sig in c(0.4, 0.8, 1.5)) {
      for (method in c("min", "max", "average")) {
        new <- multiplicity.distance(abund, tab, clust, method = method, sig = sig)
        # Legacy is given a dist object: it gives wrong cluster distances for full matrices
        classic <- apply(abund, 1, function(ab) {
          legacy$multiplicity.distance(ab, stats::as.dist(diss), clust, method = method, sig = sig)
        })
        expect_equal(unname(new), unname(classic))
        expect_equal(
          multiplicity.distance(Matrix::Matrix(abund, sparse = TRUE), tab, clust, method = method, sig = sig),
          new
        )
      }
    }
  }
})


test_that("Custom cluster distances equal the legacy implementation", {
  set.seed(7)
  n <- 20
  ids <- paste0("f", seq_len(n))
  clust <- sample(rep(c("A", "B", "C", "D"), length.out = n))
  diss <- rand_diss(n, min = 0.1, max = 1)
  abund <- matrix(runif(5 * n, 1, 100), nrow = 5, dimnames = list(NULL, ids))
  clust_ids <- c("A", "B", "C", "D")
  diss_clust <- rand_diss(4, min = 0.2, max = 1.2)

  new <- multiplicity.distance(
    abund, mat_to_table(diss, ids), clust,
    method = "custom", sig = 0.9, diss_clust = mat_to_table(diss_clust, clust_ids)
  )
  classic <- apply(abund, 1, function(ab) {
    legacy$multiplicity.distance(
      ab, diss, clust,
      method = "custom", sig = 0.9, clust_ids_order = clust_ids, diss_clust = diss_clust
    )
  })
  expect_equal(unname(new), unname(classic))
})


test_that("Distance-based multiplicity input validation", {
  ids <- c("a", "b", "c", "d")
  abund <- one_row(c(2, 3, 5, 7), ids)
  clust <- c(1, 1, 2, 2)
  diss <- matrix(0.9, 4, 4)
  diag(diss) <- 0
  tab <- mat_to_table(diss, ids)

  expect_error(multiplicity.distance(abund, tab, clust, sig = 0), "positive")
  expect_error(multiplicity.distance(abund, tab, clust, sig = c(1, 2)), "positive")
  expect_error(multiplicity.distance(abund, tab, clust, method = "median"), "not supported")
  expect_error(multiplicity.distance(abund, tab, clust, method = "custom"), "diss_clust cannot be NULL")
  expect_error(multiplicity.distance(abund, tab, c(1, 1, 2)), "same number of subunits")

  # Methods other than sigma need every pair
  expect_error(multiplicity.distance(abund, tab[-2, ], clust, method = "average"), "Missing distances")
  expect_error(multiplicity.distance(abund, tab[-2, ], clust, method = "custom", diss_clust = data.frame(1, 2, 0.5)), "Missing distances")

  # Custom cluster distances need every pair of clusters
  expect_error(
    multiplicity.distance(one_row(1:3, c("a", "b", "c")), mat_to_table(max_diss(3), c("a", "b", "c")), c(1, 2, 3),
      method = "custom", diss_clust = data.frame(1, 2, 0.5)
    ),
    "Missing distances in `diss_clust`"
  )

  # Empty samples give NA
  expect_true(is.na(multiplicity.distance(one_row(c(0, 0, 0, 0), ids), tab, clust)))
})


test_that("Test cluster distances (unit_distances)", {
    set.seed(42)
    for (n_clust in 2:10) {
        # Constructs Matrix
        diss_clust <- matrix(
            runif(min = 0.3, max = 1, n_clust * n_clust),
            nrow = n_clust,
            ncol = n_clust
        )

        # Converts to distance
        diss_clust <- (diss_clust + t(diss_clust)) / 2
        diag(diss_clust) <- 0
        diss_clust <- round(diss_clust, 2)

        clust_ids_order <- seq_len(n_clust)

        # Generates the cluster ids
        clust <- c()
        for (i in clust_ids_order) {
            clust <- c(clust, rep(i, runif(1, min = 1, max = 10)))
        }

        n_elem <- length(clust)
        n_clust <- length(clust_ids_order)
        ids <- paste0("e", seq_len(n_elem))
        named_clust <- stats::setNames(clust, ids)

        # Generates the distance
        diss_avg <- matrix(0, nrow = n_elem, ncol = n_elem)
        diss_max <- matrix(0, nrow = n_elem, ncol = n_elem)
        diss_min <- matrix(2, nrow = n_elem, ncol = n_elem)
        diag(diss_min) <- 0

        for (i in seq_len(n_clust - 1)) {
            for (j in (i + 1):n_clust) {
                # Average
                diss_avg[which(i == clust), which(j == clust)] <- diss_clust[
                    i,
                    j
                ]
                diss_avg[which(j == clust), which(i == clust)] <- diss_clust[
                    i,
                    j
                ]

                # Distorts
                if (sum(i == clust) > 1 && sum(j == clust) > 1) {
                    noise <- runif(1, min = 0.1, max = 0.3)
                    diss_avg[
                        which(i == clust)[1],
                        which(j == clust)[1]
                    ] <- diss_clust[i, j] - noise
                    diss_avg[
                        which(j == clust)[1],
                        which(i == clust)[1]
                    ] <- diss_clust[i, j] - noise
                    diss_avg[
                        which(i == clust)[2],
                        which(j == clust)[2]
                    ] <- diss_clust[i, j] + noise
                    diss_avg[
                        which(j == clust)[2],
                        which(i == clust)[2]
                    ] <- diss_clust[i, j] + noise
                }

                # Max
                diss_max[
                    which(i == clust)[1],
                    which(j == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]
                diss_max[
                    which(j == clust)[1],
                    which(i == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]

                # Min
                diss_min[
                    which(i == clust)[1],
                    which(j == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]
                diss_min[
                    which(j == clust)[1],
                    which(i == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]
            }
        }

        expected <- mat_to_table(diss_clust, as.character(clust_ids_order))
        as_matrix <- function(tab) {
            m <- matrix(0, n_clust, n_clust)
            m[cbind(as.integer(tab$ID1), as.integer(tab$ID2))] <- tab$Distance
            m[cbind(as.integer(tab$ID2), as.integer(tab$ID1))] <- tab$Distance
            m
        }

        for (method in c("average", "min", "max")) {
            diss <- switch(method, average = diss_avg, min = diss_min, max = diss_max)
            res <- unit_distances(mat_to_table(diss, ids), named_clust, method = method)
            expect_equal(as_matrix(res), diss_clust)

            # Same as the legacy implementation (given a dist object)
            expect_equal(
                as_matrix(res),
                unname(as.matrix(legacy$cluster_distance_matrix(stats::as.dist(diss), clust, clust_ids_order, method = method)))
            )
        }

        res <- unit_distances(mat_to_table(diss_avg, ids), named_clust, method = "sigma", sig = 0.7)
        expect_equal(as_matrix(res), max_diss(n_clust, 0.7))
    }

    # Validation
    clust <- c(a = 1, b = 1, c = 2)
    tab <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.1, 0.5, 0.7))
    expect_error(unit_distances(tab, unname(clust)), "named")
    expect_error(unit_distances(tab[1:2, ], clust, method = "average"), "Missing distances")
    expect_error(unit_distances(tab, clust, method = "median"), "not supported")
    expect_equal(unit_distances(tab, clust, method = "average")$Distance, 0.6)
})


test_that("sigma method caps the distances inside clusters (closed form)", {
  # Two clusters: {a, b} and {c}. With d(a, b) >= sigma every pair is at sigma,
  # so multiplicity is the ratio of inverse Simpson indices after / before
  ids <- c("a", "b", "c")
  ab <- one_row(c(2, 3, 5), ids)
  clust <- c(1, 1, 2)
  sig <- 0.5
  p <- c(2, 3, 5) / 10
  p_clust <- c(5, 5) / 10
  expected <- sum(p_clust^2) / sum(p^2)

  for (d in c(0.5, 0.7, 3)) {
    diss <- data.frame(ID1 = "a", ID2 = "b", Distance = d)
    expect_equal(unname(multiplicity.distance(ab, diss, clust, sig = sig)), expected, info = d)
  }

  # Below sigma the distance counts: Q_before = sig (1 - sum p^2) + 2 p_a p_b (d - sig)
  d <- 0.2
  q_before <- sig * (1 - sum(p^2)) + 2 * p[1] * p[2] * (d - sig)
  q_after <- sig * (1 - sum(p_clust^2))
  expect_equal(
    unname(multiplicity.distance(ab, data.frame(ID1 = "a", ID2 = "b", Distance = d), clust, sig = sig)),
    (sig - q_after) / (sig - q_before)
  )
})
