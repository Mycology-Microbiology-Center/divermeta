## Tests for distance-based functional diversity and redundancy

test_that("diversity.functional matches capped Rao formula (2 species)", {
  ab <- one_row(c(2, 1))
  sig <- 0.5
  d <- 0.7
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = d)

  # diversity.functional caps distances at sig internally
  diss_capped <- data.frame(ID1 = "1", ID2 = "2", Distance = sig)

  expect_equal(diversity.functional(ab, diss, sig), diversity.functional(ab, diss_capped, sig))

  # Closed form: D_sigma = 1 / (1 - Q/sig), Q = 2 p1 p2 min(d, sig)
  p <- c(2, 1) / 3
  Q <- 2 * p[1] * p[2] * sig
  expected <- 1 / (1 - Q / sig)
  expect_equal(unname(diversity.functional(ab, diss, sig)), expected, tolerance = 1e-12)
})


# n equally abundant subunits, all at the maximal distance: both functional
# diversities equal n
test_that("functional diversities equal n for n equally abundant, maximally distant subunits", {
  sig <- 0.5
  n <- 10
  ab <- one_row(rep(1, n))
  diss <- mat_to_table(max_diss(n, sig))

  expect_equal(round(unname(diversity.functional.traditional(ab, diss)), 3), n)
  expect_equal(unname(diversity.functional(ab, diss, sig)), n)
})


# Three identical copies of an assemblage, at the maximal distance from each
# other, have three times the diversity of one copy
test_that("diversity.functional triples for three maximally distant copies of an assemblage", {
  set.seed(42)
  sig <- 0.9
  n <- 10
  ab_unit <- runif(n, 100, 1000)
  diss_unit <- matrix(runif(n * n, min = 0.3, max = 0.6), nrow = n, ncol = n)
  diss_unit <- (diss_unit + t(diss_unit)) / 2
  diag(diss_unit) <- 0

  # Assemblage
  diss <- block_matrix(list(diss_unit, diss_unit, diss_unit), sig)
  ab <- c(ab_unit, ab_unit, ab_unit)

  expect_equal(
    unname(diversity.functional(one_row(ab), mat_to_table(diss), sig) /
      diversity.functional(one_row(ab_unit), mat_to_table(diss_unit), sig)),
    3
  )
})


test_that("diversity.functional.traditional is non-increasing in q", {
  set.seed(1)
  n <- 6
  ab <- one_row(runif(n, 1, 10))
  diss <- mat_to_table(rand_diss(n))

  q_values <- c(0, 0.5, 0.9, 1, 1.1, 2, 5)
  vals <- vapply(q_values, function(q) unname(diversity.functional.traditional(ab, diss, q)), numeric(1))
  expect_true(all(diff(vals) <= 1e-10))
  # ... and strictly decreasing for unequal abundances and distances
  expect_true(vals[1] > vals[length(vals)])
})


test_that("diversity.functional.traditional two-species across q", {
  ab <- one_row(c(3, 1))
  d <- 0.4
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = d)

  # With two species the value does not depend on q; check q = 1 continuity
  q_values <- c(0.5, 0.9, 1.1, 2)
  vals <- vapply(q_values, function(q) unname(diversity.functional.traditional(ab, diss, q)), numeric(1))

  # Ensure all finite and positive
  expect_true(all(is.finite(vals) & vals > 0))

  # q near 1 on both sides should be close
  expect_equal(vals[which(q_values == 0.9)], vals[which(q_values == 1.1)], tolerance = 1e-2)

  # q = 1 exactly: closed form exp(-sum_ij d_ij p_i p_j ln(p_i) / Q), and the legacy value
  p <- c(3, 1) / 4
  Q <- 2 * p[1] * p[2] * d
  expected <- exp(-d * p[1] * p[2] * (log(p[1]) + log(p[2])) / Q)
  expect_equal(unname(diversity.functional.traditional(ab, diss, 1)), expected, tolerance = 1e-12)
  expect_equal(
    unname(diversity.functional.traditional(ab, diss, 1)),
    legacy$diversity.functional.traditional(c(3, 1), matrix(c(0, d, d, 0), 2), 1)
  )
})


test_that("diversity.functional.traditional is 1 when Q = 0", {
  ids <- c("a", "b", "c")
  diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))
  zero <- transform(diss, Distance = 0)
  abund <- matrix(
    c(5, 0, 0,  # one present subunit
      2, 3, 1,  # all distances 0 (with `zero`)
      0, 0, 0), # empty
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), ids)
  )

  for (q in c(0, 0.5, 1, 2)) {
    res <- diversity.functional.traditional(abund, diss, q)
    expect_equal(res[["S1"]], 1)
    expect_true(is.na(res[["S3"]]))
    expect_equal(res[["S1"]], diversity.functional(abund, diss)[["S1"]])

    res <- diversity.functional.traditional(abund, zero, q)
    expect_equal(res[["S2"]], 1)
    expect_true(is.na(res[["S3"]]))
  }
})


test_that("diversity.functional.traditional is stable for q close to 1", {
  set.seed(42)
  ids <- paste0("f", 1:6)
  tab <- mat_to_table(rand_diss(6), ids)
  abund <- matrix(runif(18, 1, 10), nrow = 3, dimnames = list(paste0("S", 1:3), ids))
  at_one <- diversity.functional.traditional(abund, tab, 1)

  for (q in 1 + c(-1e-6, -1e-10, -1e-15, 1e-15, 1e-10, 1e-6)) {
    expect_equal(diversity.functional.traditional(abund, tab, q), at_one, tolerance = 1e-5)
  }
})


test_that("q must be a single non-negative finite number", {
  abund <- one_row(c(1, 2))
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = 0.2)
  for (q in list(NA_real_, Inf, NA, "1")) {
    expect_error(
      diversity.functional.traditional(abund, diss, q = q),
      "`q` must be a single non-negative numeric value"
    )
  }
})


test_that("diversity.functional and diversity.functional.traditional equal the legacy implementation", {
  set.seed(42)
  n <- 15
  ids <- paste0("f", seq_len(n))
  diss <- rand_diss(n, 0.1, 1.2)
  abund <- matrix(
    runif(6 * n, 1, 100) * (runif(6 * n) > 0.3),
    nrow = 6,
    dimnames = list(paste0("S", 1:6), ids)
  )
  tab <- mat_to_table(diss, ids)
  sparse <- Matrix::Matrix(abund, sparse = TRUE)

  for (sig in c(0.3, 0.7, 1)) {
    expected <- apply(abund, 1, function(ab) legacy$diversity.functional(ab, diss, sig))
    expect_equal(diversity.functional(abund, tab, sig), expected)
    expect_equal(diversity.functional(sparse, tab, sig), expected)
  }

  # Legacy is given only the present subunits (for q = 0, absent subunits would count)
  for (q in c(0, 0.5, 1, 2, 3)) {
    expected <- apply(abund, 1, function(ab) {
      legacy$diversity.functional.traditional(ab[ab > 0], diss[ab > 0, ab > 0], q)
    })
    expect_equal(diversity.functional.traditional(abund, tab, q), expected)
    expect_equal(diversity.functional.traditional(sparse, tab, q), expected)
  }

  # For q > 0, zero abundances do not change the legacy value
  for (q in c(0.5, 1, 2)) {
    expected <- apply(abund, 1, function(ab) legacy$diversity.functional.traditional(ab, diss, q))
    expect_equal(diversity.functional.traditional(abund, tab, q), expected)
  }

  # Validation
  expect_error(diversity.functional.traditional(abund, tab, q = -1), "non-negative")
  expect_error(diversity.functional.traditional(abund, tab, q = c(1, 2)), "single non-negative numeric")
  expect_error(diversity.functional(abund, tab, sig = 0), "positive")
})


test_that("redundancy numeric and bounds (2 species)", {
  ab <- c(2, 1)
  d <- 0.3
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = d)

  # Redundancy = (1 - sum p_i^2) - Q
  p <- ab / sum(ab)
  D <- 1 - sum(p^2)
  Q <- 2 * p[1] * p[2] * d
  expected <- D - Q

  expect_equal(unname(redundancy(one_row(ab), diss)), expected, tolerance = 1e-12)

  # Bounds: redundancy <= D and >= 0 when diss in [0,1]
  expect_true(redundancy(one_row(ab), diss) <= D + 1e-12)
  expect_true(redundancy(one_row(ab), diss) >= 0 - 1e-12)
})


test_that("redundancy input validation", {
  ab <- one_row(c(1, 2))
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = 0.2)

  expect_error(redundancy("a", diss), "must be a matrix")
  # A string is the path to a distance file
  expect_error(redundancy(ab, "not a matrix"), "not found")
  expect_error(redundancy(ab, matrix(1, 2, 3)), "Missing distances")
  expect_error(redundancy(one_row(c(-1, 2)), diss), "non-negative")

  # Empty samples give NA
  expect_true(is.na(redundancy(one_row(c(0, 0)), diss)))

  # A single subunit has no redundancy
  expect_equal(unname(redundancy(one_row(5, "a"), diss)), 0)
})

# Builds the full matrix (1 off-block) and its complete table from a list of blocks
build_blocks <- function(blocks) {
  sizes <- vapply(blocks, nrow, integer(1))
  n <- sum(sizes)
  ids <- as.character(seq_len(n))
  diss <- matrix(1, nrow = n, ncol = n)
  offset <- 0
  for (k in seq_along(blocks)) {
    idx <- offset + seq_len(sizes[k])
    diss[idx, idx] <- blocks[[k]]
    offset <- offset + sizes[k]
  }
  diag(diss) <- 0
  list(ids = ids, diss = diss, diss_frame = mat_to_table(diss, ids))
}

random_block <- function(size, min = 0, max = 1) {
  M <- matrix(runif(size^2, min, max), ncol = size, nrow = size)
  M <- (M + t(M)) / 2
  diag(M) <- 0
  M
}


test_that("redundancy numeric check (2 blocks)", {
  ids <- c("a", "b", "c", "d")
  ab <- c(2, 3, 5, 7)

  # Re = 2 sum_P p_i p_j (1 - d_ij), over the pairs P closer than 1
  p <- ab / sum(ab)
  expected <- 2 * (p[1] * p[2] * (1 - 0.2) + p[3] * p[4] * (1 - 0.15))

  diss <- matrix(1, nrow = 4, ncol = 4, dimnames = list(ids, ids))
  diag(diss) <- 0
  diss["a", "b"] <- diss["b", "a"] <- 0.2
  diss["c", "d"] <- diss["d", "c"] <- 0.15

  rb <- redundancy(one_row(ab, ids), mat_to_table(diss))

  expect_equal(unname(rb), expected, tolerance = 1e-12)
  expect_equal(unname(rb), legacy$redundancy(ab, diss), tolerance = 1e-12)
})


test_that("redundancy equals the legacy implementation for multiple blocks", {
  set.seed(42)

  for (w_ in 1:10) {
    sizes <- sample(2:10, sample(3:4, 1), replace = TRUE)
    built <- build_blocks(lapply(sizes, random_block))
    abund <- matrix(
      runif(3 * length(built$ids), min = 1, max = 100),
      nrow = 3,
      dimnames = list(NULL, built$ids)
    )

    rb <- redundancy(abund, built$diss_frame)
    classic <- apply(abund, 1, function(ab) legacy$redundancy(ab, built$diss))

    expect_equal(unname(rb), unname(classic), tolerance = 1e-12)
  }
})


test_that("redundancy caps distances greater than 1 and is never negative", {
  set.seed(7)
  built <- build_blocks(list(random_block(4, 0.5, 1.5), random_block(5, 0.5, 1.5)))
  abund <- matrix(
    runif(4 * length(built$ids), min = 1, max = 10) * (runif(4 * length(built$ids)) > 0.2),
    nrow = 4,
    dimnames = list(paste0("S", 1:4), built$ids)
  )
  capped <- transform(built$diss_frame, Distance = pmin(Distance, 1))

  for (normalize in c(FALSE, TRUE)) {
    res <- redundancy(abund, built$diss_frame, normalize = normalize)
    # Same as the distances capped at 1 beforehand
    expect_equal(res, redundancy(abund, capped, normalize = normalize), tolerance = 1e-12)
    # Never negative, and at most Simpson's diversity (1 when normalized)
    expect_true(all(res >= 0))
    upper <- if (normalize) 1 else 1 - rowSums((abund / rowSums(abund))^2)
    expect_true(all(res <= upper + 1e-12))
  }

  # Same as the legacy (matrix) implementation given the capped distances
  ab <- abund[1, ]
  expect_equal(
    unname(redundancy(one_row(ab, built$ids), built$diss_frame)),
    legacy$redundancy(ab, pmin(built$diss, 1)),
    tolerance = 1e-12
  )

  # Every distance above 1: no redundancy at all (the uncapped index would be negative)
  far <- transform(built$diss_frame, Distance = 3)
  expect_equal(unname(redundancy(abund, far)), rep(0, 4))
  expect_equal(unname(redundancy(abund, far, normalize = TRUE)), rep(0, 4))
})


test_that("average redundancy caps distances greater than 1", {
  clust <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B")
  diss <- data.frame(
    ID1 = c("a1", "a1", "a2", "b1"),
    ID2 = c("a2", "a3", "a3", "b2"),
    Distance = c(0.3, 2.5, 1.7, 4)
  )
  capped <- transform(diss, Distance = pmin(Distance, 1))
  abund <- rbind(S1 = c(5, 3, 1, 2, 2), S2 = c(1, 1, 1, 0, 4))
  colnames(abund) <- names(clust)

  for (normalize in c(FALSE, TRUE)) {
    res <- average.redundancy(abund, diss, clust, normalize = normalize)
    expect_equal(res, average.redundancy(abund, capped, clust, normalize = normalize))
    expect_true(all(res >= 0))
  }
})


test_that("redundancy handles singleton units", {
  set.seed(3)
  built <- build_blocks(list(random_block(3), matrix(0, 1, 1), random_block(4), matrix(0, 1, 1)))
  ab <- runif(length(built$ids), min = 1, max = 10)

  expect_equal(
    unname(redundancy(one_row(ab, built$ids), built$diss_frame)),
    legacy$redundancy(ab, built$diss),
    tolerance = 1e-12
  )

  # All singletons (every pair at distance 1): no redundancy
  built <- build_blocks(list(matrix(0, 1, 1), matrix(0, 1, 1), matrix(0, 1, 1)))
  expect_equal(unname(redundancy(one_row(c(1, 2, 3), built$ids), built$diss_frame)), 0)
})


test_that("normalized redundancy is 1 - Q / S (2 species)", {
  ab <- c(2, 1)
  d <- 0.3
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = d)

  p <- ab / sum(ab)
  D <- 1 - sum(p^2)
  Q <- 2 * p[1] * p[2] * d
  expected <- (D - Q) / D

  res <- unname(redundancy(one_row(ab), diss, normalize = TRUE))
  expect_equal(res, expected, tolerance = 1e-12)
  expect_equal(res, 1 - Q / D, tolerance = 1e-12)
  expect_equal(res, unname(redundancy(one_row(ab), diss)) / D, tolerance = 1e-12)

  # With two subunits, S = 2 p1 p2, so the normalized value is 1 - d
  expect_equal(res, 1 - d, tolerance = 1e-12)
})


test_that("normalized redundancy equals redundancy over Simpson's diversity", {
  set.seed(11)
  for (w_ in 1:5) {
    sizes <- sample(2:8, sample(2:4, 1), replace = TRUE)
    built <- build_blocks(lapply(sizes, random_block))
    abund <- matrix(
      runif(4 * length(built$ids), min = 1, max = 100) * (runif(4 * length(built$ids)) > 0.3),
      nrow = 4,
      dimnames = list(paste0("S", 1:4), built$ids)
    )
    abund[, 1] <- abund[, 1] + 1

    P <- abund / rowSums(abund)
    D <- 1 - rowSums(P^2)
    res <- redundancy(abund, built$diss_frame, normalize = TRUE)
    expect_equal(res, redundancy(abund, built$diss_frame) / D, tolerance = 1e-12)

    # Between 0 and 1 for distances in [0, 1]
    expect_true(all(res >= -1e-12 & res <= 1 + 1e-12))

    # Invariant to the scale of the abundances
    expect_equal(redundancy(abund * 7.3, built$diss_frame, normalize = TRUE), res, tolerance = 1e-12)
  }
})


test_that("normalized redundancy edge cases", {
  ids <- c("a", "b", "c")
  ab <- one_row(c(1, 2, 3), ids)
  zero <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = 0)
  one <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = 1)

  # All subunits identical: fully redundant. All maximally different: no redundancy
  expect_equal(unname(redundancy(ab, zero, normalize = TRUE)), 1)
  expect_equal(unname(redundancy(ab, one, normalize = TRUE)), 0)

  # A single present subunit has no redundancy, also when 49 * (1 / 49) != 1
  for (a in c(1, 3, 49, 1e7 + 1)) {
    single <- one_row(c(0, a, 0), ids)
    expect_identical(unname(redundancy(single, zero, normalize = TRUE)), 0, info = a)
    expect_identical(unname(redundancy(single, zero)), 0, info = a)
  }

  # Empty samples give NA
  expect_true(is.na(redundancy(one_row(c(0, 0, 0), ids), zero, normalize = TRUE)))

  # Several samples at once
  abund <- rbind(S1 = c(1, 2, 3), S2 = c(0, 5, 0), S3 = c(0, 0, 0))
  colnames(abund) <- ids
  expect_equal(unname(redundancy(abund, zero, normalize = TRUE)), c(1, 0, NA))

  # Validation of the flag
  expect_error(redundancy(ab, zero, normalize = NA), "single TRUE or FALSE")
  expect_error(redundancy(ab, zero, normalize = "yes"), "single TRUE or FALSE")
  expect_error(redundancy(ab, zero, normalize = c(TRUE, FALSE)), "single TRUE or FALSE")
})


test_that("pairs listed more than once must have the same distance", {
  abund <- one_row(c(1, 2, 3), c("a", "b", "c"))
  df <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.2, 0.5, 0.7))
  base <- redundancy(abund, df)

  # Same pair in both orders, or repeated, with the same distance
  df_both <- rbind(df, data.frame(ID1 = "b", ID2 = "a", Distance = 0.2))
  df_repeat <- rbind(df, df[1, ])
  expect_equal(redundancy(abund, df_both), base)
  expect_equal(redundancy(abund, df_repeat), base)

  # With a different distance
  df_conflict <- rbind(df, data.frame(ID1 = "b", ID2 = "a", Distance = 0.3))
  expect_error(redundancy(abund, df_conflict), "conflicting")

  # Self pairs and ids not in `abund` are dropped
  df_ok <- rbind(df, data.frame(ID1 = c("a", "x"), ID2 = c("a", "a"), Distance = c(0, 0.1)))
  expect_equal(redundancy(abund, df_ok), base)
})

test_that("MAD numeric check", {
  ids <- as.character(1:6)
  clust <- c(1, 2, 2, 3, 3, 3)

  diss <- matrix(
    c(
      0.0, 0.3, 0.6, 0.4, 0.3, 0.8,
      0.3, 0.0, 0.7, 0.4, 0.4, 0.7,
      0.6, 0.7, 0.0, 0.7, 0.5, 0.5,
      0.4, 0.4, 0.7, 0.0, 0.5, 0.6,
      0.3, 0.4, 0.5, 0.5, 0.0, 0.5,
      0.1, 0.5, 0.6, 0.3, 0.5, 0.6
    ),
    nrow = 6,
    byrow = TRUE
  )

  # Distances from the representatives (rows of the matrix): 3 -> 2, 4 -> 5, 4 -> 6, 2 -> 3
  diss_frame <- data.frame(
    ID1 = c("3", "4", "4", "2"),
    ID2 = c("2", "5", "6", "3"),
    Distance = c(diss[3, 2], diss[4, 5], diss[4, 6], diss[2, 3])
  )

  representatives <- c("1", "3", "4")
  names(representatives) <- c(1, 2, 3)

  val <- (1 + 0) + (1 + (0.7 + 0) / 2) + (1 + (0.0 + 0.5 + 0.6) / 3)
  val <- val / 6

  mad_val <- metagenomic.alpha.index(one_row(rep(1, 6), ids), diss_frame, clust, representatives)

  expect_equal(unname(mad_val), val, tolerance = 1e-12)
  expect_equal(unname(mad_val), legacy$metagenomic.alpha.index(clust, diss, c("1" = 1, "2" = 3, "3" = 4)))

  # Change representatives
  representatives <- c("1", "2", "4")
  names(representatives) <- c(1, 2, 3)
  expect_equal(
    metagenomic.alpha.index(one_row(rep(1, 6), ids), diss_frame, clust, representatives),
    metagenomic.alpha.index(one_row(rep(1, 6), ids), diss_frame, clust, NULL),
    tolerance = 1e-12
  )
})


test_that("MAD equals the legacy implementation on the present subunits of every sample", {
  set.seed(42)
  n <- 40
  n_clust <- 6
  ids <- paste0("g", seq_len(n))
  clust <- sample(paste0("C", seq_len(n_clust)), n, replace = TRUE)
  diss <- rand_diss(n, 0, 1)
  tab <- mat_to_table(diss, ids)
  abund <- matrix(
    rpois(5 * n, 3) * (runif(5 * n) > 0.4),
    nrow = 5,
    dimnames = list(paste0("S", 1:5), ids)
  )
  abund["S5", ] <- 0

  res <- metagenomic.alpha.index(abund, tab, clust)
  expected <- apply(abund, 1, function(ab) {
    if (sum(ab) == 0) {
      return(NA_real_)
    }
    present <- ab > 0
    legacy$metagenomic.alpha.index(clust[present], diss[present, present])
  })
  expect_equal(res, expected)
  expect_true(is.na(res[["S5"]]))

  # Sparse abundances
  expect_equal(metagenomic.alpha.index(Matrix::Matrix(abund, sparse = TRUE), tab, clust), res)

  # A representative absent from the sample is still the reference of its cluster
  ab <- c(a1 = 0, a2 = 1, a3 = 1, b1 = 1)
  d <- data.frame(ID1 = c("a1", "a1", "a2"), ID2 = c("a2", "a3", "a3"), Distance = c(0.2, 0.4, 0.9))
  reps <- c(A = "a1", B = "b1")
  cl <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B")
  expect_equal(
    unname(metagenomic.alpha.index(one_row(ab, names(ab)), d, cl, reps)),
    ((1 + (0.2 + 0.4) / 2) + (1 + 0)) / 3
  )
})


test_that("MAD input validation", {
  set.seed(1)
  n <- 100
  n_clust <- 10
  ids <- as.character(seq_len(n))
  clust <- sample(seq_len(n_clust),
    size = n,
    replace = TRUE
  )

  representatives <- sapply(seq_len(n_clust), function(cluster_name) {
    ids[sample(which(cluster_name == clust), 1)]
  })

  names(representatives) <- seq_len(n_clust)

  # Parameters
  diss <- matrix(runif(n * n, min = 0, max = 1), ncol = n, nrow = n)
  diss <- (diss + t(diss)) / 2
  diag(diss) <- 0
  tab <- mat_to_table(diss, ids)
  abund <- one_row(rep(1, n), ids)

  expect_silent(metagenomic.alpha.index(abund, tab, clust, representatives))
  expect_error(metagenomic.alpha.index(abund, "not a matrix", clust, representatives), "not found")
  expect_error(metagenomic.alpha.index(abund, tab, clust[1:5], representatives), "same number of subunits")
  expect_error(metagenomic.alpha.index(abund, tab, c(clust[1:(n - 1)], NA), representatives), "NA")

  new_representatives <- representatives
  new_representatives[1] <- "not a subunit"
  expect_error(metagenomic.alpha.index(abund, tab, clust, new_representatives), "subunits")

  # Representative from another cluster
  new_representatives <- representatives
  new_representatives[1] <- representatives[2]
  expect_error(metagenomic.alpha.index(abund, tab, clust, new_representatives), "belong")

  expect_error(metagenomic.alpha.index(abund, tab, clust, representatives[1:(n_clust - 1)]), "Missing representatives")
  expect_error(metagenomic.alpha.index(abund, tab, clust, unname(representatives)), "named")

  # Missing distance between a representative and a subunit of its cluster
  rep_1 <- representatives[["1"]]
  other <- setdiff(ids[clust == 1], rep_1)[1]
  drop <- (tab$ID1 == rep_1 & tab$ID2 == other) | (tab$ID1 == other & tab$ID2 == rep_1)
  expect_error(metagenomic.alpha.index(abund, tab[!drop, ], clust, representatives), "Missing distances")
})


test_that("MAD matches numeric cluster and subunit identifiers as numbers", {
  # Numbers of 1e5 or more are written "1e+05" by names() and colnames(), and
  # "100000" by the identifiers of the clustering
  subs <- c(100000, 200000, 300000, 400000, 500000)
  abund <- matrix(c(2, 1, 3, 1, 4, 0, 5, 1, 2, 2), nrow = 2, byrow = TRUE)
  colnames(abund) <- subs
  clust <- c(100000, 100000, 100000, 200000, 200000)
  tab <- data.frame(
    ID1 = c(100000, 100000, 400000),
    ID2 = c(200000, 300000, 500000),
    d = c(0.2, 0.6, 0.4)
  )
  reps <- c(100000, 400000)
  names(reps) <- c(100000, 200000)
  expect_identical(names(reps), c("1e+05", "2e+05"))

  # Same data with text identifiers
  abund_chr <- abund
  colnames(abund_chr) <- c("a", "b", "c", "d", "e")
  tab_chr <- data.frame(ID1 = c("a", "a", "d"), ID2 = c("b", "c", "e"), d = tab$d)
  expected <- metagenomic.alpha.index(abund_chr, tab_chr, c("A", "A", "A", "B", "B"), c(A = "a", B = "d"))

  expect_equal(unname(metagenomic.alpha.index(abund, tab, clust, reps)), unname(expected))
  # Representatives given as text, in either format
  expect_equal(unname(metagenomic.alpha.index(abund, tab, clust, stats::setNames(as.character(reps), names(reps)))), unname(expected))
  reps_txt <- stats::setNames(c("100000", "400000"), c("100000", "200000"))
  expect_equal(unname(metagenomic.alpha.index(abund, tab, clust, reps_txt)), unname(expected))
  # Extra names are ignored
  expect_equal(unname(metagenomic.alpha.index(abund, tab, clust, c(reps, X = "a"))), unname(expected))

  expect_error(metagenomic.alpha.index(abund, tab, clust, reps[1]), "Missing representatives for clusters: 200000")
  bad <- reps
  bad[2] <- 600000
  expect_error(metagenomic.alpha.index(abund, tab, clust, bad), "subunits of `abund`: 600000")
  bad[2] <- 200000
  expect_error(metagenomic.alpha.index(abund, tab, clust, bad), "belong to their cluster: 200000")
})
