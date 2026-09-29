
## Tests for inventory multiplicity

test_that("Formula matches q != 1", {


  set.seed(42)
  for(i in 1:15)
  {

    q  <-  runif(1, 0.001,10)
    if(q == 1)
      q  <- q + 1e-10

    # Three cluster
    ab_1 <- runif(4,1,100)
    ab_2 <- runif(6,1,100)
    ab_3 <- runif(15,1,100)

    # By formula
    total <- sum(ab_1) + sum(ab_2) + sum(ab_3)
    P_1 <- ab_1 / sum(ab_1)
    P_2 <- ab_2 / sum(ab_2)
    P_3 <- ab_3 / sum(ab_3)

    numenator <-  (sum(ab_1)**q)*sum(P_1**q) + (sum(ab_2)**q)*sum(P_2**q)  + (sum(ab_3)**q)*sum(P_3**q)
    denominator <- sum(ab_1)**q + sum(ab_2)**q  + sum(ab_3)**q

    by_formula <- (numenator / denominator)**(1/(1-q))

    ab <- c(ab_1, ab_2, ab_3)
    clust <- c(rep(1, length(ab_1)), rep(2, length(ab_2)), rep(3, length(ab_3)))
    by_implementation  <- unname(multiplicity.inventory(one_row(ab), clust, q))

    expect_equal(by_formula, by_implementation)

  }



})


test_that("Formula matches q = 1", {


  q <- 1
  set.seed(42)
  for(i in 1:15)
  {

    # Three cluster
    ab_1 <- runif(4,1,100)
    ab_2 <- runif(6,1,100)
    ab_3 <- runif(15,1,100)

    # By formula
    N <- sum(ab_1) + sum(ab_2) + sum(ab_3)
    N_1 <- sum(ab_1)
    N_2 <- sum(ab_2)
    N_3 <- sum(ab_3)

    ab <- c(ab_1, ab_2, ab_3)

    first <- (1/N)*(N_1*log(N_1) + N_2*log(N_2) + N_3*log(N_3))
    second <- (1/N)*sum(ab*log(ab))


    by_formula <- exp((first-second))

    clust <- c(rep(1, length(ab_1)), rep(2, length(ab_2)), rep(3, length(ab_3)))
    by_implementation  <- unname(multiplicity.inventory(one_row(ab), clust, q))

    expect_equal(by_formula, by_implementation)

  }


})



test_that("Deterministic q = 1 numeric", {

  # Abundances and clusters
  ab <- c(2, 2, 1, 1)
  clust <- c(1, 1, 2, 2)
  q <- 1

  # Expected by definition: exp(H(p)) / exp(H(p_clust))
  p <- ab / sum(ab)
  ab_cl <- tapply(ab, clust, sum)
  p_cl <- ab_cl / sum(ab_cl)

  H <- -sum(p * log(p))
  H_cl <- -sum(p_cl * log(p_cl))
  expected <- exp(H) / exp(H_cl)

  expect_equal(unname(multiplicity.inventory(one_row(ab), clust, q)), expected, tolerance = 1e-12)
})


test_that("Deterministic q = 2 numeric", {

  # Abundances and clusters (same as above) with q != 1
  ab <- c(2, 2, 1, 1)
  clust <- c(1, 1, 2, 2)
  q <- 2

  p <- ab / sum(ab)
  ab_cl <- tapply(ab, clust, sum)
  p_cl <- ab_cl / sum(ab_cl)

  # Hill numbers for q=2
  D <- (sum(p^q))^(1/(1 - q))
  D_cl <- (sum(p_cl^q))^(1/(1 - q))
  expected <- D / D_cl

  expect_equal(unname(multiplicity.inventory(one_row(ab), clust, q)), expected, tolerance = 1e-12)
})


test_that("Single cluster multiplicity equals overall diversity", {
  ab <- one_row(c(3, 4, 5))
  clust <- c(1, 1, 1)

  p <- c(3, 4, 5) / 12

  # q = 1
  expected_q1 <- exp(-sum(p * log(p)))
  expect_equal(unname(multiplicity.inventory(ab, clust, q = 1)), expected_q1, tolerance = 1e-12)

  # q = 2
  expected_q2 <- (sum(p^2))^(1 / (1 - 2))
  expect_equal(unname(multiplicity.inventory(ab, clust, q = 2)), expected_q2, tolerance = 1e-12)

  # q = 0.5
  q <- 0.5
  expected_q05 <- (sum(p^q))^(1 / (1 - q))
  expect_equal(unname(multiplicity.inventory(ab, clust, q = q)), expected_q05, tolerance = 1e-12)
})

test_that("High multiplicity: many diverse elements per cluster", {
  # Three clusters, each with 10 equally abundant elements
  # Multiplicity should equal the number of elements per cluster
  ab <- rep(10, 30)
  clust <- c(rep(1, 10), rep(2, 10), rep(3, 10))
  expect_equal(unname(multiplicity.inventory(one_row(ab), clust)), 10)
})

test_that("Low multiplicity: one element per cluster", {
  # Three clusters, each with one element
  # Multiplicity should be 1 (no diversity lost)
  ab <- rep(10, 3)
  clust <- c(1, 2, 3)
  expect_equal(unname(multiplicity.inventory(one_row(ab), clust)), 1)
})


test_that("Zero abundances are automatically removed", {
  # Elements with zero abundance should be ignored
  ab_1 <- c(10, 3, 3, 5, 7, 1, 2, 3)
  clust_1 <- c(1, 2, 3, 2, 3, 2, 1, 2)

  ab_2 <- c(10, 3, 3, 5, 7, 1, 2, 3, 0, 0, 0, 0)
  clust_2 <- c(1, 2, 3, 2, 3, 2, 1, 2, 1, 2, 3, 1)

  # Multiplicity should be equal (zero-abundance elements ignored), also for q = 0
  for (q in c(0, 1, 2)) {
    expect_equal(
      multiplicity.inventory(one_row(ab_1), clust_1, q),
      multiplicity.inventory(one_row(ab_2), clust_2, q)
    )
  }

  # A cluster absent from the sample does not count either
  ab_3 <- c(ab_1, 0, 0)
  clust_3 <- c(clust_1, 4, 4)
  for (q in c(0, 1, 2)) {
    expect_equal(
      multiplicity.inventory(one_row(ab_1), clust_1, q),
      multiplicity.inventory(one_row(ab_3), clust_3, q)
    )
  }
})

test_that("Input validation works correctly", {
  # Non-numeric abundance
  expect_error(multiplicity.inventory(one_row(c("a", "b")), c(1, 2)), "must be numeric")

  # Negative q
  expect_error(multiplicity.inventory(one_row(c(1, 2)), c(1, 2), q = -1), "must be a single non-negative")

  # Mismatched lengths
  expect_error(multiplicity.inventory(one_row(c(1, 2, 3)), c(1, 2)), "same number of subunits")

  # Named clust missing a subunit
  expect_error(multiplicity.inventory(one_row(c(1, 2, 3)), c("1" = 1, "2" = 2)), "missing some subunits")

  # All zero abundances
  expect_true(is.na(multiplicity.inventory(one_row(c(0, 0, 0)), c(1, 2, 3))))
})


test_that("Many samples equal the legacy implementation", {
  set.seed(42)
  n <- 40
  ids <- paste0("f", seq_len(n))
  clust <- stats::setNames(sample(c("x", "y", "z", "w"), n, replace = TRUE), ids)
  abund <- matrix(
    runif(10 * n, 1, 100) * (runif(10 * n) > 0.5),
    nrow = 10,
    dimnames = list(paste0("S", 1:10), ids)
  )
  abund["S10", ] <- 0
  sparse <- Matrix::Matrix(abund, sparse = TRUE)

  for (q in c(0, 0.5, 1, 2, 4)) {
    expected <- apply(abund, 1, function(ab) {
      if (sum(ab) == 0) {
        return(NA_real_)
      }
      legacy$multiplicity.inventory(ab, unname(clust), q)
    })
    res <- multiplicity.inventory(abund, clust, q)
    expect_equal(res, expected)
    expect_equal(multiplicity.inventory(sparse, clust, q), expected)
    expect_equal(multiplicity.inventory(as.data.frame(abund), clust, q), expected)

    # Named clust in a different order, and unnamed clust in column order
    expect_equal(multiplicity.inventory(abund, rev(clust), q), expected)
    expect_equal(multiplicity.inventory(abund, unname(clust), q), expected)
  }
})


test_that("Inventory multiplicity is stable for q close to 1", {
  set.seed(42)
  clust <- rep(c("A", "B", "C"), each = 4)
  abund <- matrix(runif(36, 1, 10), nrow = 3, dimnames = list(paste0("S", 1:3), paste0("f", 1:12)))
  at_one <- multiplicity.inventory(abund, clust, 1)

  for (q in 1 + c(-1e-6, -1e-10, -1e-13, -1e-15, 1e-15, 1e-10, 1e-6)) {
    expect_equal(multiplicity.inventory(abund, clust, q), at_one, tolerance = 1e-5)
  }
})


test_that("Inventory multiplicity requires a single non-negative finite q", {
  abund <- one_row(c(1, 2))
  for (q in list(NA_real_, Inf, -1, c(1, 2))) {
    expect_error(multiplicity.inventory(abund, c(1, 1), q), "`q` must be a single non-negative numeric value")
  }
})


# Clustering given as a data frame
# ----------------------------

test_that("clust can be a data frame for every function", {
  set.seed(42)
  n <- 12
  ids <- paste0("f", seq_len(n))
  clust <- stats::setNames(sample(rep(c("A", "B", "C"), length.out = n)), ids)
  abund <- matrix(
    runif(4 * n, 1, 50) * (runif(4 * n) > 0.3),
    nrow = 4,
    dimnames = list(paste0("S", 1:4), ids)
  )
  diss <- mat_to_table(rand_diss(n), ids)

  # Two columns (subunit, unit) in any order, one column with the subunits as
  # row names, one column following the columns of abund, and factors
  perm <- sample(n)
  forms <- list(
    two = data.frame(subunit = names(clust)[perm], unit = unname(clust)[perm]),
    two_factor = data.frame(subunit = factor(names(clust)), unit = factor(unname(clust))),
    rownames = data.frame(unit = unname(clust)[perm], row.names = names(clust)[perm]),
    unnamed = data.frame(unit = unname(clust)),
    matrix = matrix(clust, dimnames = list(names(clust), NULL))
  )

  for (form in names(forms)) {
    cl <- forms[[form]]
    expect_equal(multiplicity.inventory(abund, cl), multiplicity.inventory(abund, clust), info = form)
    expect_equal(
      multiplicity.distance(abund, diss, cl, sig = 0.6),
      multiplicity.distance(abund, diss, clust, sig = 0.6),
      info = form
    )
    expect_equal(
      multiplicity.distance(abund, diss, cl, method = "average"),
      multiplicity.distance(abund, diss, clust, method = "average"),
      info = form
    )
    expect_equal(
      metagenomic.alpha.index(abund, diss, cl),
      metagenomic.alpha.index(abund, diss, clust),
      info = form
    )
    expect_equal(
      relative.multiplicity(abund, diss, cl, sigma = c(A = 0.3, B = 0.6, C = 0.9)),
      relative.multiplicity(abund, diss, clust, sigma = c(A = 0.3, B = 0.6, C = 0.9)),
      info = form
    )
    expect_equal(
      average.redundancy(abund, diss, cl),
      average.redundancy(abund, diss, clust),
      info = form
    )
    expect_equal(
      divermeta(abund, diss, c("M_inventory", "M_distance"), cl),
      divermeta(abund, diss, c("M_inventory", "M_distance"), clust),
      info = form
    )
  }

  # Functions that need the subunit identifiers. The pairs of units follow the
  # order of the units in `clust`, so they are compared in a fixed order
  sorted_pairs <- function(ud) {
    key <- paste(pmin(ud$ID1, ud$ID2), pmax(ud$ID1, ud$ID2))
    stats::setNames(ud$Distance, key)[sort(key)]
  }
  for (form in c("two", "two_factor", "rownames", "matrix")) {
    cl <- forms[[form]]
    expect_equal(sorted_pairs(unit_distances(diss, cl)), sorted_pairs(unit_distances(diss, clust)), info = form)
    path <- write_diss(diss, "csv")
    expect_true(check_distance_file(path, ids, cl, scope = "within")$ok, info = form)
  }
  expect_error(unit_distances(diss, forms$unnamed), "named")

  # Unnamed values of relative multiplicity follow the order of the units in the data frame
  expect_equal(
    relative.multiplicity(abund, diss, forms$two, sigma = c(0.3, 0.6, 0.9)),
    relative.multiplicity(abund, diss, clust[forms$two$subunit], sigma = c(0.3, 0.6, 0.9))
  )

  # Other shapes are rejected with a clear message
  expect_error(
    multiplicity.inventory(abund, data.frame(a = ids, b = unname(clust), c = 1)),
    "one column .* or two columns"
  )
  expect_error(multiplicity.inventory(abund, cbind(clust, clust)), "single column")
  expect_error(multiplicity.inventory(abund, data.frame(unit = c(unname(clust)[-1], NA))), "NA")
  expect_error(
    multiplicity.inventory(abund, data.frame(subunit = ids[-1], unit = unname(clust)[-1])),
    "missing some subunits"
  )
})
