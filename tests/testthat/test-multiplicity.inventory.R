
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
    by_implementation  <- multiplicity.inventory(ab, clust, q)

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
    by_implementation  <- multiplicity.inventory(ab, clust, q)

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

  expect_equal(multiplicity.inventory(ab, clust, q), expected, tolerance = 1e-12)
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

  expect_equal(multiplicity.inventory(ab, clust, q), expected, tolerance = 1e-12)
})


test_that("Single cluster multiplicity equals overall diversity", {
  ab <- c(3, 4, 5)
  clust <- c(1, 1, 1)

  p <- ab / sum(ab)

  # q = 1
  expected_q1 <- exp(-sum(p * log(p)))
  expect_equal(multiplicity.inventory(ab, clust, q = 1), expected_q1, tolerance = 1e-12)

  # q = 2
  expected_q2 <- (sum(p^2))^(1 / (1 - 2))
  expect_equal(multiplicity.inventory(ab, clust, q = 2), expected_q2, tolerance = 1e-12)

  # q = 0.5
  q <- 0.5
  expected_q05 <- (sum(p^q))^(1 / (1 - q))
  expect_equal(multiplicity.inventory(ab, clust, q = q), expected_q05, tolerance = 1e-12)
})

test_that("High multiplicity: many diverse elements per cluster", {
  # Three clusters, each with 10 equally abundant elements
  # Multiplicity should equal the number of elements per cluster
  ab <- rep(10, 30)
  clust <- c(rep(1, 10), rep(2, 10), rep(3, 10))
  expect_equal(multiplicity.inventory(ab, clust), 10)
})

test_that("Low multiplicity: one element per cluster", {
  # Three clusters, each with one element
  # Multiplicity should be 1 (no diversity lost)
  ab <- rep(10, 3)
  clust <- c(1, 2, 3)
  expect_equal(multiplicity.inventory(ab, clust), 1)
})


test_that("Zero abundances are automatically removed", {
  # Elements with zero abundance should be ignored
  ab_1 <- c(10, 3, 3, 5, 7, 1, 2, 3)
  clust_1 <- c(1, 2, 3, 2, 3, 2, 1, 2)

  ab_2 <- c(10, 3, 3, 5, 7, 1, 2, 3, 0, 0, 0, 0)
  clust_2 <- c(1, 2, 3, 2, 3, 2, 1, 2, 1, 2, 3, 1)

  # Multiplicity should be equal (zero-abundance elements ignored)
  expect_equal(multiplicity.inventory(ab_1, clust_1), multiplicity.inventory(ab_2, clust_2))
})

test_that("Input validation works correctly", {
  # Non-numeric abundance
  expect_error(multiplicity.inventory(c("a", "b"), c(1, 2)), "must be a numeric vector")
  
  # Negative q
  expect_error(multiplicity.inventory(c(1, 2), c(1, 2), q = -1), "must be a single non-negative")
  
  # Mismatched lengths
  expect_error(multiplicity.inventory(c(1, 2, 3), c(1, 2)), "must have the same length")
  
  # All zero abundances
  expect_equal(multiplicity.inventory(c(0, 0, 0), c(1, 2, 3)), 0)
})


