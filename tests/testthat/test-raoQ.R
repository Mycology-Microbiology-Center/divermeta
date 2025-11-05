## Tests for Rao quadratic entropy

test_that("raoQuadratic deterministic numeric (2 species)", {
  ab <- c(2, 1)
  d <- 0.4
  diss <- matrix(c(0, d, d, 0), nrow = 2, ncol = 2)

  # Q = 2 * p1 * p2 * d
  p <- ab / sum(ab)
  expected <- 2 * p[1] * p[2] * d

  expect_equal(raoQuadratic(ab, diss), expected, tolerance = 1e-12)
})


test_that("raoQuadratic permutation invariance", {
  ab <- c(2, 3, 5)
  diss <- matrix(c(
    0,   0.1, 0.7,
    0.1, 0,   0.3,
    0.7, 0.3, 0
  ), nrow = 3, byrow = TRUE)

  baseline <- raoQuadratic(ab, diss)

  # Permute order consistently
  idx <- c(3, 1, 2)
  ab2 <- ab[idx]
  diss2 <- diss[idx, idx]

  expect_equal(raoQuadratic(ab2, diss2), baseline, tolerance = 1e-12)
})


test_that("raoQuadratic input validation", {
  ab <- c(1, NA)
  diss <- matrix(c(0, 0.2, 0.2, 0), 2, 2)
  expect_error(raoQuadratic(ab, diss))

  ab <- c(-1, 2)
  expect_error(raoQuadratic(ab, diss))

  ab <- c(0, 0)
  expect_error(raoQuadratic(ab, diss))

  # Non-square diss or non-conformable dims should error via matrix mult
  ab <- c(1, 2)
  diss_bad <- matrix(c(0, 0.2, 0.2, 0, 0.1, 0.1), nrow = 2)
  expect_error(raoQuadratic(ab, diss_bad))
})


test_that("raoQuadratic dist and matrix implementation equivalence", {
  ab <- runif(30, min = 10, max = 100)
  vec <- rnorm(30)
  diss <- dist(vec)

  v1 <- raoQuadratic(ab = ab, diss = diss)
  v2 <- raoQuadratic(ab = ab, diss = as.matrix(diss))


  expect_equal(v1, v2)
})
