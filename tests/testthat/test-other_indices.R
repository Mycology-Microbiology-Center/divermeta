
## Tests for distance-based functional diversity and redundancy

test_that("diversity.functional (q=1) matches capped Rao formula (2 species)", {

  ab <- c(2, 1)
  sig <- 0.5
  d <- 0.7
  diss <- matrix(c(0, d, d, 0), 2, 2)

  # diversity.functional caps distances at sig internally
  diss_capped <- matrix(c(0, sig, sig, 0), 2, 2)

  expect_equal(diversity.functional(ab, diss, sig), diversity.functional(ab, diss_capped, sig))

  # Closed form: D_sigma = 1 / (1 - Q/sig), Q = 2 p1 p2 min(d, sig)
  p <- ab / sum(ab)
  Q <- 2 * p[1] * p[2] * sig
  expected <- 1 / (1 - Q / sig)
  expect_equal(diversity.functional(ab, diss, sig), expected, tolerance = 1e-12)
})


test_that("diversity.functional.traditional two-species across q", {

  ab <- c(3, 1)
  d <- 0.4
  diss <- matrix(c(0, d, d, 0), 2, 2)

  # For two species symmetric distances, check monotonicity with q and q=1 continuity
  q_values <- c(0.5, 0.9, 1.1, 2)
  vals <- vapply(q_values, function(q) diversity.functional.traditional(ab, diss, q), numeric(1))

  # Ensure all finite and positive
  expect_true(all(is.finite(vals) & vals > 0))

  # q near 1 on both sides should be close
  expect_equal(vals[which(q_values == 0.9)], vals[which(q_values == 1.1)], tolerance = 1e-2)
})


test_that("redundancy numeric and bounds (2 species)", {

  ab <- c(2, 1)
  d <- 0.3
  diss <- matrix(c(0, d, d, 0), 2, 2)

  # Redundancy = (1 - sum p_i^2) - Q
  p <- ab / sum(ab)
  D <- 1 - sum(p^2)
  Q <- 2 * p[1] * p[2] * d
  expected <- D - Q

  expect_equal(redundancy(ab, diss), expected, tolerance = 1e-12)

  # Bounds: redundancy <= D and >= 0 when diss in [0,1]
  expect_true(redundancy(ab, diss) <= D + 1e-12)
  expect_true(redundancy(ab, diss) >= 0 - 1e-12)
})


test_that("redundancy input validation", {

  ab <- c(1, 2)
  diss <- matrix(c(0, 0.2, 0.2, 0), 2, 2)

  expect_error(redundancy("a", diss))
  expect_error(redundancy(ab, "not a matrix"))
  expect_error(redundancy(ab, matrix(1, 2, 3)))
  expect_error(redundancy(c(-1, 2), diss))
  expect_error(redundancy(c(0, 0), diss))
})


