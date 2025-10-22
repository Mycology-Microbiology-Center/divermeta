
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

