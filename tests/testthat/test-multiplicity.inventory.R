
test_that("High Multlicity Works", {

  ab <- rep(10,30)
  clust <- c(rep(1,10), rep(2,10), rep(3,10))
  # Multiplicity should be 10
  expect_equal(multiplicity.inventory(ab, clust), 10)
})


test_that("Low Multlicity Works", {

  ab <- rep(10,3)
  clust <- c(1,2,3)
  # Multiplicity should be 1
  expect_equal(multiplicity.inventory(ab, clust), 1)
})
