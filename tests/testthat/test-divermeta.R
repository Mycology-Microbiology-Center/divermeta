
## Tests for the `divermeta` function


## Function to create a test data
make_fixture <- function() {
  species <- c("f1", "f2", "f3", "f4")
  samples <- c("S1", "S2", "S_empty")
  abund <- matrix(
    c(
      10, 0,  0,   # f1
      5,  8,  0,   # f2
      0,  12, 0,   # f3
      6,  3,  0    # f4
    ),
    nrow = length(species),
    ncol = length(samples),
    byrow = TRUE,
    dimnames = list(species, samples))

  diss <- matrix(
    c(
      0.0, 0.2, 0.7, 1.0,
      0.2, 0.0, 0.6, 0.9,
      0.7, 0.6, 0.0, 0.4,
      1.0, 0.9, 0.4, 0.0
    ),
    nrow = length(species),
    ncol = length(species),
    byrow = TRUE,
    dimnames = list(species, species))

  clusters <- c(f1 = "A", f2 = "A", f3 = "B", f4 = "B")

  list(abund = abund, diss = diss, clusters = clusters)
}


test_that("single index multiplicity_inventory works and matches direct", {
  fx <- make_fixture()

  res <- divermeta(fx$abund,
    indices = c("multiplicity_inventory"),
    clusters = fx$clusters)

  expect_true(is.data.frame(res))
  expect_identical(colnames(res), c("Sample", "multiplicity_inventory"))
  expect_identical(res$Sample, colnames(fx$abund))

  expected <- apply(fx$abund, 2, function(col) {
    if (sum(col) <= 0) return(NA_real_)
    multiplicity.inventory(col, fx$clusters, q = 1)
  })
  expect_equal(res$multiplicity_inventory, as.numeric(expected))
})
