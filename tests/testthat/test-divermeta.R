
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


test_that("multiple indices and alias mapping; matches direct functions", {
  fx <- make_fixture()
  sig <- 0.8
  q <- 1

  res <- divermeta(fx$abund,
    diss = fx$diss,
    indices = c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq"),
    clusters = fx$clusters,
    q = q,
    sig = sig)

  # Column names are normalized
  expect_identical(colnames(res), c("Sample", "multiplicity_inventory", "raoQ", "FD_sigma", "redundancy", "multiplicity_distance", "FD_q"))

  # Expected values per column
  guard <- function(f) function(col) { if (sum(col) <= 0) return(NA_real_); f(col) }

  expected_mi <- apply(fx$abund, 2, guard(function(col) multiplicity.inventory(col, fx$clusters, q = q)))
  expected_rao <- apply(fx$abund, 2, guard(function(col) raoQuadratic(col, fx$diss)))
  expected_fd_sigma <- apply(fx$abund, 2, guard(function(col) diversity.functional(col, fx$diss, sig = sig)))

  ids <- rownames(fx$abund)
  ut_idx <- which(upper.tri(fx$diss), arr.ind = TRUE)
  diss_frame <- data.frame(
    ID1 = ids[ut_idx[, 1]],
    ID2 = ids[ut_idx[, 2]],
    Distance = fx$diss[upper.tri(fx$diss)],
    stringsAsFactors = FALSE
  )
  expected_md <- apply(fx$abund, 2, guard(function(col) multiplicity.distance.by_blocks(
    ids = ids,
    ab = col,
    diss_frame = diss_frame,
    clust = fx$clusters,
    sigma = sig
  )))

  expected_fd_q <- apply(fx$abund, 2, guard(function(col) diversity.functional.traditional(col, fx$diss, q = q)))

  expect_equal(res$multiplicity_inventory, as.numeric(expected_mi))
  expect_equal(res$raoQ, as.numeric(expected_rao))
  expect_equal(res$FD_sigma, as.numeric(expected_fd_sigma))
  expect_equal(res$multiplicity_distance, as.numeric(expected_md))
  expect_equal(res$FD_q, as.numeric(expected_fd_q))
})


