## Tests for the `divermeta` function


## Function to create a test data
make_fixture <- function() {
  species <- c("f1", "f2", "f3", "f4")
  samples <- c("S1", "S2", "S_empty")
  abund <- matrix(
    c(
      10, 0,  0, # f1
      5,  8,  0, # f2
      0,  12, 0, # f3
      6,  3,  0 # f4
    ),
    nrow = length(species),
    ncol = length(samples),
    byrow = TRUE,
    dimnames = list(species, samples)
  )

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
    dimnames = list(species, species)
  )

  clusters <- c(f1 = "A", f2 = "A", f3 = "B", f4 = "B")

  list(abund = abund, diss = diss, clusters = clusters)
}


test_that("single index multiplicity_inventory works and matches direct", {
  fx <- make_fixture()

  res <- divermeta(fx$abund,
    indices = c("multiplicity_inventory"),
    clusters = fx$clusters
  )

  expect_true(is.data.frame(res))
  expect_identical(colnames(res), c("Sample", "multiplicity_inventory"))
  expect_identical(res$Sample, colnames(fx$abund))

  expected <- apply(fx$abund, 2, function(col) {
    if (sum(col) <= 0) {
      return(NA_real_)
    }
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
    sig = sig
  )

  # Column names are normalized
  expect_identical(colnames(res), c("Sample", "multiplicity_inventory", "raoQ", "FD_sigma", "redundancy", "multiplicity_distance", "FD_q"))

  # Expected values per column
  guard <- function(f) {
    function(col) {
      if (sum(col) <= 0) {
        return(NA_real_)
      }
      f(col)
    }
  }

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
  expected_md <- apply(fx$abund, 2, guard(function(col) {
    multiplicity.distance.by_blocks(
      ids = ids,
      ab = col,
      diss_frame = diss_frame,
      clust = fx$clusters,
      sigma = sig
    )
  }))

  expected_fd_q <- apply(fx$abund, 2, guard(function(col) diversity.functional.traditional(col, fx$diss, q = q)))

  expect_equal(res$multiplicity_inventory, as.numeric(expected_mi))
  expect_equal(res$raoQ, as.numeric(expected_rao))
  expect_equal(res$FD_sigma, as.numeric(expected_fd_sigma))
  expect_equal(res$multiplicity_distance, as.numeric(expected_md))
  expect_equal(res$FD_q, as.numeric(expected_fd_q))

  res_norm <- divermeta(fx$abund,
    diss = fx$diss,
    indices = c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq"),
    clusters = fx$clusters,
    q = q,
    sig = sig,
    normalize = TRUE
  )

  # Checks normalization
  for (ind in c("multiplicity_inventory", "raoQ", "FD_sigma", "redundancy", "multiplicity_distance", "FD_q")) {
    expect_equal(max(res_norm[[ind]], na.rm = TRUE), 1.0)
  }
})



test_that("multiple indices and alias mapping; checks for dist object equivalence", {
  fx <- make_fixture()
  sig <- 0.8
  q <- 1

  res_m <- divermeta(fx$abund,
    diss = fx$diss,
    indices = c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq"),
    clusters = fx$clusters,
    q = q,
    sig = sig
  )

  res_d <- divermeta(fx$abund,
    diss = as.dist(fx$diss),
    indices = c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq"),
    clusters = fx$clusters,
    q = q,
    sig = sig
  )

  expect_equal(res_m$multiplicity_inventory, res_d$multiplicity_inventory)
  expect_equal(res_m$raoQ, res_d$raoQ)
  expect_equal(res_m$FD_sigma, res_d$FD_sigma)
  expect_equal(res_m$multiplicity_distance, res_d$multiplicity_distance)
  expect_equal(res_m$FD_q, res_d$FD_q)
})



test_that("errors when diss required but missing; and when clusters required but missing", {
  fx <- make_fixture()

  # diss-required indices
  expect_error(divermeta(fx$abund, indices = c("raoQ")))

  # clusters-required indices
  expect_error(divermeta(fx$abund, indices = c("multiplicity_inventory")))
})


test_that("unsupported index label errors clearly", {
  fx <- make_fixture()
  expect_error(divermeta(fx$abund, indices = c("not_an_index")))
})


test_that("non-numeric abund errors; and diss misalignment errors when no names", {
  fx <- make_fixture()

  bad_abund <- data.frame(
    f1 = c("1", "2", "3"),
    f2 = c("4", "5", "6"),
    f3 = c("7", "8", "9")
  )
  # Transpose to get rows as features; still character matrix
  bad_abund <- as.matrix(t(bad_abund))
  expect_error(divermeta(bad_abund, indices = c("multiplicity_inventory"), clusters = fx$clusters))

  # No names on diss; mismatched dims with abund should error
  abund2 <- fx$abund
  diss_bad <- matrix(0, nrow = 5, ncol = 5)
  expect_error(divermeta(abund2, diss = diss_bad, indices = c("raoQ")))
})


test_that("alignment by names with warnings and dropping/excess features", {
  fx <- make_fixture()

  # Add an extra feature to diss; reorder rows/cols
  diss2 <- rbind(fx$diss, extra = c(0.3, 0.3, 0.3, 0.3))
  diss2 <- cbind(diss2, extra = c(0.3, 0.3, 0.3, 0.3, 0.0))
  diss2 <- diss2[c("f3", "f1", "extra", "f4", "f2"), c("f3", "f1", "extra", "f4", "f2")]

  # Expect warnings about missing features in abund, and in diss
  expect_warning({
    res <- divermeta(fx$abund, diss = diss2, indices = c("raoQ"))
    expect_true(all(is.finite(res$raoQ[1:2])))
    expect_true(is.na(res$raoQ[3])) # empty column should be NA
  })
})


test_that("named clusters align regardless of order", {
  fx <- make_fixture()
  cl1 <- fx$clusters
  cl2 <- fx$clusters[c("f4", "f3", "f2", "f1")] # different order, still named

  r1 <- divermeta(fx$abund, indices = c("multiplicity_inventory"), clusters = cl1)
  r2 <- divermeta(fx$abund, indices = c("multiplicity_inventory"), clusters = cl2)

  expect_equal(r1$multiplicity_inventory, r2$multiplicity_inventory)
})


test_that("all-zero sample columns produce NA for applicable indices", {
  fx <- make_fixture()

  res <- divermeta(fx$abund, diss = fx$diss, indices = c("raoQ", "FD_sigma", "multiplicity_inventory"), clusters = fx$clusters, sig = 0.8)

  # Third sample is all zeros
  expect_true(is.na(res$raoQ[3]))
  expect_true(is.na(res$FD_sigma[3]))
  expect_true(is.na(res$multiplicity_inventory[3]))
})
