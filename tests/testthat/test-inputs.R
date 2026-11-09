## Tests for the validation of the shared inputs


# Small study: two units of two subunits and an empty sample
input_study <- function() {
  ids <- c("a1", "a2", "b1", "b2")
  abund <- matrix(
    c(
      4, 1, 2, 3,
      0, 5, 1, 0,
      0, 0, 0, 0
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), ids)
  )
  set.seed(7)
  list(
    abund = abund,
    clust = c(a1 = "A", a2 = "A", b1 = "B", b2 = "B"),
    diss = mat_to_table(rand_diss(4), ids)
  )
}


# Non-finite values
# ----------------------------

test_that("infinite abundances are errors", {
  st <- input_study()
  bad <- st$abund
  bad[1, 2] <- Inf

  expect_error(raoQuadratic(bad, st$diss), "`abund` contains infinite values")
  expect_error(multiplicity.inventory(bad, st$clust), "`abund` contains infinite values")
  expect_error(divermeta(bad, clust = st$clust), "`abund` contains infinite values")
  expect_error(
    raoQuadratic(Matrix::Matrix(bad, sparse = TRUE), st$diss),
    "`abund` contains infinite values"
  )
})


test_that("infinite distances are errors, in memory and in files", {
  st <- input_study()
  bad <- st$diss
  bad$Distance[2] <- Inf

  expect_error(raoQuadratic(st$abund, bad), "`diss` contains infinite distances")
  expect_error(multiplicity.distance(st$abund, bad, st$clust, method = "min"), "infinite distances")
  expect_error(unchecked(raoQuadratic(st$abund, write_diss(bad, "csv"))), "infinite distances")
  expect_error(
    raoQuadratic(st$abund, write_diss(bad, "csv"), check_large_distance_file = TRUE),
    "did not pass the check"
  )

  diss_clust <- data.frame(ID1 = "A", ID2 = "B", Distance = Inf)
  expect_error(
    multiplicity.distance(st$abund, st$diss, st$clust, method = "custom", diss_clust = diss_clust),
    "`diss_clust` contains infinite distances"
  )
})


# Named clust
# ----------------------------

test_that("subunits of a named clust that are not in abund are dropped", {
  st <- input_study()
  extra <- c(st$clust, z1 = "A", z2 = "C")

  expect_equal(multiplicity.inventory(st$abund, extra), multiplicity.inventory(st$abund, st$clust))
  for (method in c("sigma", "average")) {
    expect_equal(
      multiplicity.distance(st$abund, st$diss, extra, method = method),
      multiplicity.distance(st$abund, st$diss, st$clust, method = method),
      info = method
    )
  }
  expect_equal(
    metagenomic.alpha.index(st$abund, st$diss, extra),
    metagenomic.alpha.index(st$abund, st$diss, st$clust)
  )
  indices <- c("multiplicity_inventory", "multiplicity_distance", "raoQ")
  expect_equal(
    divermeta(st$abund, st$diss, indices, extra),
    divermeta(st$abund, st$diss, indices, st$clust)
  )
})


# No pair needed
# ----------------------------

test_that("a single subunit needs no pair, in memory and in files", {
  ab <- one_row(5, "a")
  # Rows of other subunits only
  diss <- data.frame(ID1 = c("x", "x"), ID2 = c("y", "z"), Distance = c(0.3, 0.6))
  path <- write_diss(diss, "csv")

  expect_equal(unname(raoQuadratic(ab, diss)), 0)
  expect_equal(unname(unchecked(raoQuadratic(ab, path))), 0)
  expect_equal(unname(raoQuadratic(ab, path, check_large_distance_file = TRUE)), 0)
  expect_equal(unname(unchecked(multiplicity.distance(ab, path, c(a = "A"), method = "min"))), 1)

  # A file of other subunits is still an error when pairs are needed
  ab2 <- one_row(c(5, 1), c("a", "b"))
  expect_error(unchecked(raoQuadratic(ab2, path)), "No row of the distance file pairs two subunits")
})
