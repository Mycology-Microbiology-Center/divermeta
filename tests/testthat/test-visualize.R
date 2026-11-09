## Tests for the plotting functions

skip_if_not_installed("ggplot2")


test_that("visualize_abundance draws every value in its sample and subunit", {
  abund <- matrix(
    c(
      1, 0, 3, 4, 5,
      6, 7, 0, 9, 10,
      11, 12, 13, 0, 15
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), paste0("f", 1:5))
  )
  p <- visualize_abundance(abund)
  expect_s3_class(p, "ggplot")

  dat <- p$data
  expect_equal(nrow(dat), length(abund))
  expect_equal(dat$value, abund[cbind(as.character(dat$row), as.character(dat$col))])
  # First sample on top, subunits in column order
  expect_identical(levels(dat$row), rev(rownames(abund)))
  expect_identical(levels(dat$col), colnames(abund))

  # Same data from a data frame and a sparse matrix
  expect_equal(visualize_abundance(as.data.frame(abund))$data, dat)
  expect_equal(visualize_abundance(Matrix::Matrix(abund, sparse = TRUE))$data, dat)

  # Transformed sizes
  expect_equal(visualize_abundance(abund, transform = "sqrt")$data$value, sqrt(dat$value))
})


test_that("visualize_abundance labels", {
  abund <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("a", "b", "c")))
  dat <- visualize_abundance(abund)$data
  expect_identical(levels(dat$row), c("2", "1"))

  dat <- visualize_abundance(abund, sample.labels = c("x", "y"), subunit.labels = c("A", "B", "C"))$data
  expect_identical(levels(dat$row), c("y", "x"))
  expect_identical(levels(dat$col), c("A", "B", "C"))
  expect_equal(dat$value[dat$row == "y" & dat$col == "C"], 6)

  expect_error(visualize_abundance(abund, sample.labels = "x"), "one label")
  expect_error(visualize_abundance(abund, subunit.labels = c("A", NA, "C")), "one label")
  expect_error(visualize_abundance(unname(abund)), "column names")
  expect_error(visualize_abundance(abund * -1), "non-negative")
})


test_that("visualize_dist draws the table as a symmetric matrix", {
  diss <- data.frame(
    ID1 = c("a", "a", "b", "c"),
    ID2 = c("b", "c", "c", "d"),
    Distance = c(0.2, 0.7, 0.6, 0.4)
  )
  p <- visualize_dist(diss)
  expect_s3_class(p, "ggplot")

  dat <- p$data
  value <- function(r, c) dat$mag[dat$row == r & dat$col == c]
  expect_identical(levels(dat$row), rev(c("a", "b", "c", "d")))
  expect_identical(levels(dat$col), c("a", "b", "c", "d"))
  for (k in seq_len(nrow(diss))) {
    expect_equal(value(diss$ID1[k], diss$ID2[k]), diss$Distance[k])
    expect_equal(value(diss$ID2[k], diss$ID1[k]), diss$Distance[k])
  }
  expect_equal(value("b", "b"), 0)
  # Pairs not listed are not drawn
  expect_length(value("a", "d"), 0)
  expect_equal(nrow(dat), 4 + 2 * nrow(diss))

  # Order and subset of the subunits, and labels
  dat <- visualize_dist(diss, ids = c("c", "a"), labels = c("C", "A"))$data
  expect_identical(levels(dat$col), c("C", "A"))
  expect_equal(dat$mag[dat$row == "C" & dat$col == "A"], 0.7)
  expect_equal(nrow(dat), 4)
})


test_that("visualize_dist validates its input", {
  diss <- data.frame(ID1 = c("a", "a"), ID2 = c("b", "c"), Distance = c(0.2, 0.7))
  expect_error(visualize_dist(write_diss(diss, "csv")), "not a file")
  expect_error(visualize_dist(diss[, 1:2]), "three columns")
  expect_error(visualize_dist(transform(diss, Distance = c(-1, 0.7))), "non-negative")
  expect_error(visualize_dist(rbind(diss, data.frame(ID1 = "b", ID2 = "a", Distance = 0.9))), "conflicting")
  expect_error(visualize_dist(diss, labels = "x"), "one label")
  expect_error(visualize_dist(diss, ids = c("a", "a")), "unique")
})
