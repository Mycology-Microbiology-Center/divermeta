## Tests for Rao quadratic entropy

test_that("raoQuadratic deterministic numeric (2 species)", {
  ab <- c(2, 1)
  d <- 0.4
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = d)

  # Q = 2 * p1 * p2 * d
  p <- ab / sum(ab)
  expected <- 2 * p[1] * p[2] * d

  expect_equal(unname(raoQuadratic(one_row(ab), diss)), expected, tolerance = 1e-12)
})


test_that("raoQuadratic permutation invariance", {
  ab <- c(2, 3, 5)
  diss <- matrix(c(
    0,   0.1, 0.7,
    0.1, 0,   0.3,
    0.7, 0.3, 0
  ), nrow = 3, byrow = TRUE)
  ids <- c("a", "b", "c")

  baseline <- raoQuadratic(one_row(ab, ids), mat_to_table(diss, ids))

  # Permute order consistently
  idx <- c(3, 1, 2)
  ab2 <- ab[idx]
  diss2 <- diss[idx, idx]

  expect_equal(raoQuadratic(one_row(ab2, ids[idx]), mat_to_table(diss2, ids[idx])), baseline, tolerance = 1e-12)
})


test_that("raoQuadratic input validation", {
  diss <- data.frame(ID1 = "1", ID2 = "2", Distance = 0.2)
  expect_error(raoQuadratic(one_row(c(1, NA)), diss), "NA")

  expect_error(raoQuadratic(one_row(c(-1, 2)), diss), "non-negative")

  # Empty samples give NA
  expect_true(is.na(raoQuadratic(one_row(c(0, 0)), diss)))

  # Malformed distance tables
  ab <- one_row(c(1, 2))
  expect_error(raoQuadratic(ab, data.frame(ID1 = "1", ID2 = "2")), "three columns")
  expect_error(raoQuadratic(ab, list(1, 2, 3)), "three columns")
  # A string is the path to a distance file
  expect_error(raoQuadratic(ab, "not a table"), "not found")
  expect_error(raoQuadratic(ab, data.frame(ID1 = "1", ID2 = "2", Distance = "a")), "numeric")
  expect_error(raoQuadratic(ab, data.frame(ID1 = "1", ID2 = "2", Distance = NA_real_)), "NA")
  expect_error(raoQuadratic(ab, data.frame(ID1 = "1", ID2 = "2", Distance = -0.1)), "non-negative")

  # Abundance table without subunit identifiers
  expect_error(raoQuadratic(matrix(c(1, 2), nrow = 1), diss), "column names")
})


test_that("raoQuadratic equals the legacy implementation for many samples", {
  set.seed(42)
  n <- 30
  ids <- paste0("f", seq_len(n))
  diss <- as.matrix(dist(rnorm(n)))
  abund <- matrix(
    runif(8 * n, min = 10, max = 100) * (runif(8 * n) > 0.3),
    nrow = 8,
    dimnames = list(paste0("S", 1:8), ids)
  )

  res <- raoQuadratic(abund, mat_to_table(diss, ids))
  expected <- apply(abund, 1, function(ab) legacy$raoQuadratic(ab, diss))

  expect_equal(res, expected)
  expect_identical(names(res), rownames(abund))

  # Sparse abundances
  expect_equal(raoQuadratic(Matrix::Matrix(abund, sparse = TRUE), mat_to_table(diss, ids)), res)
})


test_that("Distance table orientation, order, duplicates and extra rows", {
  set.seed(1)
  n <- 6
  ids <- letters[1:n]
  diss <- rand_diss(n)
  abund <- matrix(runif(3 * n, 1, 10), nrow = 3, dimnames = list(NULL, ids))
  tab <- mat_to_table(diss, ids)
  base <- raoQuadratic(abund, tab)

  # Random orientation and row order
  swap <- runif(nrow(tab)) > 0.5
  mixed <- tab
  mixed[swap, c("ID1", "ID2")] <- mixed[swap, c("ID2", "ID1")]
  mixed <- mixed[sample(nrow(mixed)), ]
  expect_equal(raoQuadratic(abund, mixed), base)

  # Both orientations with the same distance
  flipped <- tab
  flipped[, c("ID1", "ID2")] <- flipped[, c("ID2", "ID1")]
  expect_equal(raoQuadratic(abund, rbind(tab, flipped)), base)

  # Self pairs and unknown identifiers are ignored
  extra <- data.frame(ID1 = c("a", "zz"), ID2 = c("a", "b"), Distance = c(0, 0.3))
  expect_equal(raoQuadratic(abund, rbind(tab, extra)), base)

  # Conflicting distances for the same pair
  conflict <- tab[1, ]
  conflict$Distance <- conflict$Distance + 0.1
  expect_error(raoQuadratic(abund, rbind(tab, conflict)), "conflicting")

  # Missing pair
  expect_error(raoQuadratic(abund, tab[-3, ]), "Missing distances.*1 pair")

  # Numeric identifiers
  num_abund <- abund
  colnames(num_abund) <- sprintf("%.0f", (1:6) * 1e5)
  num_tab <- tab
  num_tab$ID1 <- match(tab$ID1, ids) * 1e5
  num_tab$ID2 <- match(tab$ID2, ids) * 1e5
  expect_equal(unname(raoQuadratic(num_abund, num_tab)), unname(base))
})


test_that("Numeric identifiers of 1e5 or more match the column names", {
  set.seed(42)
  num_ids <- c(99999, 100000, 100001)
  abund <- matrix(runif(6, 1, 10), nrow = 2, dimnames = list(c("S1", "S2"), NULL))
  # colnames() stores 100000 as "1e+05"
  colnames(abund) <- num_ids
  diss_num <- data.frame(ID1 = num_ids[c(1, 1, 2)], ID2 = num_ids[c(2, 3, 3)], Distance = c(0.2, 0.5, 0.7))

  abund_chr <- abund
  colnames(abund_chr) <- c("x", "y", "z")
  diss_chr <- data.frame(ID1 = c("x", "x", "y"), ID2 = c("y", "z", "z"), Distance = diss_num$Distance)

  expect_identical(colnames(abund)[2], "1e+05")
  expect_equal(raoQuadratic(abund, diss_num), raoQuadratic(abund_chr, diss_chr))
})
