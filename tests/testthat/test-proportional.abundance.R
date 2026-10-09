# Tests for proportional abundance

# Helpers
# ----------------------------

# Random abundance table (samples in rows) with zeros
make_table <- function(n_samples, n_elements, zero_prob = 0.4) {
  matrix(
    runif(n_samples * n_elements, 1, 100) * (runif(n_samples * n_elements) > zero_prob),
    nrow = n_samples,
    dimnames = list(paste0("S", seq_len(n_samples)), paste0("e", seq_len(n_elements)))
  )
}

# Naive implementation: element by element and sample by sample
manual_pa <- function(abund, percentile = 100, cap_at_one = FALSE) {
  ref <- apply(abund, 2, function(x) {
    if (all(x == 0)) NA_real_ else unname(stats::quantile(x[x > 0], percentile / 100, type = 7))
  })
  vapply(seq_len(nrow(abund)), function(i) {
    if (sum(abund[i, ]) == 0) {
      return(0)
    }
    present <- abund[i, ] > 0
    ratio <- abund[i, present] / ref[present]
    if (cap_at_one) ratio <- pmin(ratio, 1)
    sum(ratio)
  }, numeric(1))
}


# Tests
# ----------------------------

test_that("Formula by hand with the maximum as reference", {
  abund <- matrix(
    c(
      10, 0, 4,
      5, 8, 2,
      0, 2, 8
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), c("e1", "e2", "e3"))
  )
  # References: e1 = 10, e2 = 8, e3 = 8
  expected <- c(S1 = 10 / 10 + 4 / 8, S2 = 5 / 10 + 8 / 8 + 2 / 8, S3 = 2 / 8 + 8 / 8)
  expect_equal(proportional.abundance(abund), expected)
  expect_equal(proportional.abundance(abund, percentile = 100), expected)
})


test_that("Formula by hand with a lower percentile", {
  abund <- matrix(
    c(
      1, 2,
      2, 0,
      3, 6,
      10, 4
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("S", 1:4), c("e1", "e2"))
  )
  # Median of the nonzero abundances: e1 = median(1, 2, 3, 10) = 2.5, e2 = median(2, 6, 4) = 4
  expected <- c(
    S1 = 1 / 2.5 + 2 / 4,
    S2 = 2 / 2.5,
    S3 = 3 / 2.5 + 6 / 4,
    S4 = 10 / 2.5 + 4 / 4
  )
  expect_equal(proportional.abundance(abund, percentile = 50), expected)

  # Capped at one: the ratios above 1 become 1
  expected_cap <- c(S1 = 1 / 2.5 + 2 / 4, S2 = 2 / 2.5, S3 = 1 + 1, S4 = 1 + 1)
  expect_equal(proportional.abundance(abund, percentile = 50, cap_at_one = TRUE), expected_cap)
})


test_that("Equals the naive implementation on random tables", {
  set.seed(42)
  for (w_ in 1:10) {
    abund <- make_table(sample(2:12, 1), sample(1:9, 1))
    for (percentile in c(100, 95, 75, 50, 10, 0.5)) {
      for (cap_at_one in c(FALSE, TRUE)) {
        expect_equal(
          unname(proportional.abundance(abund, percentile, cap_at_one)),
          manual_pa(abund, percentile, cap_at_one),
          tolerance = 1e-12,
          info = paste(w_, percentile, cap_at_one)
        )
      }
    }
  }
})


test_that("Range of the values", {
  set.seed(7)
  abund <- make_table(10, 6, zero_prob = 0.3)
  n_present <- rowSums(abund > 0)

  # The maximum as reference: every ratio is in (0, 1]
  res <- proportional.abundance(abund)
  expect_true(all(res > 0 & res <= n_present + 1e-12, na.rm = TRUE))

  # Lower percentiles can exceed the number of elements present, unless capped
  low <- proportional.abundance(abund, percentile = 20)
  expect_true(any(low > n_present))
  capped <- proportional.abundance(abund, percentile = 20, cap_at_one = TRUE)
  expect_true(all(capped <= n_present + 1e-12))
  expect_true(all(capped <= low + 1e-12))

  # Capping does nothing with the maximum as reference
  expect_equal(proportional.abundance(abund, cap_at_one = TRUE), res)
})


test_that("Not increasing in the percentile", {
  set.seed(3)
  abund <- make_table(9, 5)
  percentiles <- c(5, 25, 50, 75, 90, 100)
  values <- vapply(percentiles, function(p) unname(proportional.abundance(abund, p)), numeric(nrow(abund)))
  expect_true(all(diff(t(values)) <= 1e-12, na.rm = TRUE))
})


test_that("Zeros are not part of the reference", {
  set.seed(11)
  abund <- make_table(6, 4, zero_prob = 0.5)
  base <- proportional.abundance(abund, percentile = 40)

  # Extra empty samples change neither the reference nor the other samples
  extended <- rbind(abund, Z1 = 0, Z2 = 0)
  expect_equal(proportional.abundance(extended, percentile = 40)[rownames(abund)], base)

  # A column with a lot of zeros has the same reference as without them
  sparse_col <- matrix(c(0, 0, 0, 0, 5, 0, 0, 10), ncol = 1, dimnames = list(paste0("S", 1:8), "e1"))
  expect_equal(
    unname(proportional.abundance(sparse_col, percentile = 100))[c(5, 8)],
    c(0.5, 1)
  )
  expect_equal(
    unname(proportional.abundance(sparse_col, percentile = 50))[c(5, 8)],
    c(5 / 7.5, 10 / 7.5)
  )
})


test_that("Elements absent from every sample contribute 0", {
  set.seed(5)
  abund <- make_table(5, 4, zero_prob = 0.2)
  base <- proportional.abundance(abund, percentile = 60)

  extended <- cbind(abund, x1 = 0, x2 = 0)
  res <- proportional.abundance(extended, percentile = 60)
  expect_false(anyNA(res))
  expect_true(all(is.finite(res)))
  expect_equal(res, base)
})


test_that("A single sample counts the elements present", {
  abund <- one_row(c(3, 0, 7, 1, 0), paste0("e", 1:5))
  expect_equal(unname(proportional.abundance(abund)), 3)
  expect_equal(unname(proportional.abundance(abund, percentile = 30)), 3)

  # A single element
  col <- matrix(c(2, 4, 0, 8), ncol = 1, dimnames = list(paste0("S", 1:4), "e1"))
  expect_equal(unname(proportional.abundance(col)), c(0.25, 0.5, 0, 1))
})


test_that("Elements with all the same abundance have ratio 1", {
  abund <- matrix(
    c(
      4, 0, 3,
      4, 5, 0,
      4, 5, 3
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("S", 1:3), c("e1", "e2", "e3"))
  )
  for (percentile in c(100, 50, 1)) {
    expect_equal(unname(proportional.abundance(abund, percentile)), c(2, 2, 3))
  }
})


test_that("Empty samples are 0", {
  set.seed(2)
  abund <- make_table(5, 4, zero_prob = 0.2)
  abund["S3", ] <- 0
  res <- proportional.abundance(abund, percentile = 70)
  expect_equal(unname(res["S3"]), 0)
  expect_false(anyNA(res[names(res) != "S3"]))

  # The empty sample does not change the others
  expect_equal(res[-3], proportional.abundance(abund[-3, ], percentile = 70))

  # All samples empty
  expect_equal(unname(proportional.abundance(abund * 0)), rep(0, 5))
})


test_that("Invariant to the scale of the abundances and to the order", {
  set.seed(8)
  abund <- make_table(7, 5)
  for (percentile in c(100, 35)) {
    base <- proportional.abundance(abund, percentile)
    expect_equal(proportional.abundance(abund * 13.7, percentile), base)

    # Columns in any order
    expect_equal(proportional.abundance(abund[, rev(colnames(abund))], percentile), base)

    # Rows in any order, keeping the names
    perm <- sample(nrow(abund))
    expect_equal(proportional.abundance(abund[perm, ], percentile), base[perm])
  }
})


test_that("Matrix, data.frame and sparse abundance formats give the same result", {
  set.seed(21)
  abund <- make_table(8, 6)
  base <- proportional.abundance(abund, percentile = 60)

  expect_equal(proportional.abundance(as.data.frame(abund), percentile = 60), base)
  expect_equal(proportional.abundance(Matrix::Matrix(abund, sparse = TRUE), percentile = 60), base)

  # Named by the samples, or by their numbers
  expect_identical(names(base), rownames(abund))
  no_names <- abund
  rownames(no_names) <- NULL
  expect_identical(names(proportional.abundance(no_names)), as.character(1:8))
})


test_that("Does not depend on the elements being units or subunits", {
  # The function receives only the abundances and the percentile
  expect_identical(names(formals(proportional.abundance)), c("abund", "percentile", "cap_at_one"))

  set.seed(4)
  abund <- make_table(6, 8)
  clust <- rep(c("A", "B", "C", "D"), each = 2)
  units <- t(rowsum(t(abund), clust))
  expect_equal(unname(proportional.abundance(units, percentile = 80)), manual_pa(units, 80))
})


test_that("Input validation", {
  abund <- one_row(c(1, 2, 3), paste0("e", 1:3))

  for (bad in list(0, -5, 100.1, 101, NA_real_, Inf, -Inf, NaN, "a", TRUE, c(50, 60), numeric(0), NULL)) {
    expect_error(proportional.abundance(abund, percentile = bad), "percentile", info = format(bad))
  }
  expect_no_error(proportional.abundance(abund, percentile = 100))
  expect_no_error(proportional.abundance(abund, percentile = 1e-6))

  expect_error(proportional.abundance(abund, cap_at_one = NA), "single TRUE or FALSE")
  expect_error(proportional.abundance(abund, cap_at_one = "yes"), "single TRUE or FALSE")
  expect_error(proportional.abundance(abund, cap_at_one = 1), "single TRUE or FALSE")

  expect_error(proportional.abundance(one_row(c(-1, 2, 3), paste0("e", 1:3))), "non-negative")
  expect_error(proportional.abundance(one_row(c(NA, 2, 3), paste0("e", 1:3))), "NA")
  expect_error(proportional.abundance(one_row(c(Inf, 2, 3), paste0("e", 1:3))), "infinite")
  expect_error(proportional.abundance("a"), "must be a matrix")
  expect_error(proportional.abundance(matrix(1:6, 2)), "column names")
})
