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
manual_pa <- function(abund, percentile = 100, cap_at = NULL) {
  ref <- apply(abund, 2, function(x) {
    if (all(x == 0)) NA_real_ else unname(stats::quantile(x[x > 0], percentile / 100, type = 7))
  })
  vapply(seq_len(nrow(abund)), function(i) {
    if (sum(abund[i, ]) == 0) {
      return(0)
    }
    present <- abund[i, ] > 0
    ratio <- abund[i, present] / ref[present]
    if (!is.null(cap_at)) ratio <- pmin(ratio, cap_at)
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
  expect_equal(proportional.abundance(abund, percentile = 50, cap_at = 1), expected_cap)

  # Capped at another value: S3 = (1.2, 1.5) and S4 = (4, 1) become (1.2, 1.5) and (1.5, 1)
  expected_cap <- c(S1 = 1 / 2.5 + 2 / 4, S2 = 2 / 2.5, S3 = 3 / 2.5 + 1.5, S4 = 1.5 + 1)
  expect_equal(proportional.abundance(abund, percentile = 50, cap_at = 1.5), expected_cap)

  # Below one, also the ratios under 1 are capped
  expected_cap <- c(S1 = 0.4 + 0.5, S2 = 0.5, S3 = 0.5 + 0.5, S4 = 0.5 + 0.5)
  expect_equal(proportional.abundance(abund, percentile = 50, cap_at = 0.5), expected_cap)
})


test_that("Equals the naive implementation on random tables", {
  set.seed(42)
  for (w_ in 1:10) {
    abund <- make_table(sample(2:12, 1), sample(1:9, 1))
    for (percentile in c(100, 95, 75, 50, 10, 0.5)) {
      for (cap_at in list(NULL, 1, 0.5, 3)) {
        expect_equal(
          unname(proportional.abundance(abund, percentile, cap_at)),
          manual_pa(abund, percentile, cap_at),
          tolerance = 1e-12,
          info = paste(w_, percentile, format(cap_at))
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
  capped <- proportional.abundance(abund, percentile = 20, cap_at = 1)
  expect_true(all(capped <= n_present + 1e-12))
  expect_true(all(capped <= low + 1e-12))

  # Capping does nothing with the maximum as reference
  expect_equal(proportional.abundance(abund, cap_at = 1), res)
  expect_equal(proportional.abundance(abund, cap_at = 7), res)
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
  expect_identical(names(formals(proportional.abundance)), c("abund", "percentile", "cap_at"))

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

  for (bad in list(NA, "yes", TRUE, 0, -1, Inf, NaN, c(1, 2), numeric(0))) {
    expect_error(proportional.abundance(abund, cap_at = bad), "`cap_at` must be NULL", info = format(bad))
  }
  expect_no_error(proportional.abundance(abund, cap_at = NULL))
  expect_no_error(proportional.abundance(abund, cap_at = 0.1))

  expect_error(proportional.abundance(one_row(c(-1, 2, 3), paste0("e", 1:3))), "non-negative")
  expect_error(proportional.abundance(one_row(c(NA, 2, 3), paste0("e", 1:3))), "NA")
  expect_error(proportional.abundance(one_row(c(Inf, 2, 3), paste0("e", 1:3))), "infinite")
  expect_error(proportional.abundance("a"), "must be a matrix")
  expect_error(proportional.abundance(matrix(1:6, 2)), "column names")
})


# Known reference abundances
# ----------------------------

# Naive implementation with given references: a present element with a zero
# reference contributes cap_at
manual_pa_ref <- function(abund, ref_abu, cap_at = NULL) {
  vapply(seq_len(nrow(abund)), function(i) {
    total <- 0
    for (j in seq_len(ncol(abund))) {
      x <- abund[i, j]
      if (x == 0) next
      ratio <- if (ref_abu[j] == 0) Inf else x / ref_abu[j]
      if (!is.null(cap_at)) ratio <- min(ratio, cap_at)
      total <- total + ratio
    }
    total
  }, numeric(1))
}

pa_ref_table <- function() {
  matrix(
    c(
      10, 0, 4,
      5, 8, 2,
      0, 2, 8
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), c("e1", "e2", "e3"))
  )
}


test_that("ref_abu: formula by hand", {
  abund <- pa_ref_table()
  ref_abu <- c(e1 = 5, e2 = 4, e3 = 8)

  # Ratios: S1 = (2, -, 0.5), S2 = (1, 2, 0.25), S3 = (-, 0.5, 1)
  expected <- c(S1 = 2 + 0.5, S2 = 1 + 2 + 0.25, S3 = 0.5 + 1)
  expect_equal(proportional.abundance.ref_abu(abund, ref_abu), expected)

  expected_cap <- c(S1 = 1 + 0.5, S2 = 1 + 1 + 0.25, S3 = 0.5 + 1)
  expect_equal(proportional.abundance.ref_abu(abund, ref_abu, cap_at = 1), expected_cap)

  expected_cap <- c(S1 = 0.5 + 0.5, S2 = 0.5 + 0.5 + 0.25, S3 = 0.5 + 0.5)
  expect_equal(proportional.abundance.ref_abu(abund, ref_abu, cap_at = 0.5), expected_cap)
})


test_that("ref_abu: the percentiles as references give proportional.abundance", {
  set.seed(31)
  for (w_ in 1:8) {
    abund <- make_table(sample(2:12, 1), sample(1:9, 1))
    for (percentile in c(100, 50, 10)) {
      # Elements absent from every sample contribute 0 with any reference
      ref_abu <- apply(abund, 2, function(x) {
        if (all(x == 0)) 1 else unname(stats::quantile(x[x > 0], percentile / 100, type = 7))
      })
      for (cap_at in list(NULL, 1, 0.5)) {
        expect_equal(
          proportional.abundance.ref_abu(abund, ref_abu, cap_at),
          proportional.abundance(abund, percentile, cap_at),
          tolerance = 1e-12,
          info = paste(w_, percentile, format(cap_at))
        )
      }
    }
  }
})


test_that("ref_abu: equals the naive implementation on random tables", {
  set.seed(32)
  for (w_ in 1:10) {
    abund <- make_table(sample(2:12, 1), sample(1:9, 1))
    ref_abu <- runif(ncol(abund), 0.5, 120)
    for (cap_at in list(NULL, 1, 0.5, 3)) {
      expect_equal(
        unname(proportional.abundance.ref_abu(abund, ref_abu, cap_at)),
        manual_pa_ref(abund, ref_abu, cap_at),
        tolerance = 1e-12,
        info = paste(w_, format(cap_at))
      )
    }

    # Some zero references, allowed with cap_at
    ref_abu[sample(length(ref_abu), ceiling(length(ref_abu) / 2))] <- 0
    for (cap_at in list(1, 0.5, 3)) {
      res <- proportional.abundance.ref_abu(abund, ref_abu, cap_at)
      expect_true(all(is.finite(res)))
      expect_equal(unname(res), manual_pa_ref(abund, ref_abu, cap_at), tolerance = 1e-12,
                   info = paste("zeros", w_, cap_at))
    }
  }
})


test_that("ref_abu: zero references give cap_at where the element is present", {
  abund <- pa_ref_table()

  # e2 has no reference: present in S2 and S3, absent from S1
  ref_abu <- c(e1 = 5, e2 = 0, e3 = 8)
  expect_equal(
    proportional.abundance.ref_abu(abund, ref_abu, cap_at = 1.5),
    c(S1 = 1.5 + 0.5, S2 = 1 + 1.5 + 0.25, S3 = 1.5 + 1)
  )
  expect_equal(
    proportional.abundance.ref_abu(abund, ref_abu, cap_at = 0.4),
    c(S1 = 0.4 + 0.4, S2 = 0.4 + 0.4 + 0.25, S3 = 0.4 + 0.4)
  )

  # Every reference zero: cap_at times the number of elements present
  for (cap_at in c(0.3, 1, 2)) {
    expect_equal(
      proportional.abundance.ref_abu(abund, c(0, 0, 0), cap_at),
      cap_at * rowSums(abund > 0)
    )
  }

  # Empty samples are 0, also with zero references
  empty <- rbind(abund, Z = 0)
  res <- proportional.abundance.ref_abu(empty, ref_abu, cap_at = 1)
  expect_equal(unname(res["Z"]), 0)
  expect_false(anyNA(res))
  expect_equal(unname(proportional.abundance.ref_abu(abund * 0, c(0, 0, 0), cap_at = 1)), rep(0, 3))
})


test_that("ref_abu: matching the references to the elements", {
  abund <- pa_ref_table()
  ref_abu <- c(e1 = 5, e2 = 4, e3 = 8)
  base <- proportional.abundance.ref_abu(abund, ref_abu)

  # Unnamed values follow the columns
  expect_equal(proportional.abundance.ref_abu(abund, unname(ref_abu)), base)

  # Named in any order, extra names ignored
  expect_equal(proportional.abundance.ref_abu(abund, ref_abu[c(3, 1, 2)]), base)
  expect_equal(proportional.abundance.ref_abu(abund, c(ref_abu, x9 = 0.001)), base)

  # Columns in any order with named references
  expect_equal(proportional.abundance.ref_abu(abund[, c(3, 1, 2)], ref_abu), base)

  expect_error(proportional.abundance.ref_abu(abund, ref_abu[1:2]), "missing values for elements: e3")
  expect_error(proportional.abundance.ref_abu(abund, c(5, 4)), "one value per element")
  expect_error(proportional.abundance.ref_abu(abund, c(5, 4, 8, 1)), "one value per element")
})


test_that("ref_abu: formats, scale and absent elements", {
  set.seed(33)
  abund <- make_table(8, 6)
  ref_abu <- runif(6, 1, 100)
  base <- proportional.abundance.ref_abu(abund, ref_abu, cap_at = 2)

  expect_equal(proportional.abundance.ref_abu(as.data.frame(abund), ref_abu, cap_at = 2), base)
  expect_equal(
    proportional.abundance.ref_abu(Matrix::Matrix(abund, sparse = TRUE), ref_abu, cap_at = 2),
    base
  )
  expect_identical(names(base), rownames(abund))

  # Same constant on the abundances and the references
  expect_equal(proportional.abundance.ref_abu(abund * 7.3, ref_abu * 7.3, cap_at = 2), base)

  # Elements absent from every sample contribute 0, with any reference
  extended <- cbind(abund, x1 = 0, x2 = 0)
  expect_equal(proportional.abundance.ref_abu(extended, c(ref_abu, 3, 0), cap_at = 2), base)
})


test_that("ref_abu: input validation", {
  abund <- one_row(c(1, 2, 3), paste0("e", 1:3))
  pa <- function(ref_abu, ...) proportional.abundance.ref_abu(abund, ref_abu, ...)

  expect_identical(
    names(formals(proportional.abundance.ref_abu)), c("abund", "ref_abu", "cap_at")
  )

  # Zero references only with cap_at
  expect_error(pa(c(1, 0, 2)), "strictly positive.*allowed only with `cap_at`")
  expect_no_error(pa(c(1, 0, 2), cap_at = 1))
  expect_error(pa(c(1, -1, 2), cap_at = 1), "non-negative")
  expect_error(pa(c(1, -1, 2)), "strictly positive")

  for (bad in list(c(1, NA, 2), c(1, Inf, 2), c(1, NaN, 2), "a", numeric(0), NULL, TRUE)) {
    expect_error(pa(bad), "`ref_abu` must contain", info = format(bad))
    expect_error(pa(bad, cap_at = 1), "`ref_abu` must contain", info = format(bad))
  }
  expect_error(proportional.abundance.ref_abu(abund), "ref_abu")

  for (bad in list(NA, "yes", TRUE, 0, -1, Inf, NaN, c(1, 2), numeric(0))) {
    expect_error(pa(c(1, 2, 3), cap_at = bad), "`cap_at` must be NULL", info = format(bad))
  }

  expect_error(proportional.abundance.ref_abu(one_row(c(-1, 2, 3), paste0("e", 1:3)), c(1, 2, 3)), "non-negative")
  expect_error(proportional.abundance.ref_abu("a", 1), "must be a matrix")
  expect_error(proportional.abundance.ref_abu(matrix(1:6, 2), c(1, 2, 3)), "column names")
})
