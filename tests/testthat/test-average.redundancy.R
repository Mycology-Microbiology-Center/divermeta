# Tests for average redundancy

# Helpers
# ----------------------------

# Random study with units of the given sizes (samples in rows, subunits in
# columns) and only the distances within units
make_units <- function(sizes, n_samples, zero_prob = 0.4) {
  unit_ids <- paste0("U", seq_along(sizes))
  sub_ids <- lapply(seq_along(sizes), function(u) paste0(unit_ids[u], "_", seq_len(sizes[u])))
  clust <- stats::setNames(rep(unit_ids, sizes), unlist(sub_ids))

  mats <- lapply(sizes, rand_diss)
  diss_frame <- do.call(rbind, lapply(seq_along(sizes), function(u) mat_to_table(mats[[u]], sub_ids[[u]])))

  n <- sum(sizes)
  abund <- t(matrix(
    runif(n * n_samples, 1, 100) * (runif(n * n_samples) > zero_prob),
    nrow = n,
    dimnames = list(names(clust), paste0("S", seq_len(n_samples)))
  ))

  list(clust = clust, unit_ids = unit_ids, sub_ids = sub_ids, mats = mats, diss_frame = diss_frame, abund = abund)
}

# Average redundancy of every sample from redundancy() of every unit alone
manual_ar <- function(study, normalize = FALSE, include_absent = FALSE) {
  vapply(rownames(study$abund), function(s) {
    per_unit <- vapply(seq_along(study$sub_ids), function(u) {
      ids <- study$sub_ids[[u]]
      ab <- study$abund[s, ids]
      if (sum(ab) == 0) {
        return(NA_real_)
      }
      unname(redundancy(one_row(ab, ids), mat_to_table(study$mats[[u]], ids), normalize = normalize))
    }, numeric(1))
    present <- !is.na(per_unit)
    if (!any(present)) {
      return(NA_real_)
    }
    if (include_absent) {
      sum(per_unit[present]) / length(per_unit)
    } else {
      mean(per_unit[present])
    }
  }, numeric(1))
}


# Tests
# ----------------------------

test_that("Equals the average of the redundancy of every unit", {
  set.seed(42)
  for (w_ in 1:8) {
    study <- make_units(sample(1:8, sample(2:5, 1), replace = TRUE), n_samples = 6)
    for (normalize in c(FALSE, TRUE)) {
      for (include_absent in c(FALSE, TRUE)) {
        expect_equal(
          average.redundancy(
            study$abund, study$diss_frame, study$clust,
            normalize = normalize, include_absent = include_absent
          ),
          manual_ar(study, normalize, include_absent),
          tolerance = 1e-12,
          info = paste(w_, normalize, include_absent)
        )
      }
    }
  }
})


test_that("A single unit equals redundancy", {
  set.seed(1)
  study <- make_units(7, n_samples = 5, zero_prob = 0.3)
  for (normalize in c(FALSE, TRUE)) {
    expect_equal(
      average.redundancy(study$abund, study$diss_frame, study$clust, normalize = normalize),
      redundancy(study$abund, study$diss_frame, normalize = normalize),
      tolerance = 1e-12
    )
  }
})


test_that("Closed form for two units of two subunits", {
  clust <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
  diss <- data.frame(ID1 = c("a1", "b1"), ID2 = c("a2", "b2"), Distance = c(0.2, 0.7))
  abund <- one_row(c(3, 1, 2, 2), names(clust))

  # Re = 2 p1 p2 (1 - d), and normalized Re = 1 - d
  re_a <- 2 * (3 / 4) * (1 / 4) * (1 - 0.2)
  re_b <- 2 * (1 / 2) * (1 / 2) * (1 - 0.7)
  expect_equal(unname(average.redundancy(abund, diss, clust)), mean(c(re_a, re_b)))
  expect_equal(unname(average.redundancy(abund, diss, clust, normalize = TRUE)), mean(c(0.8, 0.3)))
})


test_that("Only the distances within units are needed", {
  set.seed(5)
  study <- make_units(c(3, 4, 2), n_samples = 4)
  base <- average.redundancy(study$abund, study$diss_frame, study$clust)

  # Distances between units are ignored
  all_ids <- names(study$clust)
  full <- matrix(runif(length(all_ids)^2), length(all_ids), dimnames = list(all_ids, all_ids))
  full <- (full + t(full)) / 2
  for (u in seq_along(study$sub_ids)) {
    ids <- study$sub_ids[[u]]
    full[ids, ids] <- study$mats[[u]]
  }
  diag(full) <- 0
  expect_equal(average.redundancy(study$abund, mat_to_table(full, all_ids), study$clust), base)

  # Orientation of the pairs does not matter
  flipped <- study$diss_frame[, c(2, 1, 3)]
  expect_equal(average.redundancy(study$abund, flipped, study$clust), base)

  # A missing pair within a unit is an error
  expect_error(
    average.redundancy(study$abund, study$diss_frame[-1, ], study$clust),
    "Missing distances"
  )
})


test_that("Units with a single present subunit contribute 0", {
  clust <- c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B", c1 = "C")
  diss <- data.frame(
    ID1 = c("a1", "a1", "a2", "b1"),
    ID2 = c("a2", "a3", "a3", "b2"),
    Distance = c(0, 0, 0, 0)
  )

  # A is fully redundant (all distances 0), B and C have one present subunit
  abund <- one_row(c(1, 2, 3, 0, 49, 7), names(clust))
  re_a <- 1 - sum((1:3 / 6)^2)
  expect_equal(unname(average.redundancy(abund, diss, clust)), re_a / 3)
  expect_equal(unname(average.redundancy(abund, diss, clust, normalize = TRUE)), 1 / 3)

  # Only singletons: no redundancy
  abund <- one_row(c(0, 5, 0, 49, 0, 7), names(clust))
  expect_identical(unname(average.redundancy(abund, diss, clust)), 0)
  expect_identical(unname(average.redundancy(abund, diss, clust, normalize = TRUE)), 0)
})


test_that("Absent units are ignored by default", {
  set.seed(9)
  study <- make_units(c(4, 3, 5), n_samples = 1, zero_prob = 0)
  base <- average.redundancy(study$abund, study$diss_frame, study$clust)

  # Declaring extra absent units does not change the value
  clust <- c(study$clust, x1 = "X", x2 = "X", y1 = "Y")
  diss <- rbind(study$diss_frame, data.frame(ID1 = "x1", ID2 = "x2", Distance = 0.5))
  abund <- cbind(study$abund, x1 = 0, x2 = 0, y1 = 0)
  expect_equal(average.redundancy(abund, diss, clust), base)

  # With include_absent they count as zero
  expect_equal(average.redundancy(abund, diss, clust, include_absent = TRUE), base * 3 / 5)

  # Subunits named in clust but missing from abund have zero abundance
  expect_equal(average.redundancy(study$abund, diss, clust), base)
  expect_equal(average.redundancy(study$abund, diss, clust, include_absent = TRUE), base * 3 / 5)

  # An empty sample is NA, also with include_absent
  empty <- study$abund * 0
  expect_true(is.na(average.redundancy(empty, study$diss_frame, study$clust)))
  expect_true(is.na(average.redundancy(empty, study$diss_frame, study$clust, include_absent = TRUE)))
})


test_that("Several copies of a unit give the same value as one unit", {
  set.seed(3)
  n <- 6
  ab <- runif(n, 1, 100)
  d <- rand_diss(n)

  one <- make_units(n, 1)
  one_ab <- one_row(ab, one$sub_ids[[1]])
  one_diss <- mat_to_table(d, one$sub_ids[[1]])

  three <- make_units(c(n, n, n), 1)
  three_ab <- one_row(rep(ab, 3), names(three$clust))
  three_diss <- do.call(rbind, lapply(three$sub_ids, function(ids) mat_to_table(d, ids)))

  for (normalize in c(FALSE, TRUE)) {
    expect_equal(
      unname(average.redundancy(three_ab, three_diss, three$clust, normalize = normalize)),
      unname(average.redundancy(one_ab, one_diss, one$clust, normalize = normalize))
    )
  }
})


test_that("Invariant to the scale of the abundances", {
  set.seed(8)
  study <- make_units(c(3, 6, 4), n_samples = 5)
  for (normalize in c(FALSE, TRUE)) {
    expect_equal(
      average.redundancy(study$abund * 13.7, study$diss_frame, study$clust, normalize = normalize),
      average.redundancy(study$abund, study$diss_frame, study$clust, normalize = normalize)
    )
  }
})


test_that("Matrix, data.frame and sparse abundance formats give the same result", {
  set.seed(21)
  study <- make_units(c(4, 2, 5), n_samples = 6)
  base <- average.redundancy(study$abund, study$diss_frame, study$clust, normalize = TRUE)

  expect_equal(
    average.redundancy(as.data.frame(study$abund), study$diss_frame, study$clust, normalize = TRUE),
    base
  )
  expect_equal(
    average.redundancy(Matrix::Matrix(study$abund, sparse = TRUE), study$diss_frame, study$clust, normalize = TRUE),
    base
  )

  # Named clust in any order, and unnamed clust following the columns
  expect_equal(
    average.redundancy(study$abund, study$diss_frame, rev(study$clust), normalize = TRUE),
    base
  )
  expect_equal(
    average.redundancy(study$abund, study$diss_frame, unname(study$clust), normalize = TRUE),
    base
  )
})


test_that("Input validation", {
  clust <- c(a1 = "A", a2 = "A", b1 = "B")
  diss <- data.frame(ID1 = "a1", ID2 = "a2", Distance = 0.3)
  abund <- one_row(c(1, 2, 3), names(clust))

  expect_error(average.redundancy(abund, diss, clust, normalize = NA), "single TRUE or FALSE")
  expect_error(average.redundancy(abund, diss, clust, normalize = "yes"), "single TRUE or FALSE")
  expect_error(average.redundancy(abund, diss, clust, include_absent = 1), "single TRUE or FALSE")
  expect_error(average.redundancy(one_row(c(-1, 2, 3), names(clust)), diss, clust), "non-negative")
  expect_error(average.redundancy(abund, diss, clust[1:2]), "missing some subunits")
  expect_error(average.redundancy(abund, diss, unname(clust)[1:2]), "same number of subunits")
  expect_error(average.redundancy("a", diss, clust), "must be a matrix")
})
