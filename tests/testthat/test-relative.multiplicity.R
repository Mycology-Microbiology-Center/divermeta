# Tests for relative multiplicity

# Helpers
# ----------------------------

# Random study with units of the given sizes (samples in rows, subunits in columns)
make_study <- function(sizes, n_samples, zero_prob = 0.4) {
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

  list(
    clust = clust,
    unit_ids = unit_ids,
    sub_ids = sub_ids,
    mats = mats,
    diss_frame = diss_frame,
    abund = abund
  )
}

# Single sample given as lists of per-unit abundances and matrices, as the
# legacy implementation takes it
single_sample <- function(ab, diss) {
  sub_ids <- lapply(seq_along(ab), function(u) paste0("U", u, "_", seq_along(ab[[u]])))
  clust <- stats::setNames(rep(paste0("U", seq_along(ab)), lengths(ab)), unlist(sub_ids))
  list(
    abund = one_row(unlist(ab), names(clust)),
    diss = do.call(rbind, lapply(seq_along(ab), function(u) mat_to_table(as.matrix(diss[[u]]), sub_ids[[u]]))),
    clust = clust
  )
}

# Reference diversities of per-unit reference compositions (legacy implementation)
legacy_ref_div <- function(ab_ref, diss_ref, sigma = 1) {
  sigma <- rep(sigma, length.out = length(ab_ref))
  vapply(seq_along(ab_ref), function(i) {
    legacy$.rm_unit_diversity(ab_ref[[i]], diss_ref[[i]], sigma[i])
  }, numeric(1))
}

# Relative multiplicity of a single sample with arbitrary reference compositions
single_rm <- function(ab, diss, ab_ref, diss_ref, sigma = 1, ...) {
  x <- single_sample(ab, diss)
  unname(relative.multiplicity.ref_div(
    x$abund, x$diss, x$clust,
    ref_div = legacy_ref_div(ab_ref, diss_ref, sigma),
    sigma = sigma, ...
  ))
}

# Relative multiplicity of every sample computed with the legacy single sample function
manual_rm <- function(study, sigma = 1, cap_at_one = FALSE,
                      assume_max = FALSE, assume_homogeneous = FALSE) {
  sigma <- rep(sigma, length.out = length(study$unit_ids))
  pooled <- colSums(study$abund)

  ab_ref <- lapply(study$sub_ids, function(ids) {
    if (assume_homogeneous) rep(1, length(ids)) else pooled[ids]
  })
  diss_ref <- lapply(seq_along(study$sub_ids), function(u) {
    if (assume_max) max_diss(length(study$sub_ids[[u]]), sigma[u]) else study$mats[[u]]
  })

  sapply(rownames(study$abund), function(s) {
    ab <- lapply(study$sub_ids, function(ids) study$abund[s, ids])
    legacy$relative.multiplicity(ab, study$mats, ab_ref, diss_ref, sigma = sigma, cap_at_one = cap_at_one)
  })
}


# Single sample
# ----------------------------

# Maximally different subunits with equal abundances have diversity equal to
# their number, so the index reduces to the scenarios of the manuscript
test_that("Manuscript scenarios", {
  sizes <- c(4, 8, 16)
  diss <- lapply(sizes, max_diss)
  rm <- function(ab, ...) {
    x <- single_sample(ab, diss)
    unname(relative.multiplicity(
      x$abund, x$diss, x$clust,
      assume_max_reference_distance = TRUE,
      assume_homogeneous_abundance = TRUE,
      ...
    ))
  }

  # A: half of every unit
  ab <- lapply(sizes, function(n) rep(c(1, 0), n / 2))
  expect_equal(rm(ab), 0.5)

  # B: the same number of subunits in every unit
  ab <- lapply(sizes, function(n) c(rep(1, 4), rep(0, n - 4)))
  expect_equal(rm(ab), mean(c(1, 0.5, 0.25)))

  # C: only the largest unit, complete. Absent units are ignored by default
  ab <- list(rep(0, 4), rep(0, 8), rep(1, 16))
  expect_equal(rm(ab), 1)

  # With include_absent, absent units count as zero
  expect_equal(rm(ab, include_absent = TRUE), 1 / 3)

  # A and B have every unit present, so include_absent does not change them
  ab <- lapply(sizes, function(n) c(rep(1, 4), rep(0, n - 4)))
  expect_equal(rm(ab, include_absent = TRUE), rm(ab))

  # Same values as the legacy implementation
  ab_ref <- lapply(sizes, function(n) rep(1, n))
  ab <- list(rep(0, 4), rep(0, 8), rep(1, 16))
  expect_equal(rm(ab, include_absent = TRUE), legacy$relative.multiplicity(ab, diss, ab_ref, diss, include_absent = TRUE))
})


test_that("Absent units are ignored by default", {
  set.seed(42)
  sizes <- c(3, 5, 7)
  ab <- lapply(sizes, function(n) runif(n, 1, 100))
  diss <- lapply(sizes, rand_diss)
  ab_ref <- lapply(sizes, function(n) runif(n, 1, 100))
  diss_ref <- lapply(sizes, rand_diss)
  base <- single_rm(ab, diss, ab_ref, diss_ref, sigma = 0.6)
  expect_equal(base, legacy$relative.multiplicity(ab, diss, ab_ref, diss_ref, sigma = 0.6))

  # Declaring extra absent units does not change the value
  ab_ext <- c(ab, list(rep(0, 4), rep(0, 2)))
  diss_ext <- c(diss, list(rand_diss(4), rand_diss(2)))
  ab_ref_ext <- c(ab_ref, list(runif(4, 1, 100), runif(2, 1, 100)))
  diss_ref_ext <- c(diss_ref, list(rand_diss(4), rand_diss(2)))
  expect_equal(single_rm(ab_ext, diss_ext, ab_ref_ext, diss_ref_ext, sigma = 0.6), base)

  # With include_absent they count as zero
  expect_equal(
    single_rm(ab_ext, diss_ext, ab_ref_ext, diss_ref_ext, sigma = 0.6, include_absent = TRUE),
    base * 3 / 5
  )

  # No unit present: an empty sample is NA, also with include_absent
  ab_none <- lapply(sizes, function(n) rep(0, n))
  x <- single_sample(ab_none, diss)
  expect_true(is.na(relative.multiplicity.ref_div(x$abund, x$diss, x$clust, c(2, 3, 4))))
  expect_true(is.na(relative.multiplicity.ref_div(x$abund, x$diss, x$clust, c(2, 3, 4), include_absent = TRUE)))

  # Combined with cap_at_one
  n <- 5
  d <- max_diss(n)
  ab_cap <- list(rep(1, n), c(1, 1, 0, 0, 0), rep(0, n))
  ab_ref_cap <- list(c(1000, 1, 1, 1, 1), rep(1, n), rep(1, n))
  diss_cap <- list(d, d, d)
  expect_equal(
    single_rm(ab_cap, diss_cap, ab_ref_cap, diss_cap, cap_at_one = TRUE),
    mean(c(1, 2 / 5))
  )
  expect_equal(
    single_rm(ab_cap, diss_cap, ab_ref_cap, diss_cap, cap_at_one = TRUE, include_absent = TRUE),
    mean(c(1, 2 / 5, 0))
  )
})


test_that("Equals the average ratio of functional diversities", {
  set.seed(42)
  for (w_ in 1:10) {
    k <- 4
    sig <- runif(1, 0.3, 1)
    sizes <- sample(2:10, k, replace = TRUE)
    ref_sizes <- sizes + sample(0:5, k, replace = TRUE)

    ab <- lapply(sizes, function(n) runif(n, 1, 100))
    diss <- lapply(sizes, rand_diss)
    ab_ref <- lapply(ref_sizes, function(n) runif(n, 1, 100))
    diss_ref <- lapply(ref_sizes, rand_diss)

    expected <- mean(sapply(seq_len(k), function(i) {
      legacy$diversity.functional(ab[[i]], diss[[i]], sig) / legacy$diversity.functional(ab_ref[[i]], diss_ref[[i]], sig)
    }))

    expect_equal(single_rm(ab, diss, ab_ref, diss_ref, sigma = sig), expected)
    expect_equal(single_rm(ab, diss, ab_ref, diss_ref, sigma = sig), legacy$relative.multiplicity(ab, diss, ab_ref, diss_ref, sigma = sig))
  }
})


test_that("Several copies of a unit give the same value as one unit", {
  set.seed(42)
  n <- 8
  ab <- runif(n, 1, 100)
  diss <- rand_diss(n)
  ab_ref <- runif(n + 4, 1, 100)
  diss_ref <- rand_diss(n + 4)

  single <- single_rm(list(ab), list(diss), list(ab_ref), list(diss_ref), sigma = 0.7)

  for (m in c(2, 3, 7)) {
    copies <- single_rm(
      rep(list(ab), m), rep(list(diss), m),
      rep(list(ab_ref), m), rep(list(diss_ref), m),
      sigma = 0.7
    )
    expect_equal(copies, single)
  }
})


test_that("Sample equal to the reference gives one", {
  set.seed(42)
  sizes <- c(3, 6, 9)
  ab <- lapply(sizes, function(n) runif(n, 1, 100))
  diss <- lapply(sizes, rand_diss)

  expect_equal(single_rm(ab, diss, ab, diss, sigma = 0.5), 1)

  # The pooled reference of a single sample is the sample itself
  x <- single_sample(ab, diss)
  expect_equal(unname(relative.multiplicity(x$abund, x$diss, x$clust, sigma = 0.5)), 1)
})


test_that("cap_at_one caps units that are more diverse than their reference", {
  n <- 5
  diss <- max_diss(n)

  # Very uneven reference: much less diverse than an even sample
  ab_ref <- c(1000, 1, 1, 1, 1)
  ab <- rep(1, n)
  ratio <- n / legacy$diversity.functional(ab_ref, diss)

  expect_gt(single_rm(list(ab), list(diss), list(ab_ref), list(diss)), 1)
  expect_equal(single_rm(list(ab), list(diss), list(ab_ref), list(diss)), ratio)
  expect_equal(single_rm(list(ab), list(diss), list(ab_ref), list(diss), cap_at_one = TRUE), 1)

  # Mix of a capped unit and one below one
  ab_list <- list(ab, c(1, 1, 0, 0, 0))
  ab_ref_list <- list(ab_ref, rep(1, n))
  diss_list <- list(diss, diss)
  expect_equal(
    single_rm(ab_list, diss_list, ab_ref_list, diss_list, cap_at_one = TRUE),
    mean(c(1, 2 / 5))
  )
  expect_equal(
    single_rm(ab_list, diss_list, ab_ref_list, diss_list, cap_at_one = FALSE),
    mean(c(ratio, 2 / 5))
  )
})


test_that("sigma can be given per unit", {
  set.seed(42)
  study <- make_study(c(4, 6, 8), n_samples = 5)

  # Repeated value equals scalar
  expect_equal(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = rep(0.6, 3)),
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = 0.6)
  )

  # Different value per unit
  sig <- c(0.3, 0.6, 0.9)
  expect_equal(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = sig),
    manual_rm(study, sig)
  )

  # Distances above sigma are capped
  capped <- study$diss_frame
  capped$Distance <- pmin(capped$Distance, 0.3)
  expect_equal(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = 0.3),
    relative.multiplicity(study$abund, capped, study$clust, sigma = 0.3)
  )
})


test_that("Invariances", {
  set.seed(42)
  study <- make_study(c(5, 7, 9), n_samples = 6, zero_prob = 0.3)
  sig <- c(U1 = 0.4, U2 = 0.7, U3 = 1)
  rm <- function(abund, clust = study$clust, diss = study$diss_frame) {
    relative.multiplicity(abund, diss, clust, sigma = sig)
  }
  base <- rm(study$abund)

  # Permuting subunits (columns)
  perm <- sample(ncol(study$abund))
  expect_equal(rm(study$abund[, perm]), base)
  expect_equal(rm(study$abund[, perm], unname(study$clust[perm])), base)

  # Permuting units
  expect_equal(rm(study$abund, rev(study$clust)), base)

  # Permuting samples
  expect_equal(rm(study$abund[c(3, 1, 2, 6, 5, 4), ])[names(base)], base)

  # Scaling the abundances of every sample (the pooled reference scales too)
  expect_equal(rm(study$abund * 13), base)

  # Scaling the abundances of a single sample, with fixed reference diversities
  ref_div <- c(2, 3, 4)
  scaled <- study$abund
  scaled[2, ] <- scaled[2, ] * 13
  expect_equal(
    relative.multiplicity.ref_div(scaled, study$diss_frame, study$clust, ref_div, sigma = sig),
    relative.multiplicity.ref_div(study$abund, study$diss_frame, study$clust, ref_div, sigma = sig)
  )

  # Removing subunits with zero abundance in every sample
  study2 <- study
  study2$abund[, 1] <- 0
  keep <- colnames(study2$abund)[-1]
  expect_equal(
    relative.multiplicity(study2$abund, study$diss_frame, study$clust, sigma = sig),
    relative.multiplicity(study2$abund[, keep], study$diss_frame, study$clust[keep], sigma = sig)
  )
})


test_that("Input validation", {
  set.seed(1)
  study <- make_study(c(2, 3), n_samples = 2, zero_prob = 0)
  abund <- study$abund
  diss <- study$diss_frame
  clust <- study$clust

  expect_error(relative.multiplicity(abund, diss, clust, sigma = 0), "strictly positive")
  expect_error(relative.multiplicity(abund, diss, clust, sigma = -1), "strictly positive")
  expect_error(relative.multiplicity(abund, diss, clust, sigma = c(1, 1, 1)), "one value per unit")
  abund_na <- abund
  abund_na[1, 1] <- NA
  expect_error(relative.multiplicity(abund_na, diss, clust), "NA")
  abund_neg <- abund
  abund_neg[1, 1] <- -1
  expect_error(relative.multiplicity(abund_neg, diss, clust), "non-negative")
  expect_error(relative.multiplicity(abund, diss, clust, cap_at_one = NA), "cap_at_one")
  expect_error(relative.multiplicity(abund, diss, clust, include_absent = 1), "include_absent")
  expect_error(relative.multiplicity.ref_div(abund, diss, clust, c(1, 1), include_absent = NA), "include_absent")
  expect_error(relative.multiplicity(abund, diss, unname(clust)[1:4]), "same number of subunits")
  expect_error(relative.multiplicity(abund, diss[, 1:2], clust), "three columns")
  # Every sample empty: all NA
  expect_equal(
    relative.multiplicity(matrix(0, 2, 5, dimnames = list(c("S1", "S2"), names(clust))), diss, clust),
    c(S1 = NA_real_, S2 = NA_real_)
  )
})


# From reference diversities
# ----------------------------

test_that("ref_div version equals the full version", {
  set.seed(42)
  for (w_ in 1:5) {
    sizes <- sample(2:10, 4, replace = TRUE)
    ab <- lapply(sizes, function(n) runif(n, 1, 100) * (runif(n) > 0.3))
    diss <- lapply(sizes, rand_diss)
    ab_ref <- lapply(sizes, function(n) runif(n, 1, 100))
    diss_ref <- lapply(sizes, rand_diss)
    sig <- runif(4, 0.3, 1)

    for (cap in c(FALSE, TRUE)) {
      expect_equal(
        single_rm(ab, diss, ab_ref, diss_ref, sigma = sig, cap_at_one = cap),
        legacy$relative.multiplicity(ab, diss, ab_ref, diss_ref, sigma = sig, cap_at_one = cap)
      )
    }
  }
})


test_that("ref_div equal to the number of subunits gives the manuscript scenarios", {
  sizes <- c(4, 8, 16)
  diss <- lapply(sizes, max_diss)
  rm <- function(ab, ...) {
    x <- single_sample(ab, diss)
    unname(relative.multiplicity.ref_div(x$abund, x$diss, x$clust, sizes, ...))
  }

  ab <- lapply(sizes, function(n) rep(c(1, 0), n / 2))
  expect_equal(rm(ab), 0.5)

  ab <- lapply(sizes, function(n) c(rep(1, 4), rep(0, n - 4)))
  expect_equal(rm(ab), mean(c(1, 0.5, 0.25)))

  ab <- list(rep(0, 4), rep(0, 8), rep(1, 16))
  expect_equal(rm(ab), 1)
  expect_equal(rm(ab, include_absent = TRUE), 1 / 3)
})


test_that("ref_div input validation", {
  x <- single_sample(list(c(1, 2), c(3, 4, 5)), list(max_diss(2), max_diss(3)))
  rm <- function(ref_div) relative.multiplicity.ref_div(x$abund, x$diss, x$clust, ref_div)

  expect_error(rm(c(1, 0)), "strictly positive")
  expect_error(rm(c(1, -2)), "strictly positive")
  expect_error(rm(c(1, NA)), "strictly positive")
  expect_error(rm(c(1, Inf)), "strictly positive")
  expect_error(rm(c(1, 2, 3)), "one value per unit")
  expect_error(rm(2), "one value per unit")
})


# Many samples
# ----------------------------

test_that("Matrix, data.frame and sparse abundance formats give the same result", {
  set.seed(42)
  study <- make_study(c(4, 6, 3, 8), n_samples = 6)

  for (assume_max in c(FALSE, TRUE)) {
    for (assume_homogeneous in c(FALSE, TRUE)) {
      rm <- function(abund) {
        relative.multiplicity(
          abund, study$diss_frame, study$clust,
          sigma = 0.6,
          assume_max_reference_distance = assume_max,
          assume_homogeneous_abundance = assume_homogeneous
        )
      }
      wide_res <- rm(study$abund)
      expect_equal(rm(as.data.frame(study$abund)), wide_res)
      expect_equal(rm(Matrix::Matrix(study$abund, sparse = TRUE)), wide_res)
    }
  }

  # Numeric identifiers
  num_clust <- stats::setNames(study$clust, seq_along(study$clust))
  num_frame <- study$diss_frame
  num_frame$ID1 <- match(num_frame$ID1, names(study$clust))
  num_frame$ID2 <- match(num_frame$ID2, names(study$clust))
  num_abund <- study$abund
  colnames(num_abund) <- seq_along(study$clust)
  expect_equal(
    relative.multiplicity(num_abund, num_frame, num_clust),
    relative.multiplicity(study$abund, study$diss_frame, study$clust)
  )
})


test_that("Many samples match the legacy single sample function", {
  set.seed(42)
  for (w_ in 1:3) {
    study <- make_study(sample(2:8, 4, replace = TRUE), n_samples = 5)
    sig <- runif(4, 0.3, 1)

    for (assume_max in c(FALSE, TRUE)) {
      for (assume_homogeneous in c(FALSE, TRUE)) {
        for (cap in c(FALSE, TRUE)) {
          res <- relative.multiplicity(
            study$abund, study$diss_frame, study$clust,
            sigma = sig, cap_at_one = cap,
            assume_max_reference_distance = assume_max,
            assume_homogeneous_abundance = assume_homogeneous
          )
          expected <- manual_rm(study, sig, cap, assume_max, assume_homogeneous)
          expect_equal(res, expected)
        }
      }
    }
  }
})


test_that("Several copies of a unit give the same value in every sample", {
  set.seed(42)
  study <- make_study(5, n_samples = 4)

  # Duplicate the unit with new subunit ids and a new unit label
  copy_ids <- paste0("copy_", names(study$clust))
  clust2 <- c(study$clust, stats::setNames(rep("U_copy", length(copy_ids)), copy_ids))
  frame_copy <- study$diss_frame
  frame_copy$ID1 <- paste0("copy_", frame_copy$ID1)
  frame_copy$ID2 <- paste0("copy_", frame_copy$ID2)
  diss_frame2 <- rbind(study$diss_frame, frame_copy)
  abund2 <- cbind(study$abund, study$abund)
  colnames(abund2) <- names(clust2)

  for (assume_max in c(FALSE, TRUE)) {
    for (assume_homogeneous in c(FALSE, TRUE)) {
      expect_equal(
        relative.multiplicity(
          abund2, diss_frame2, clust2,
          assume_max_reference_distance = assume_max,
          assume_homogeneous_abundance = assume_homogeneous
        ),
        relative.multiplicity(
          study$abund, study$diss_frame, study$clust,
          assume_max_reference_distance = assume_max,
          assume_homogeneous_abundance = assume_homogeneous
        )
      )
    }
  }
})


test_that("A single sample with pooled reference gives one", {
  set.seed(42)
  study <- make_study(c(4, 6, 3), n_samples = 1, zero_prob = 0.2)

  expect_equal(
    unname(relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = 0.5)),
    1
  )
})


test_that("Homogeneous and maximum distance reference is the number of subunits", {
  set.seed(42)
  study <- make_study(c(4, 8, 16), n_samples = 5)

  res <- relative.multiplicity(
    study$abund, study$diss_frame, study$clust,
    assume_max_reference_distance = TRUE,
    assume_homogeneous_abundance = TRUE
  )
  expect_true(all(res <= 1))

  expect_equal(
    res,
    relative.multiplicity.ref_div(study$abund, study$diss_frame, study$clust, ref_div = c(4, 8, 16))
  )

  # Manuscript scenario A: half of every unit, maximally different, equal abundances
  sizes <- c(4, 8, 16)
  ids <- lapply(seq_along(sizes), function(u) paste0("U", u, "_", seq_len(sizes[u])))
  clust <- stats::setNames(rep(paste0("U", 1:3), sizes), unlist(ids))
  diss_frame <- do.call(rbind, lapply(seq_along(sizes), function(u) mat_to_table(max_diss(sizes[u]), ids[[u]])))
  abund <- matrix(
    unlist(lapply(sizes, function(n) rep(c(1, 0), n / 2))),
    nrow = 1,
    dimnames = list("A", names(clust))
  )
  expect_equal(
    unname(relative.multiplicity(
      abund, diss_frame, clust,
      assume_max_reference_distance = TRUE,
      assume_homogeneous_abundance = TRUE
    )),
    0.5
  )

  # Subunits listed in clust but not in abund count in the homogeneous reference
  half <- abund[, abund[1, ] > 0, drop = FALSE]
  expect_equal(
    unname(relative.multiplicity(
      half, diss_frame, clust,
      assume_max_reference_distance = TRUE,
      assume_homogeneous_abundance = TRUE
    )),
    0.5
  )
})


test_that("Distance table is robust to orientation and extra rows", {
  set.seed(42)
  study <- make_study(c(4, 6, 3), n_samples = 4)
  base <- relative.multiplicity(study$abund, study$diss_frame, study$clust)

  # Pairs in both orientations
  flipped <- study$diss_frame
  flipped[, c("ID1", "ID2")] <- flipped[, c("ID2", "ID1")]
  expect_equal(relative.multiplicity(study$abund, rbind(study$diss_frame, flipped), study$clust), base)

  # Random orientation per row, and other column names
  swap <- runif(nrow(study$diss_frame)) > 0.5
  mixed <- study$diss_frame
  mixed[swap, c("ID1", "ID2")] <- mixed[swap, c("ID2", "ID1")]
  colnames(mixed) <- c("a", "b", "d")
  expect_equal(relative.multiplicity(study$abund, mixed, study$clust), base)

  # Cross unit rows, self pairs and unknown subunits are ignored
  extra <- data.frame(
    ID1 = c(study$sub_ids[[1]][1], study$sub_ids[[2]][1], "unknown"),
    ID2 = c(study$sub_ids[[2]][1], study$sub_ids[[2]][1], study$sub_ids[[3]][1]),
    Distance = c(0.01, 0.5, 0.01)
  )
  expect_equal(relative.multiplicity(study$abund, rbind(study$diss_frame, extra), study$clust), base)

  # Missing within unit distance
  expect_error(
    relative.multiplicity(study$abund, study$diss_frame[-1, ], study$clust),
    "Missing distances.*unit: U1"
  )

  # Conflicting distances for the same pair
  conflict <- study$diss_frame[1, ]
  conflict$Distance <- conflict$Distance + 0.1
  expect_error(
    relative.multiplicity(study$abund, rbind(study$diss_frame, conflict), study$clust),
    "conflicting"
  )
})


test_that("Many samples edge cases", {
  set.seed(42)
  study <- make_study(c(4, 6, 3), n_samples = 4)

  # Empty sample is NA
  abund <- study$abund
  abund["S2", ] <- 0
  res <- relative.multiplicity(abund, study$diss_frame, study$clust)
  expect_true(is.na(res[["S2"]]))
  expect_false(any(is.na(res[c("S1", "S3", "S4")])))

  # A unit absent from a sample is ignored by default, and counts as zero with include_absent
  abund <- study$abund
  abund["S1", study$sub_ids[[2]]] <- 0
  ref_div <- c(U1 = 2, U2 = 3, U3 = 1.5)
  ab <- lapply(study$sub_ids, function(ids) abund["S1", ids])
  r1 <- legacy$diversity.functional(ab[[1]][ab[[1]] > 0], study$mats[[1]][ab[[1]] > 0, ab[[1]] > 0]) / 2
  r3 <- legacy$diversity.functional(ab[[3]][ab[[3]] > 0], study$mats[[3]][ab[[3]] > 0, ab[[3]] > 0]) / 1.5

  res <- relative.multiplicity.ref_div(abund, study$diss_frame, study$clust, ref_div)
  expect_equal(res[["S1"]], legacy$relative.multiplicity.ref_div(ab, study$mats, ref_div))
  expect_equal(res[["S1"]], mean(c(r1, r3)))

  res <- relative.multiplicity.ref_div(
    abund, study$diss_frame, study$clust, ref_div,
    include_absent = TRUE
  )
  expect_equal(res[["S1"]], legacy$relative.multiplicity.ref_div(ab, study$mats, ref_div, include_absent = TRUE))
  expect_equal(res[["S1"]], mean(c(r1, 0, r3)))

  # Subunit not in clust
  abund <- cbind(study$abund, extra = rep(1, nrow(study$abund)))
  expect_error(relative.multiplicity(abund, study$diss_frame, study$clust), "missing some subunits")
  abund <- cbind(study$abund, extra = rep(0, nrow(study$abund)))
  expect_error(relative.multiplicity(abund, study$diss_frame, study$clust), "missing some subunits")

  # Named sigma in a different order than the units
  sig <- c(U1 = 0.3, U2 = 0.6, U3 = 0.9)
  expect_equal(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = sig[c(3, 1, 2)]),
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = unname(sig))
  )
  expect_error(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, sigma = c(U1 = 1, U2 = 1)),
    "missing values for units"
  )

  # Unit with no abundance in the study: ignored silently by default
  abund <- study$abund
  abund[, study$sub_ids[[3]]] <- 0
  keep <- names(study$clust) %in% c(study$sub_ids[[1]], study$sub_ids[[2]])
  expect_no_warning(
    res <- relative.multiplicity(abund, study$diss_frame, study$clust)
  )
  expect_equal(
    res,
    relative.multiplicity(abund[, keep], study$diss_frame, study$clust[keep])
  )

  # With include_absent and the pooled reference, it is dropped with a warning
  expect_warning(
    res <- relative.multiplicity(abund, study$diss_frame, study$clust, include_absent = TRUE),
    "dropped: U3"
  )
  expect_equal(
    res,
    relative.multiplicity(abund[, keep], study$diss_frame, study$clust[keep], include_absent = TRUE)
  )

  # With the homogeneous reference it is ignored by default, and counts as zero with include_absent
  res_h_kept <- relative.multiplicity(
    abund[, keep], study$diss_frame, study$clust[keep],
    assume_homogeneous_abundance = TRUE
  )
  res_h <- relative.multiplicity(abund, study$diss_frame, study$clust, assume_homogeneous_abundance = TRUE)
  expect_equal(res_h, res_h_kept)
  res_h <- relative.multiplicity(
    abund, study$diss_frame, study$clust,
    assume_homogeneous_abundance = TRUE, include_absent = TRUE
  )
  expect_equal(res_h, res_h_kept * 2 / 3)

  # Flags must be logical
  expect_error(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, assume_max_reference_distance = "yes"),
    "assume_max_reference_distance"
  )
  expect_error(
    relative.multiplicity(study$abund, study$diss_frame, study$clust, include_absent = "yes"),
    "include_absent"
  )
  expect_error(
    relative.multiplicity.ref_div(study$abund, study$diss_frame, study$clust, c(1, 1, 1), include_absent = NA),
    "include_absent"
  )
})


# Many samples from reference diversities
# ----------------------------

test_that("Many samples ref_div version equals the full version", {
  set.seed(42)
  study <- make_study(c(4, 6, 3, 5), n_samples = 5)
  sig <- c(0.4, 0.6, 0.8, 1)
  pooled <- colSums(study$abund)

  # Pooled abundances and full matrices
  ref_div <- sapply(1:4, function(u) {
    legacy$diversity.functional(pooled[study$sub_ids[[u]]], study$mats[[u]], sig[u])
  })
  for (cap in c(FALSE, TRUE)) {
    expect_equal(
      relative.multiplicity.ref_div(
        study$abund, study$diss_frame, study$clust, ref_div,
        sigma = sig, cap_at_one = cap
      ),
      relative.multiplicity(
        study$abund, study$diss_frame, study$clust,
        sigma = sig, cap_at_one = cap
      )
    )
  }

  # Homogeneous abundances and maximum distances
  expect_equal(
    relative.multiplicity.ref_div(
      study$abund, study$diss_frame, study$clust, c(4, 6, 3, 5),
      sigma = sig
    ),
    relative.multiplicity(
      study$abund, study$diss_frame, study$clust,
      sigma = sig,
      assume_max_reference_distance = TRUE,
      assume_homogeneous_abundance = TRUE
    )
  )
})


test_that("Named ref_div is matched by unit", {
  set.seed(42)
  study <- make_study(c(4, 6, 3), n_samples = 4)
  ref_div <- c(U1 = 2, U2 = 3.5, U3 = 1.2)
  rm <- function(ref_div) relative.multiplicity.ref_div(study$abund, study$diss_frame, study$clust, ref_div)
  base <- rm(unname(ref_div))

  expect_equal(rm(ref_div[c(2, 3, 1)]), base)
  expect_equal(rm(c(ref_div, U9 = 10)), base)
  expect_error(rm(ref_div[1:2]), "missing values for units: U3")
  expect_error(rm(c(1, 2)), "one value per unit")
  expect_error(rm(c(1, 0, 2)), "strictly positive")
})


test_that("Unnamed sigma and ref_div follow unique(clust), not the columns of abund", {
  # unique(clust) is B, A, while the columns of abund start with A
  clust <- c(b1 = "B", b2 = "B", a1 = "A", a2 = "A", a3 = "A")
  abund <- matrix(
    c(1, 2, 3, 4, 5,
      2, 2, 1, 1, 1),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("a1", "a2", "a3", "b1", "b2"))
  )
  diss <- data.frame(
    ID1 = c("a1", "a1", "a2", "b1"),
    ID2 = c("a2", "a3", "a3", "b2"),
    Distance = c(0.3, 0.8, 0.6, 0.4)
  )

  expect_equal(
    relative.multiplicity.ref_div(abund, diss, clust, ref_div = c(10, 1)),
    relative.multiplicity.ref_div(abund, diss, clust, ref_div = c(B = 10, A = 1))
  )
  expect_equal(
    relative.multiplicity(abund, diss, clust, sigma = c(0.2, 1)),
    relative.multiplicity(abund, diss, clust, sigma = c(B = 0.2, A = 1))
  )
  # The named values differ between units, so the order matters
  expect_false(isTRUE(all.equal(
    relative.multiplicity(abund, diss, clust, sigma = c(B = 0.2, A = 1)),
    relative.multiplicity(abund, diss, clust, sigma = c(A = 0.2, B = 1))
  )))
})


test_that("Unnamed sigma follows a named clust given in reversed order", {
  set.seed(42)
  study <- make_study(c(4, 6, 3), n_samples = 5)
  clust <- rev(study$clust)
  sigma <- c(U1 = 0.3, U2 = 0.6, U3 = 0.9)
  sigma_unnamed <- unname(sigma[unique(clust)])

  expect_equal(
    relative.multiplicity(study$abund, study$diss_frame, clust, sigma = sigma_unnamed),
    relative.multiplicity(study$abund, study$diss_frame, clust, sigma = sigma)
  )
  expect_equal(
    relative.multiplicity.ref_div(study$abund, study$diss_frame, clust, c(2, 3, 4), sigma = sigma_unnamed),
    relative.multiplicity.ref_div(
      study$abund, study$diss_frame, clust, stats::setNames(c(2, 3, 4), unique(clust)),
      sigma = sigma
    )
  )
})


test_that("Unnamed sigma follows unique(clust) with subunits missing from abund", {
  set.seed(42)
  study <- make_study(c(4, 6, 3), n_samples = 5)
  # Unit U4 only has subunits missing from abund, and is listed first
  clust <- c(c(U4_1 = "U4", U4_2 = "U4"), study$clust)
  diss <- rbind(study$diss_frame, data.frame(ID1 = "U4_1", ID2 = "U4_2", Distance = 0.5))
  sigma <- c(U4 = 0.2, U1 = 0.3, U2 = 0.6, U3 = 0.9)

  for (include_absent in c(FALSE, TRUE)) {
    expect_equal(
      relative.multiplicity(
        study$abund, diss, clust, sigma = unname(sigma),
        include_absent = include_absent, assume_homogeneous_abundance = TRUE
      ),
      relative.multiplicity(
        study$abund, diss, clust, sigma = rev(sigma),
        include_absent = include_absent, assume_homogeneous_abundance = TRUE
      )
    )
  }
})


test_that("Named sigma and ref_div match numeric unit identifiers as numbers", {
  set.seed(7)
  study <- make_study(c(4, 6, 3), n_samples = 4)
  # Units of 1e5 or more: names() writes them "1e+05", the clustering "100000"
  units <- c(U1 = 100000, U2 = 200000, U3 = 300000)
  clust <- stats::setNames(unname(units[study$clust]), names(study$clust))

  sigma <- c(0.8, 1, 0.6)
  ref_div <- c(2, 3.5, 1.2)
  named <- function(x, nms) stats::setNames(x, nms)
  expect_identical(names(named(sigma, units)), c("1e+05", "2e+05", "3e+05"))

  rm <- function(clust, sigma) relative.multiplicity(study$abund, study$diss_frame, clust, sigma)
  rm_ref <- function(clust, ref_div) {
    relative.multiplicity.ref_div(study$abund, study$diss_frame, clust, ref_div, sigma = 0.8)
  }
  base <- rm(study$clust, named(sigma, study$unit_ids))
  base_ref <- rm_ref(study$clust, named(ref_div, study$unit_ids))

  expect_equal(rm(clust, named(sigma, units)), base)
  expect_equal(rm(clust, named(sigma, units)[c(3, 1, 2)]), base)
  expect_equal(rm(clust, named(sigma, c("100000", "200000", "300000"))), base)
  expect_equal(rm_ref(clust, named(ref_div, units)[c(2, 3, 1)]), base_ref)
  expect_equal(rm_ref(clust, c(named(ref_div, units), U9 = 10)), base_ref)

  expect_error(rm(clust, named(sigma, units)[1:2]), "missing values for units: 300000")
  expect_error(rm_ref(clust, named(ref_div, c(units[1:2], 400000))), "missing values for units: 300000")
})
