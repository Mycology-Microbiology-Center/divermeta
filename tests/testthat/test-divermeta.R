## Tests for the `divermeta` function


## Function to create a test data
make_fixture <- function() {
  species <- c("f1", "f2", "f3", "f4")
  samples <- c("S1", "S2", "S_empty")
  abund <- matrix(
    c(
      10, 5, 0, 6, # S1
      0,  8, 12, 3, # S2
      0,  0, 0, 0 # S_empty
    ),
    nrow = length(samples),
    ncol = length(species),
    byrow = TRUE,
    dimnames = list(samples, species)
  )

  diss_mtx <- matrix(
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

  clust <- c(f1 = "A", f2 = "A", f3 = "B", f4 = "B")

  list(abund = abund, diss = mat_to_table(diss_mtx), diss_mtx = diss_mtx, clust = clust)
}

all_indices <- c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq")
normalized <- c("multiplicity_inventory", "raoQ", "FD_sigma", "redundancy", "multiplicity_distance", "FD_q")

# With the indices that average over the units, which the legacy divermeta lacks
unit_indices <- c("RM", "AR")
every_index <- c(all_indices, unit_indices)
every_normalized <- c(normalized, "relative_multiplicity", "average_redundancy")


test_that("single index multiplicity_inventory works and matches direct", {
  fx <- make_fixture()

  res <- divermeta(fx$abund,
    indices = c("multiplicity_inventory"),
    clust = fx$clust
  )

  expect_true(is.data.frame(res))
  expect_identical(colnames(res), c("Sample", "multiplicity_inventory"))
  expect_identical(res$Sample, rownames(fx$abund))

  expected <- apply(fx$abund, 1, function(row) {
    if (sum(row) <= 0) {
      return(NA_real_)
    }
    legacy$multiplicity.inventory(row, fx$clust, q = 1)
  })
  expect_equal(res$multiplicity_inventory, as.numeric(expected))
  expect_equal(res$multiplicity_inventory, unname(multiplicity.inventory(fx$abund, fx$clust, q = 1)))
})




test_that("multiple indices and alias mapping; matches direct functions", {
  fx <- make_fixture()
  sig <- 0.8
  q <- 1

  res <- divermeta(fx$abund,
    diss = fx$diss,
    indices = every_index,
    clust = fx$clust,
    q = q,
    sig = sig
  )

  # Column names are normalized
  expect_identical(colnames(res), c("Sample", every_normalized))

  expect_equal(res$multiplicity_inventory, unname(multiplicity.inventory(fx$abund, fx$clust, q = q)))
  expect_equal(res$raoQ, unname(raoQuadratic(fx$abund, fx$diss)))
  expect_equal(res$FD_sigma, unname(diversity.functional(fx$abund, fx$diss, sig = sig)))
  expect_equal(res$redundancy, unname(redundancy(fx$abund, fx$diss)))
  expect_equal(res$multiplicity_distance, unname(multiplicity.distance(fx$abund, fx$diss, fx$clust, method = "sigma", sig = sig)))
  expect_equal(res$FD_q, unname(diversity.functional.traditional(fx$abund, fx$diss, q = q)))
  expect_equal(res$relative_multiplicity, unname(relative.multiplicity(fx$abund, fx$diss, fx$clust, sigma = sig)))
  expect_equal(res$average_redundancy, unname(average.redundancy(fx$abund, fx$diss, fx$clust)))

  res_norm <- divermeta(fx$abund,
    diss = fx$diss,
    indices = every_index,
    clust = fx$clust,
    q = q,
    sig = sig,
    normalize = TRUE
  )

  # Checks normalization
  for (ind in every_normalized) {
    expect_equal(max(res_norm[[ind]], na.rm = TRUE), 1.0)
  }
})


test_that("matches the legacy divermeta", {
  fx <- make_fixture()
  q <- 2

  # Every distance between subunits of different clusters (0.6 or more) is at
  # least sigma, so the legacy and new sigma methods agree
  sig <- 0.5
  res <- divermeta(fx$abund, diss = fx$diss, indices = all_indices, clust = fx$clust, q = q, sig = sig)
  old <- legacy$divermeta(t(fx$abund), diss = fx$diss_mtx, indices = all_indices, clusters = fx$clust, q = q, sig = sig)
  expect_equal(res, old)

  # Linkage methods (legacy is given a dist object: it is wrong for full matrices)
  for (method in c("min", "max", "average")) {
    res <- divermeta(fx$abund, diss = fx$diss, indices = "M_distance", clust = fx$clust, sig = 0.8, method = method)
    old <- legacy$divermeta(t(fx$abund), diss = stats::as.dist(fx$diss_mtx), indices = "M_distance", clusters = fx$clust, sig = 0.8, method = method)
    expect_equal(res, old)
  }

  # Custom cluster distances
  diss_clust <- data.frame(ID1 = "A", ID2 = "B", Distance = 0.55)
  res <- divermeta(fx$abund, diss = fx$diss, indices = "M_distance", clust = fx$clust, sig = 0.8, method = "custom", diss_clust = diss_clust)
  expected <- apply(fx$abund[1:2, ], 1, function(row) {
    legacy$multiplicity.distance(row, fx$diss_mtx, fx$clust, method = "custom", sig = 0.8, clust_ids_order = c("A", "B"), diss_clust = matrix(c(0, 0.55, 0.55, 0), 2))
  })
  expect_equal(res$multiplicity_distance, c(unname(expected), NA))
})



test_that("matrix, data.frame and sparse abundances give the same result", {
  fx <- make_fixture()
  sig <- 0.8
  q <- 1

  res_m <- divermeta(fx$abund, diss = fx$diss, indices = every_index, clust = fx$clust, q = q, sig = sig)
  res_df <- divermeta(as.data.frame(fx$abund), diss = fx$diss, indices = every_index, clust = fx$clust, q = q, sig = sig)
  res_sp <- divermeta(Matrix::Matrix(fx$abund, sparse = TRUE), diss = fx$diss, indices = every_index, clust = fx$clust, q = q, sig = sig)

  expect_equal(res_df, res_m)
  expect_equal(res_sp, res_m)
})



test_that("errors when diss required but missing; and when clust required but missing", {
  fx <- make_fixture()

  # diss-required indices
  expect_error(divermeta(fx$abund, indices = c("raoQ")), "diss")

  # clust-required indices
  expect_error(divermeta(fx$abund, indices = c("multiplicity_inventory")), "clust")

  for (ind in c(unit_indices, "relative_multiplicity", "average_redundancy")) {
    expect_error(divermeta(fx$abund, indices = ind, clust = fx$clust), "diss")
    expect_error(divermeta(fx$abund, diss = fx$diss, indices = ind), "clust")
  }
})


test_that("unsupported index label errors clearly", {
  fx <- make_fixture()
  expect_error(divermeta(fx$abund, indices = c("not_an_index")), "Unsupported")
})


test_that("non-numeric abund errors; and incomplete diss errors", {
  fx <- make_fixture()

  bad_abund <- matrix(as.character(1:12), nrow = 3, dimnames = dimnames(fx$abund))
  expect_error(divermeta(bad_abund, indices = c("multiplicity_inventory"), clust = fx$clust), "numeric")

  # Subunit of abund missing from diss
  expect_error(divermeta(fx$abund, diss = fx$diss[fx$diss$ID1 != "f1", ], indices = c("raoQ")), "Missing distances")
})


test_that("extra subunits and order of diss do not matter", {
  fx <- make_fixture()

  # Add an extra subunit to diss and shuffle its rows and orientation
  extra <- data.frame(ID1 = "extra", ID2 = c("f1", "f2", "f3", "f4"), Distance = 0.3)
  diss2 <- rbind(fx$diss, extra)
  diss2 <- diss2[c(3, 8, 1, 6, 2, 10, 4, 7, 5, 9), ]
  diss2[c(1, 4), c("ID1", "ID2")] <- diss2[c(1, 4), c("ID2", "ID1")]

  res <- divermeta(fx$abund, diss = diss2, indices = c("raoQ"))
  expect_equal(res, divermeta(fx$abund, diss = fx$diss, indices = c("raoQ")))
  expect_true(all(is.finite(res$raoQ[1:2])))
  expect_true(is.na(res$raoQ[3])) # empty sample should be NA
})


test_that("named clust align regardless of order", {
  fx <- make_fixture()
  cl1 <- fx$clust
  cl2 <- fx$clust[c("f4", "f3", "f2", "f1")] # different order, still named

  r1 <- divermeta(fx$abund, indices = c("multiplicity_inventory"), clust = cl1)
  r2 <- divermeta(fx$abund, indices = c("multiplicity_inventory"), clust = cl2)

  expect_equal(r1$multiplicity_inventory, r2$multiplicity_inventory)
})


test_that("all-zero samples produce NA for applicable indices", {
  fx <- make_fixture()

  res <- divermeta(
    fx$abund, diss = fx$diss, indices = every_index, clust = fx$clust, sig = 0.8, include_absent = TRUE
  )

  # Third sample is all zeros
  for (ind in every_normalized) {
    expect_true(is.na(res[[ind]][3]))
  }
})


test_that("divermeta requires a single non-negative finite q", {
  fx <- make_fixture()
  for (q in list(NA_real_, Inf)) {
    for (ind in c("FD_q", "M_inventory")) {
      expect_error(
        divermeta(fx$abund, diss = fx$diss, indices = ind, clust = fx$clust, q = q),
        "`q` must be a single non-negative numeric value"
      )
    }
  }
})


test_that("repeated indices and aliases give repeated columns", {
  fx <- make_fixture()
  res <- divermeta(
    fx$abund, fx$diss,
    c("raoQ", "raoQ", "M_inventory", "multiplicity_inventory"),
    fx$clust
  )
  expect_identical(
    colnames(res),
    c("Sample", "raoQ", "raoQ.1", "multiplicity_inventory", "multiplicity_inventory.1")
  )
  expect_equal(res$raoQ.1, res$raoQ)
  expect_equal(res$multiplicity_inventory.1, res$multiplicity_inventory)
})


test_that("normalized columns are between 0 and 1, also with distances above 1", {
  fx <- make_fixture()
  far <- transform(fx$diss, Distance = Distance * 3)
  res <- divermeta(fx$abund, far, every_index, fx$clust, sig = 0.8, normalize = TRUE)
  for (ind in every_normalized) {
    vals <- res[[ind]][!is.na(res[[ind]])]
    expect_true(all(vals >= 0 & vals <= 1), info = ind)
    expect_equal(max(vals), 1, info = ind)
  }
})


test_that("the arguments of distance-based multiplicity are validated", {
  fx <- make_fixture()
  expect_error(
    divermeta(fx$abund, fx$diss, "M_distance", fx$clust, method = "custom"),
    "diss_clust cannot be NULL"
  )
  expect_error(divermeta(fx$abund, fx$diss, "M_distance", fx$clust, method = "median"), "not supported")
  expect_error(divermeta(fx$abund, fx$diss, "M_distance", fx$clust, sig = 0), "positive")
  expect_error(
    divermeta(fx$abund, fx$diss, "M_distance", fx$clust, method = "custom", diss_clust = data.frame("A", "C", 0.5)),
    "Missing distances in `diss_clust`"
  )
  # Not used by the other indices
  expect_no_error(divermeta(fx$abund, fx$diss, "raoQ", method = "median", sig = 0))
})


test_that("results do not depend on how the pairs are split into internal chunks", {
  set.seed(1)
  n <- 30
  ids <- paste0("f", seq_len(n))
  clust <- stats::setNames(sample(c("A", "B", "C", "D"), n, replace = TRUE), ids)
  tab <- mat_to_table(rand_diss(n, 0.05, 1.2), ids)
  abund <- matrix(
    runif(7 * n, 1, 50) * (runif(7 * n) > 0.3),
    nrow = 7,
    dimnames = list(paste0("S", 1:7), ids)
  )
  run <- function() {
    list(
      divermeta(abund, tab, c("raoQ", "FD_sigma", "FDq", "redundancy", "M_distance", "RM", "AR"), clust, q = 1, sig = 0.7),
      divermeta(abund, tab, c("FDq", "M_distance", "RM", "AR"), clust, q = 2, sig = 0.7, method = "average"),
      divermeta(abund, tab, c("RM", "AR"), clust, sig = 0.6, assume_homogeneous_abundance = TRUE),
      relative.multiplicity(abund, tab, clust, sigma = 0.6),
      average.redundancy(abund, tab, clust, normalize = TRUE)
    )
  }
  expected <- run()

  # A few cells per internal chunk: many chunks of pairs for every read
  local_mocked_bindings(.chunk_cells = 13)
  expect_equal(run(), expected)
})


test_that("normalize must be a single TRUE or FALSE", {
  fx <- make_fixture()
  for (normalize in list("yes", NA, NULL, c(TRUE, TRUE), 1)) {
    expect_error(
      divermeta(fx$abund, fx$diss, "raoQ", normalize = normalize),
      "`normalize` must be a single TRUE or FALSE value"
    )
  }
})


test_that("normalizing a study of empty samples gives NA without a warning", {
  fx <- make_fixture()
  empty <- fx$abund * 0
  expect_no_warning(
    res <- divermeta(empty, fx$diss, every_index, fx$clust, sig = 0.8, normalize = TRUE)
  )
  for (ind in every_normalized) {
    expect_true(all(is.na(res[[ind]])), info = ind)
  }
})


test_that("the Sample column holds the row names of abund, or the row numbers", {
  fx <- make_fixture()
  res <- divermeta(fx$abund, fx$diss, "raoQ")
  expect_identical(res$Sample, rownames(fx$abund))
  expect_identical(colnames(res), c("Sample", "raoQ"))

  unnamed <- unname(fx$abund)
  colnames(unnamed) <- colnames(fx$abund)
  expect_identical(divermeta(unnamed, fx$diss, "raoQ")$Sample, c("1", "2", "3"))
})


## Relative multiplicity and average redundancy

# A study with units of several subunits, some of them absent from some samples
make_unit_study <- function() {
  set.seed(7)
  n <- 14
  ids <- paste0("f", seq_len(n))
  clust <- stats::setNames(rep(c("A", "B", "C", "D"), c(5, 4, 3, 2)), ids)
  abund <- matrix(
    runif(6 * n, 1, 40) * (runif(6 * n) > 0.35),
    nrow = 6,
    dimnames = list(paste0("S", 1:6), ids)
  )
  abund[1, clust == "C"] <- 0 # a unit absent from a sample
  abund[6, ] <- 0 # an empty sample
  list(abund = abund, diss = mat_to_table(rand_diss(n, 0.05, 1.2), ids), clust = clust)
}

# Every combination of the flags of relative multiplicity
rm_flags <- expand.grid(
  cap_at_one = c(FALSE, TRUE),
  include_absent = c(FALSE, TRUE),
  assume_max_reference_distance = c(FALSE, TRUE),
  assume_homogeneous_abundance = c(FALSE, TRUE)
)


test_that("relative multiplicity and average redundancy match the direct functions", {
  st <- make_unit_study()
  sig <- 0.7

  for (k in seq_len(nrow(rm_flags))) {
    flags <- as.list(rm_flags[k, ])
    info <- paste(names(flags)[unlist(flags)], collapse = ", ")
    rm <- do.call(relative.multiplicity, c(list(st$abund, st$diss, st$clust, sigma = sig), flags))
    ar <- average.redundancy(st$abund, st$diss, st$clust, include_absent = flags$include_absent)

    # Alone, where only the pairs inside the units are read
    res <- do.call(divermeta, c(
      list(st$abund, st$diss, c("relative_multiplicity", "average_redundancy"), st$clust, sig = sig),
      flags
    ))
    expect_identical(colnames(res), c("Sample", "relative_multiplicity", "average_redundancy"))
    expect_equal(res$relative_multiplicity, unname(rm), info = info)
    expect_equal(res$average_redundancy, unname(ar), info = info)

    # With indices that read every pair
    res <- do.call(divermeta, c(
      list(st$abund, st$diss, c("raoQ", "AR", "M_distance", "RM"), st$clust, sig = sig, method = "average"),
      flags
    ))
    expect_equal(res$relative_multiplicity, unname(rm), info = info)
    expect_equal(res$average_redundancy, unname(ar), info = info)
    expect_equal(res$raoQ, unname(raoQuadratic(st$abund, st$diss)), info = info)
    expect_equal(
      res$multiplicity_distance,
      unname(multiplicity.distance(st$abund, st$diss, st$clust, method = "average", sig = sig)),
      info = info
    )
  }
})


test_that("subunits only listed in clust count in relative multiplicity and average redundancy", {
  st <- make_unit_study()
  sig <- 0.7

  # Leave out a subunit of unit A and the whole unit D from the abundances
  gone <- c("f2", names(st$clust)[st$clust == "D"])
  abund <- st$abund[, !(colnames(st$abund) %in% gone)]
  others <- c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq")

  for (k in seq_len(nrow(rm_flags))) {
    flags <- as.list(rm_flags[k, ])
    info <- paste(names(flags)[unlist(flags)], collapse = ", ")
    run <- function(expr) if (flags$include_absent) suppressWarnings(expr) else expr
    rm <- run(do.call(relative.multiplicity, c(list(abund, st$diss, st$clust, sigma = sig), flags)))
    ar <- average.redundancy(abund, st$diss, st$clust, include_absent = flags$include_absent)

    res <- run(do.call(divermeta, c(
      list(abund, st$diss, c(others, unit_indices), st$clust, sig = sig),
      flags
    )))
    expect_equal(res$relative_multiplicity, unname(rm), info = info)
    expect_equal(res$average_redundancy, unname(ar), info = info)

    # The other indices are those of the subunits in abund
    expect_equal(
      res[c("Sample", normalized)],
      divermeta(abund, st$diss, others, st$clust[colnames(abund)], sig = sig),
      info = info
    )
  }

  # The unit left out counts when absent units are included
  expect_false(isTRUE(all.equal(
    divermeta(abund, st$diss, "AR", st$clust, include_absent = TRUE)$average_redundancy,
    divermeta(abund, st$diss, "AR", st$clust[colnames(abund)], include_absent = TRUE)$average_redundancy
  )))

  # Units with no reference are dropped with a warning, as in the direct function
  expect_warning(
    divermeta(abund, st$diss, "RM", st$clust, include_absent = TRUE),
    "have no reference and are dropped: D"
  )
})


test_that("the homogeneous reference with listed distances needs every pair", {
  st <- make_unit_study()
  abund <- st$abund
  abund[, "f2"] <- 0
  partial <- st$diss[st$diss$ID1 != "f2" & st$diss$ID2 != "f2", ]

  # Pairs of a subunit absent from every sample may be missing, with a warning
  expect_warning(res <- divermeta(abund, partial, c("RM", "AR"), st$clust), "f2")
  expect_equal(res, divermeta(abund, st$diss, c("RM", "AR"), st$clust))

  expect_error(
    divermeta(abund, partial, c("RM", "AR"), st$clust, assume_homogeneous_abundance = TRUE),
    "Missing distances"
  )
  expect_warning(
    divermeta(
      abund, partial, c("RM", "AR"), st$clust,
      assume_homogeneous_abundance = TRUE, assume_max_reference_distance = TRUE
    ),
    "f2"
  )
})


test_that("the flags of relative multiplicity and average redundancy are validated", {
  fx <- make_fixture()
  flags <- c("include_absent", "cap_at_one", "assume_max_reference_distance", "assume_homogeneous_abundance")
  for (flag in flags) {
    args <- list(fx$abund, fx$diss, "RM", fx$clust)
    args[[flag]] <- "yes"
    expect_error(do.call(divermeta, args), paste0("`", flag, "` must be a single TRUE or FALSE value"))

    # Not used by the other indices
    args[[3]] <- "raoQ"
    expect_no_error(do.call(divermeta, args))
  }

  expect_error(
    divermeta(fx$abund, fx$diss, "AR", fx$clust, include_absent = NA),
    "`include_absent` must be a single TRUE or FALSE value"
  )
  expect_no_error(divermeta(fx$abund, fx$diss, "AR", fx$clust, cap_at_one = "yes"))
  expect_error(divermeta(fx$abund, fx$diss, "RM", fx$clust, sig = 0), "positive")
  expect_error(divermeta(fx$abund, fx$diss, "RM", fx$clust, sig = c(1, 1)))
})
