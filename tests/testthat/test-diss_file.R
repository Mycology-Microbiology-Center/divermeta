## Tests for distances read from files (in chunks): every index must give the
## same result as with the distance table in memory


test_that("every file format is read like the table in memory", {
  st <- make_file_study()
  expected <- raoQuadratic(st$abund, st$diss)

  for (format in file_formats) {
    path <- write_diss(st$diss, format)
    for (chunk_size in c(1, 7, 1e6)) {
      res <- unchecked(raoQuadratic(st$abund, path, chunk_size = chunk_size, header = has_header(format)))
      expect_equal(res, expected, info = paste(format, chunk_size))
    }
  }
})


test_that("zstd files are read through the command line tool", {
  skip_if(!nzchar(Sys.which("zstd")), "zstd is not installed")
  st <- make_file_study()
  path <- write_diss(st$diss, "tsv")
  zst <- paste0(path, ".zst")
  system2("zstd", c("-q", shQuote(path), "-o", shQuote(zst)))
  on.exit(unlink(zst), add = TRUE)

  expect_equal(unchecked(raoQuadratic(st$abund, zst, chunk_size = 7)), raoQuadratic(st$abund, st$diss))
})


test_that("several files are read one after the other", {
  st <- make_file_study()
  half <- seq_len(nrow(st$diss)) <= nrow(st$diss) / 2
  paths <- c(write_diss(st$diss[half, ], "csv"), write_diss(st$diss[!half, ], "tsv.gz"))

  expect_equal(
    unchecked(redundancy(st$abund, paths, chunk_size = 5)),
    redundancy(st$abund, st$diss)
  )
  expect_equal(
    redundancy(st$abund, paths, chunk_size = 5, check_large_distance_file = TRUE),
    redundancy(st$abund, st$diss)
  )
})


test_that("diversity indices from a file equal the in-memory ones", {
  st <- make_file_study()
  path <- write_diss(st$diss, "csv.gz")

  for (sig in c(0.3, 0.8, 2)) {
    expect_equal(
      unchecked(diversity.functional(st$abund, path, sig = sig, chunk_size = 7)),
      diversity.functional(st$abund, st$diss, sig = sig)
    )
  }
  for (q in c(0, 0.5, 1, 2, 3)) {
    expect_equal(
      unchecked(diversity.functional.traditional(st$abund, path, q = q, chunk_size = 7)),
      diversity.functional.traditional(st$abund, st$diss, q = q)
    )
  }
  expect_equal(
    unchecked(redundancy(st$abund, path, chunk_size = 7)),
    redundancy(st$abund, st$diss)
  )

  # Sparse abundances
  sparse <- Matrix::Matrix(st$abund, sparse = TRUE)
  expect_equal(
    unchecked(raoQuadratic(sparse, path, chunk_size = 7)),
    raoQuadratic(st$abund, st$diss)
  )
})


test_that("distance-based multiplicity from a file equals the in-memory one", {
  st <- make_file_study()
  path <- write_diss(st$diss, "tsv")

  for (method in c("sigma", "min", "max", "average")) {
    for (sig in c(0.4, 1)) {
      expect_equal(
        unchecked(multiplicity.distance(st$abund, path, st$clust, method = method, sig = sig, chunk_size = 9)),
        multiplicity.distance(st$abund, st$diss, st$clust, method = method, sig = sig),
        info = paste(method, sig)
      )
    }
  }

  diss_clust <- unit_distances(st$diss, st$clust, method = "average")
  expect_equal(
    unchecked(multiplicity.distance(st$abund, path, st$clust, method = "custom", diss_clust = diss_clust, chunk_size = 9)),
    multiplicity.distance(st$abund, st$diss, st$clust, method = "custom", diss_clust = diss_clust)
  )

  # With sigma, a file with only the pairs inside units is enough
  within <- st$clust[st$diss$ID1] == st$clust[st$diss$ID2]
  path_within <- write_diss(st$diss[within, ], "csv")
  expect_equal(
    unchecked(multiplicity.distance(st$abund, path_within, st$clust, sig = 0.7, chunk_size = 4)),
    multiplicity.distance(st$abund, st$diss, st$clust, sig = 0.7)
  )
})


test_that("multiplicity from a file equals the legacy implementation", {
  st <- make_file_study()
  path <- write_diss(st$diss, "csv")
  ids <- colnames(st$abund)
  diss_mtx <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  diss_mtx[cbind(st$diss$ID1, st$diss$ID2)] <- st$diss$Distance
  diss_mtx <- diss_mtx + t(diss_mtx)

  for (method in c("min", "max", "average")) {
    res <- unchecked(multiplicity.distance(st$abund, path, st$clust, method = method, sig = 0.8, chunk_size = 6))
    # Legacy is given a dist object: it gives wrong cluster distances for full matrices
    expected <- apply(st$abund[1:5, ], 1, function(ab) {
      legacy$multiplicity.distance(ab, stats::as.dist(diss_mtx), unname(st$clust), method = method, sig = 0.8)
    })
    expect_equal(res[1:5], expected, info = method)
  }
})


test_that("relative multiplicity from a file equals the in-memory one", {
  st <- make_file_study()
  path <- write_diss(st$diss, "tsv.gz")

  for (max_dist in c(FALSE, TRUE)) {
    for (homogeneous in c(FALSE, TRUE)) {
      for (include_absent in c(FALSE, TRUE)) {
        expect_equal(
          unchecked(relative.multiplicity(
            st$abund, path, st$clust,
            sigma = 0.9,
            include_absent = include_absent,
            assume_max_reference_distance = max_dist,
            assume_homogeneous_abundance = homogeneous,
            chunk_size = 8
          )),
          relative.multiplicity(
            st$abund, st$diss, st$clust,
            sigma = 0.9,
            include_absent = include_absent,
            assume_max_reference_distance = max_dist,
            assume_homogeneous_abundance = homogeneous
          )
        )
      }
    }
  }

  # Per-unit sigma and a clustering with subunits missing from abund
  sigma <- c(A = 0.5, B = 0.8, C = 1, D = 0.3)
  clust <- c(st$clust, extra1 = "B", extra2 = "B")
  diss <- rbind(
    st$diss,
    data.frame(ID1 = "extra1", ID2 = c(names(st$clust)[st$clust == "B"], "extra2"), Distance = 0.4),
    data.frame(ID1 = "extra2", ID2 = names(st$clust)[st$clust == "B"], Distance = 0.6)
  )
  path <- write_diss(diss, "csv")
  expect_equal(
    unchecked(relative.multiplicity(st$abund, path, clust, sigma = sigma, chunk_size = 5)),
    relative.multiplicity(st$abund, diss, clust, sigma = sigma)
  )

  ref_div <- c(A = 2, B = 3, C = 2.5, D = 1.5)
  expect_equal(
    unchecked(relative.multiplicity.ref_div(st$abund, path, clust, ref_div, sigma = sigma, chunk_size = 5)),
    relative.multiplicity.ref_div(st$abund, diss, clust, ref_div, sigma = sigma)
  )
})


test_that("average redundancy from a file equals the in-memory one", {
  st <- make_file_study()
  path <- write_diss(st$diss, "tsv.gz")

  for (normalize in c(FALSE, TRUE)) {
    for (include_absent in c(FALSE, TRUE)) {
      expected <- average.redundancy(
        st$abund, st$diss, st$clust,
        normalize = normalize, include_absent = include_absent
      )
      expect_equal(
        unchecked(average.redundancy(
          st$abund, path, st$clust,
          normalize = normalize, include_absent = include_absent, chunk_size = 8
        )),
        expected
      )
      expect_equal(
        average.redundancy(
          st$abund, path, st$clust,
          normalize = normalize, include_absent = include_absent,
          chunk_size = 8, check_large_distance_file = TRUE
        ),
        expected
      )
    }
  }
})


test_that("MAD from a file equals the in-memory one and needs no check", {
  st <- make_file_study()
  path <- write_diss(st$diss, "csv")

  # MAD finds its own missing pairs: no warning
  expect_no_warning(res <- metagenomic.alpha.index(st$abund, path, st$clust, chunk_size = 7))
  expect_equal(res, metagenomic.alpha.index(st$abund, st$diss, st$clust))

  reps <- tapply(names(st$clust), st$clust, function(x) x[2])
  expect_equal(
    metagenomic.alpha.index(st$abund, path, st$clust, representatives = reps, chunk_size = 7),
    metagenomic.alpha.index(st$abund, st$diss, st$clust, representatives = reps)
  )

  # Only the pairs with the representatives are needed
  rep_pairs <- st$diss$ID1 %in% reps | st$diss$ID2 %in% reps
  path_reps <- write_diss(st$diss[rep_pairs, ], "tsv")
  expect_equal(
    metagenomic.alpha.index(st$abund, path_reps, st$clust, representatives = reps, chunk_size = 3),
    metagenomic.alpha.index(st$abund, st$diss, st$clust, representatives = reps)
  )

  # A missing pair, and a pair listed again with a different distance
  expect_error(
    metagenomic.alpha.index(st$abund, path_reps, st$clust, chunk_size = 3),
    "Missing distances"
  )
  # (pairs MAD does not use, e.g. between units, are ignored)
  k <- which(rep_pairs & st$clust[st$diss$ID1] == st$clust[st$diss$ID2])[1]
  conflicting <- rbind(st$diss, transform(st$diss[k, ], Distance = Distance + 0.1))
  expect_error(
    metagenomic.alpha.index(st$abund, write_diss(conflicting, "csv"), st$clust, representatives = reps, chunk_size = 3),
    "conflicting"
  )
  # The same distance again does not change the result
  repeated <- rbind(st$diss, st$diss[k, ])
  expect_equal(
    metagenomic.alpha.index(st$abund, write_diss(repeated, "csv"), st$clust, representatives = reps),
    metagenomic.alpha.index(st$abund, st$diss, st$clust, representatives = reps)
  )
})


test_that("unit distances from a file equal the in-memory ones", {
  st <- make_file_study()
  path <- write_diss(st$diss, "tsv")

  for (method in c("min", "max", "average")) {
    expect_equal(
      unchecked(unit_distances(path, st$clust, method = method, chunk_size = 6)),
      unit_distances(st$diss, st$clust, method = method)
    )
  }

  # Sigma does not read the distances
  expect_no_warning(res <- unit_distances(path, st$clust, method = "sigma", sig = 0.5))
  expect_equal(res, unit_distances(st$diss, st$clust, method = "sigma", sig = 0.5))

  # Only the pairs between units are needed and checked
  between <- st$clust[st$diss$ID1] != st$clust[st$diss$ID2]
  path_between <- write_diss(st$diss[between, ], "csv")
  expect_equal(
    unit_distances(path_between, st$clust, chunk_size = 6, check_large_distance_file = TRUE),
    unit_distances(st$diss, st$clust)
  )
})


test_that("unit distances from a file stop when a pair of units has no distance", {
  # A-B is listed twice and A-C is missing, so the row count still matches
  cl <- c(a1 = "A", b1 = "B", c1 = "C")
  path <- write_diss(
    data.frame(ID1 = c("a1", "a1", "b1"), ID2 = c("b1", "b1", "c1"), Distance = c(0.5, 0.5, 0.7)),
    "tsv",
    shuffle = FALSE
  )
  for (method in c("min", "max", "average")) {
    expect_error(
      suppressWarnings(unit_distances(path, cl, method)),
      "No distances in .diss. between the subunits"
    )
  }
})


test_that("unit distances with sigma do not need the distances", {
  cl <- c(a1 = "A", a2 = "A", b1 = "B", c1 = "C")
  res <- unit_distances(NULL, cl, "sigma", sig = 0.4)
  expect_equal(res$ID1, c("A", "A", "B"))
  expect_equal(res$ID2, c("B", "C", "C"))
  expect_equal(res$Distance, rep(0.4, 3))
})


test_that("divermeta reads the file once for all indices", {
  st <- make_file_study()
  path <- write_diss(st$diss, "csv.gz")
  indices <- c("M_inventory", "raoQ", "FD_sigma", "redundancy", "M_distance", "FDq")

  opened <- 0
  local_mocked_bindings(.open_diss_file = function(path) {
    opened <<- opened + 1
    gzfile(path, "r")
  })

  for (method in c("sigma", "average")) {
    opened <- 0
    res <- unchecked(divermeta(st$abund, path, indices, st$clust, q = 2, sig = 0.7, method = method, chunk_size = 11))
    expect_equal(res, divermeta(st$abund, st$diss, indices, st$clust, q = 2, sig = 0.7, method = method))
    # One opening to check the columns, one to read the rows
    expect_equal(opened, 2)
  }

  # With the check, the file is read once more
  opened <- 0
  res <- divermeta(st$abund, path, indices, st$clust, chunk_size = 11, check_large_distance_file = TRUE)
  expect_equal(res, divermeta(st$abund, st$diss, indices, st$clust))
  expect_equal(opened, 3)

  # Relative multiplicity reads the samples and the reference together
  opened <- 0
  unchecked(relative.multiplicity(st$abund, path, st$clust, chunk_size = 11))
  expect_equal(opened, 2)
})


test_that("the warning is given only for unchecked files", {
  st <- make_file_study()
  path <- write_diss(st$diss, "csv")

  expect_warning(raoQuadratic(st$abund, path), "check_large_distance_file")
  expect_no_warning(raoQuadratic(st$abund, path, check_large_distance_file = TRUE))
  expect_no_warning(raoQuadratic(st$abund, st$diss))
  expect_no_warning(raoQuadratic(st$abund, st$diss, check_large_distance_file = TRUE))

  # divermeta warns once
  w <- testthat::capture_warnings(
    divermeta(st$abund, path, c("raoQ", "FD_sigma", "redundancy"))
  )
  expect_length(w, 1)
})


test_that("problems with the file are errors", {
  st <- make_file_study(n = 8)
  ab <- st$abund

  expect_error(raoQuadratic(ab, tempfile(fileext = ".csv")), "not found")
  expect_error(raoQuadratic(ab, tempdir()), "not found")
  expect_error(raoQuadratic(ab, character(0)), "path")
  expect_error(raoQuadratic(ab, write_diss(st$diss, "csv"), chunk_size = 0), "chunk_size")
  expect_error(
    raoQuadratic(ab, write_diss(st$diss, "csv"), check_large_distance_file = NA),
    "check_large_distance_file"
  )

  # Fewer than three columns
  two <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("ID1,ID2", "f1,f2"), two)
  expect_error(raoQuadratic(ab, two), "three columns")

  # Invalid distances
  bad <- st$diss
  bad$Distance <- as.character(bad$Distance)
  bad$Distance[3] <- "abc"
  expect_error(unchecked(raoQuadratic(ab, write_diss(bad, "tsv"))), "numeric")
  bad$Distance[3] <- "NA"
  expect_error(unchecked(raoQuadratic(ab, write_diss(bad, "tsv"))), "NA distances")
  bad$Distance[3] <- "-0.2"
  expect_error(unchecked(raoQuadratic(ab, write_diss(bad, "tsv"))), "non-negative")

  # No identifier of abund
  other <- transform(st$diss, ID1 = paste0("x", ID1), ID2 = paste0("x", ID2))
  expect_error(unchecked(raoQuadratic(ab, write_diss(other, "csv"))), "No row of the distance file")
})


test_that("the number of rows must match the number of pairs", {
  st <- make_file_study(n = 12)
  # Every subunit present in some sample, so every pair is needed
  ab <- st$abund + 1

  # A missing pair
  expect_error(
    unchecked(raoQuadratic(ab, write_diss(st$diss[-5, ], "csv"))),
    "Missing distances in `diss` for 1 pair"
  )
  # A duplicated pair (e.g. in both orientations)
  flipped <- transform(st$diss[5, ], ID1 = ID2, ID2 = ID1)
  expect_error(
    unchecked(raoQuadratic(ab, write_diss(rbind(st$diss, flipped), "csv", shuffle = FALSE))),
    "more than once"
  )
  # Extra rows with unknown identifiers and self pairs are ignored
  extra <- rbind(
    st$diss,
    data.frame(ID1 = c("x1", "f1"), ID2 = c("f1", "f1"), Distance = c(0.3, 0))
  )
  expect_equal(unchecked(raoQuadratic(ab, write_diss(extra, "csv"))), raoQuadratic(ab, st$diss))

  # Pairs inside units for sigma and relative multiplicity
  within <- which(st$clust[st$diss$ID1] == st$clust[st$diss$ID2])
  path <- write_diss(st$diss[-within[1], ], "csv")
  expect_error(unchecked(multiplicity.distance(ab, path, st$clust)), "in unit")
  expect_error(unchecked(relative.multiplicity(ab, path, st$clust)), "in unit")
  path <- write_diss(rbind(st$diss, st$diss[within[1], ]), "csv")
  expect_error(unchecked(multiplicity.distance(ab, path, st$clust)), "more than once")

  # A linkage method with no distance between two units (hidden by a duplicate)
  clust <- c(st$clust[1:11], f12 = "E")
  e_pairs <- which(st$diss$ID1 == "f12" | st$diss$ID2 == "f12")
  hidden <- rbind(st$diss[-e_pairs, ], st$diss[rep(1, length(e_pairs)), ])
  expect_error(
    unchecked(multiplicity.distance(ab, write_diss(hidden, "csv"), clust, method = "min")),
    "No distances in `diss` between the subunits"
  )
  # ... which the check finds
  expect_error(
    multiplicity.distance(ab, write_diss(hidden, "csv"), clust, method = "min", check_large_distance_file = TRUE),
    "did not pass the check"
  )
})


# Rows outside the scope
# ----------------------------

test_that("invalid distances between units are ignored by the indices that only use pairs inside units", {
  st <- make_file_study(n = 16)
  ab <- st$abund
  within <- st$clust[st$diss$ID1] == st$clust[st$diss$ID2]
  between <- which(!within)

  # An NA and a negative distance between units, and a pair between units
  # listed again with a different distance
  bad <- st$diss
  bad$Distance[between[1:2]] <- c(NA, -0.5)
  bad <- rbind(bad, transform(st$diss[between[3], ], Distance = Distance + 0.1))

  expected <- list(
    md = multiplicity.distance(ab, st$diss, st$clust, sig = 0.7),
    rm = relative.multiplicity(ab, st$diss, st$clust, sigma = 0.7),
    ar = average.redundancy(ab, st$diss, st$clust),
    dm = divermeta(ab, st$diss, c("M_inventory", "M_distance"), st$clust, sig = 0.7)
  )

  # In memory
  expect_equal(multiplicity.distance(ab, bad, st$clust, sig = 0.7), expected$md)
  expect_equal(relative.multiplicity(ab, bad, st$clust, sigma = 0.7), expected$rm)
  expect_equal(average.redundancy(ab, bad, st$clust), expected$ar)
  expect_equal(divermeta(ab, bad, c("M_inventory", "M_distance"), st$clust, sig = 0.7), expected$dm)

  # From a file (without the duplicated row, which the row count would reject),
  # with and without the check
  path <- write_diss(bad[-nrow(bad), ], "csv")
  expect_true(check_distance_file(path, colnames(ab), st$clust, "within")$ok)
  for (check in c(FALSE, TRUE)) {
    run <- function(expr) if (check) expr else unchecked(expr)
    expect_equal(
      run(multiplicity.distance(ab, path, st$clust, sig = 0.7, chunk_size = 5, check_large_distance_file = check)),
      expected$md
    )
    expect_equal(
      run(relative.multiplicity(ab, path, st$clust, sigma = 0.7, chunk_size = 5, check_large_distance_file = check)),
      expected$rm
    )
    expect_equal(
      run(average.redundancy(ab, path, st$clust, chunk_size = 5, check_large_distance_file = check)),
      expected$ar
    )
    expect_equal(
      run(divermeta(ab, path, c("M_inventory", "M_distance"), st$clust, sig = 0.7, chunk_size = 5,
        check_large_distance_file = check
      )),
      expected$dm
    )
  }

  # Indices that use every pair still stop
  expect_error(raoQuadratic(ab, bad), "NA distances")
  expect_error(multiplicity.distance(ab, bad, st$clust, method = "average"), "NA distances")
  expect_error(divermeta(ab, bad, c("raoQ", "M_distance"), st$clust), "NA distances")
  expect_error(unchecked(raoQuadratic(ab, path)), "NA distances")
  expect_error(raoQuadratic(ab, path, check_large_distance_file = TRUE), "Invalid distances: 2")

  # Invalid distances inside units still stop
  bad_in <- st$diss
  bad_in$Distance[which(within)[1]] <- NA
  expect_error(multiplicity.distance(ab, bad_in, st$clust), "NA distances")
  expect_error(relative.multiplicity(ab, bad_in, st$clust), "NA distances")
  bad_in$Distance[which(within)[1]] <- -1
  expect_error(average.redundancy(ab, bad_in, st$clust), "non-negative")
  expect_error(unchecked(average.redundancy(ab, write_diss(bad_in, "tsv"), st$clust)), "non-negative")
})


test_that("invalid distances inside units are ignored by unit_distances", {
  st <- make_file_study(n = 16)
  within <- which(st$clust[st$diss$ID1] == st$clust[st$diss$ID2])
  bad <- st$diss
  bad$Distance[within[1:2]] <- c(NA, -0.5)

  for (method in c("min", "max", "average")) {
    expected <- unit_distances(st$diss, st$clust, method = method)
    expect_equal(unit_distances(bad, st$clust, method = method), expected)
    expect_equal(
      unit_distances(write_diss(bad, "csv"), st$clust, method = method, check_large_distance_file = TRUE),
      expected
    )
  }
})


test_that("a file whose rows are all outside the scope is accepted when no pair is needed", {
  # Every unit has a single subunit: sigma multiplicity needs no pair
  ids <- c("a", "b", "c")
  clust <- c(a = "A", b = "B", c = "C")
  ab <- one_row(c(1, 2, 3), ids)
  diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.2, NA, 0.4))
  path <- write_diss(diss, "csv")

  expect_equal(unname(multiplicity.distance(ab, diss, clust)), 1)
  expect_equal(unname(unchecked(multiplicity.distance(ab, path, clust))), 1)
})


# Header of the files
# ----------------------------

test_that("files are expected to have a header unless header = FALSE", {
  st <- make_file_study(n = 10)
  ab <- st$abund
  expected <- raoQuadratic(ab, st$diss)
  with_header <- write_diss(st$diss, "csv")
  no_header <- write_diss(st$diss, "tab")

  # Default: header expected
  expect_equal(unchecked(raoQuadratic(ab, with_header)), expected)
  expect_error(raoQuadratic(ab, no_header), "looks like data.*header = FALSE")
  expect_error(check_distance_file(no_header, colnames(ab)), "header = FALSE")
  expect_error(metagenomic.alpha.index(ab, no_header, st$clust), "header = FALSE")
  expect_error(divermeta(ab, no_header, c("raoQ", "M_inventory"), st$clust), "header = FALSE")
  expect_error(unit_distances(no_header, st$clust), "header = FALSE")

  # header = FALSE reads the first line as data
  expect_equal(unchecked(raoQuadratic(ab, no_header, header = FALSE)), expected)
  expect_equal(raoQuadratic(ab, no_header, header = FALSE, check_large_distance_file = TRUE), expected)
  expect_true(check_distance_file(no_header, colnames(ab), header = FALSE)$ok)
  expect_equal(
    metagenomic.alpha.index(ab, no_header, st$clust, header = FALSE),
    metagenomic.alpha.index(ab, st$diss, st$clust)
  )
  expect_equal(
    unchecked(divermeta(ab, no_header, c("raoQ", "M_distance"), st$clust, header = FALSE)),
    divermeta(ab, st$diss, c("raoQ", "M_distance"), st$clust)
  )
  expect_equal(
    unchecked(unit_distances(no_header, st$clust, header = FALSE)),
    unit_distances(st$diss, st$clust)
  )
  expect_equal(
    unchecked(relative.multiplicity(ab, no_header, st$clust, header = FALSE)),
    relative.multiplicity(ab, st$diss, st$clust)
  )
  expect_equal(
    unchecked(average.redundancy(ab, no_header, st$clust, header = FALSE)),
    average.redundancy(ab, st$diss, st$clust)
  )

  # ... and a header line is then an invalid row
  expect_error(unchecked(raoQuadratic(ab, with_header, header = FALSE)), "numeric")

  # Compressed file without a header, split in two files
  half <- seq_len(nrow(st$diss)) <= nrow(st$diss) / 2
  gz <- withr::local_tempfile(fileext = ".csv.gz")
  con <- gzfile(gz, "w")
  utils::write.table(st$diss[half, ], con, sep = ",", row.names = FALSE, col.names = FALSE)
  close(con)
  plain <- write_diss(st$diss[!half, ], "txt")
  expect_equal(unchecked(raoQuadratic(ab, c(gz, plain), header = FALSE, chunk_size = 3)), expected)

  # The flag is validated
  expect_error(raoQuadratic(ab, with_header, header = NA), "`header`")
  expect_error(check_distance_file(with_header, colnames(ab), header = "yes"), "`header`")
})


test_that("a first row with an invalid distance is never dropped silently", {
  ids <- c("a1", "a2", "b1", "b2")
  ab <- one_row(c(1, 2, 3, 4), ids)
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("a1,a2,oops", "a1,b1,0.5", "a1,b2,0.5", "a2,b1,0.5", "a2,b2,0.5", "b1,b2,0.3"), path)

  # Read as a header, it lacks the pair a1-a2; read as data, its distance is invalid
  expect_error(unchecked(raoQuadratic(ab, path)), "Missing distances.*1 pair")
  expect_error(unchecked(raoQuadratic(ab, path, header = FALSE)), "numeric")
  res <- check_distance_file(path, ids, header = FALSE)
  expect_equal(res$n_invalid, 1)
  expect_match(res$examples$invalid, "a1-a2 \\(distance: oops, row 1")
})
