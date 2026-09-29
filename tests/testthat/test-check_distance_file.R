## Tests for the disk-based check of distance files


# Complete table of n subunits in units of `sizes`
check_study <- function(sizes = c(3, 4, 5)) {
  n <- sum(sizes)
  ids <- paste0("g", seq_len(n))
  clust <- stats::setNames(rep(paste0("U", seq_along(sizes)), sizes), ids)
  diss <- mat_to_table(rand_diss(n), ids)
  list(ids = ids, clust = clust, diss = diss, within = clust[diss$ID1] == clust[diss$ID2])
}


test_that("a complete file passes for every scope", {
  set.seed(3)
  st <- check_study()
  path <- write_diss(st$diss, "csv.gz")

  res <- check_distance_file(path, st$ids)
  expect_s3_class(res, "divermeta_distance_check")
  expect_true(res$ok)
  expect_equal(res$n_rows, 66)
  expect_equal(res$n_expected, 66)
  expect_equal(res$n_used, 66)
  expect_equal(res$n_missing, 0)
  expect_equal(res$n_duplicated, 0)
  expect_output(print(res), "OK")

  res <- check_distance_file(path, st$ids, st$clust, scope = "within")
  expect_true(res$ok)
  expect_equal(res$n_expected, 3 + 6 + 10)
  expect_equal(res$n_ignored, 66 - 19)

  res <- check_distance_file(path, st$ids, st$clust, scope = "between")
  expect_true(res$ok)
  expect_equal(res$n_expected, 66 - 19)

  # Only the pairs of the scope are needed
  expect_true(check_distance_file(write_diss(st$diss[st$within, ], "tsv"), st$ids, st$clust, "within")$ok)
  expect_false(check_distance_file(write_diss(st$diss[st$within, ], "tsv"), st$ids)$ok)
})


test_that("missing, duplicated, conflicting and invalid rows are counted", {
  set.seed(4)
  st <- check_study()
  diss <- st$diss
  diss$Distance <- as.character(diss$Distance)

  bad <- rbind(
    diss[-c(2, 10, 30), ], # 3 missing pairs
    transform(diss[5, ], ID1 = ID2, ID2 = ID1), # duplicated, other orientation
    diss[6, ], # duplicated
    transform(diss[7, ], Distance = "0.999"), # conflicting
    data.frame(ID1 = "g1", ID2 = "g12", Distance = "abc"), # non-numeric
    data.frame(ID1 = "g1", ID2 = "g11", Distance = "-0.5"), # negative
    data.frame(ID1 = c("x", "g3"), ID2 = c("g1", "g3"), Distance = "0.2") # ignored
  )
  path <- write_diss(bad, "tsv")

  for (max_pairs in c(1e7, 5, 1)) {
    res <- check_distance_file(path, st$ids, max_pairs_in_memory = max_pairs, chunk_size = 7)
    expect_false(res$ok)
    expect_equal(res$n_rows, 63 + 3 + 2 + 2)
    expect_equal(res$n_invalid, 2)
    expect_equal(res$n_ignored, 2)
    # The invalid rows are not used (g1-g12 and g1-g11 are listed again, validly)
    expect_equal(res$n_missing, 3)
    expect_equal(res$n_duplicated, 3)
    expect_equal(res$n_conflicting, 1)
    expect_false(is.null(res$examples$missing))
    expect_match(res$examples$conflicting, "0.999")
  }
  expect_output(print(res), "FAILED")

  # The example of a missing pair is a real missing pair
  res <- check_distance_file(write_diss(st$diss[-17, ], "csv"), st$ids)
  expect_equal(res$n_missing, 1)
  expect_equal(res$examples$missing, paste0(st$diss$ID1[17], "-", st$diss$ID2[17]))
})


test_that("missing pairs are found for every scope", {
  set.seed(5)
  st <- check_study()
  within <- which(st$within)
  between <- which(!st$within)
  path <- write_diss(st$diss[-c(within[2], between[3], between[4]), ], "csv")

  res <- check_distance_file(path, st$ids, st$clust, "within", max_pairs_in_memory = 2)
  expect_equal(res$n_missing, 1)
  expect_equal(res$examples$missing, paste0(st$diss$ID1[within[2]], "-", st$diss$ID2[within[2]]))

  res <- check_distance_file(path, st$ids, st$clust, "between", max_pairs_in_memory = 2)
  expect_equal(res$n_missing, 2)

  res <- check_distance_file(path, st$ids, max_pairs_in_memory = 2)
  expect_equal(res$n_missing, 3)
})


test_that("the temporary files are deleted", {
  set.seed(6)
  st <- check_study()
  tmp <- withr::local_tempdir()
  path <- write_diss(rbind(st$diss, st$diss[1:3, ]), "csv")

  res <- check_distance_file(path, st$ids, max_pairs_in_memory = 4, tmp_dir = tmp)
  expect_equal(res$n_duplicated, 3)
  expect_length(list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE), 0)
})


test_that("the arguments are validated", {
  set.seed(7)
  st <- check_study()
  path <- write_diss(st$diss, "csv")

  expect_error(check_distance_file(st$diss, st$ids), "path")
  expect_error(check_distance_file(tempfile(), st$ids), "not found")
  expect_error(check_distance_file(path, character(0)), "ids")
  expect_error(check_distance_file(path, c("a", "a")), "unique")
  expect_error(check_distance_file(path, st$ids, scope = "some"), "scope")
  expect_error(check_distance_file(path, st$ids, scope = "within"), "clust")
  expect_error(check_distance_file(path, st$ids, st$clust[1:3], scope = "within"), "clust")
  expect_error(check_distance_file(path, st$ids, max_pairs_in_memory = 0), "max_pairs_in_memory")
  expect_error(check_distance_file(path, st$ids, chunk_size = -1), "chunk_size")
  expect_error(check_distance_file(path, st$ids, tmp_dir = tempfile()), "tmp_dir")
})


test_that("the indices stop when the check fails", {
  set.seed(8)
  st <- check_study()
  abund <- matrix(
    runif(3 * length(st$ids), 1, 10),
    nrow = 3,
    dimnames = list(c("S1", "S2", "S3"), st$ids)
  )

  # A missing pair and a duplicated one: the number of rows is right
  cancelled <- rbind(st$diss[-4, ], st$diss[9, ])
  path <- write_diss(cancelled, "csv")

  # Without the check, the result is silently wrong
  res <- unchecked(raoQuadratic(abund, path))
  expect_false(isTRUE(all.equal(res, raoQuadratic(abund, st$diss))))

  # With the check, it stops
  for (fun in list(raoQuadratic, diversity.functional, redundancy, diversity.functional.traditional)) {
    expect_error(fun(abund, path, check_large_distance_file = TRUE), "Missing pairs: 1")
  }
  expect_error(
    divermeta(abund, path, c("raoQ", "M_distance"), st$clust, check_large_distance_file = TRUE),
    "did not pass the check"
  )

  # Inside units
  within <- which(st$within)
  path <- write_diss(rbind(st$diss[-within[1], ], st$diss[within[2], ]), "csv")
  expect_error(
    multiplicity.distance(abund, path, st$clust, check_large_distance_file = TRUE),
    "pairs of subunits of the same unit"
  )
  expect_error(
    relative.multiplicity(abund, path, st$clust, check_large_distance_file = TRUE),
    "Duplicated pairs: 1"
  )
  # ... which the linkage methods need too
  expect_error(
    multiplicity.distance(abund, path, st$clust, method = "average", check_large_distance_file = TRUE),
    "every pair of subunits"
  )
  # ... but not the distances between units
  expect_no_warning(unit_distances(path, st$clust, check_large_distance_file = TRUE))
})


test_that("several files are checked together, with rows numbered in each file", {
  set.seed(5)
  st <- check_study()
  diss <- st$diss
  diss$Distance <- as.character(diss$Distance)
  half <- seq_len(30)
  path1 <- write_diss(diss[half, ], "csv", shuffle = FALSE)
  path2 <- write_diss(diss[-half, ], "tsv", shuffle = FALSE)

  res <- check_distance_file(c(path1, path2), st$ids)
  expect_true(res$ok)
  expect_equal(res$n_rows, nrow(diss))
  expect_equal(res$n_used, nrow(diss))

  # The second row of the second file has an invalid distance
  bad2 <- diss[-half, ]
  bad2$Distance[2] <- "-1"
  res <- check_distance_file(c(path1, write_diss(bad2, "tsv", shuffle = FALSE)), st$ids)
  expect_false(res$ok)
  expect_equal(res$n_invalid, 1)
  expect_match(res$examples$invalid, "row 2 of ", fixed = TRUE)
  expect_match(res$examples$invalid, "\\.tsv\\)$")

  # A pair missing from the first file and listed twice in the second
  res <- check_distance_file(
    c(write_diss(diss[half[-1], ], "csv"), write_diss(rbind(diss[-half, ], diss[31, ]), "tsv")),
    st$ids
  )
  expect_equal(res$n_missing, 1)
  expect_equal(res$n_duplicated, 1)
})


test_that("a row with an invalid distance is invalid, not also missing", {
  set.seed(6)
  st <- check_study()
  diss <- st$diss
  diss$Distance <- as.character(diss$Distance)

  for (bad_value in c("NA", "-0.3", "Inf", "abc")) {
    bad <- diss
    bad$Distance[4] <- bad_value
    res <- check_distance_file(write_diss(bad, "tsv"), st$ids)
    expect_false(res$ok, info = bad_value)
    expect_equal(res$n_invalid, 1, info = bad_value)
    expect_equal(res$n_missing, 0, info = bad_value)
    expect_equal(res$n_duplicated, 0, info = bad_value)
    expect_equal(res$n_used, nrow(diss) - 1, info = bad_value)
  }
})
