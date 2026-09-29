## Tests for subunits with zero abundance in every sample, which need no distances


# Units A (a1, a2, a3), B (b1, b2, b3) and C (c1). a3 and b3 have zero abundance
# in every sample; the others are present in at least one
absent_study <- function() {
  ids <- c("a1", "a2", "a3", "b1", "b2", "b3", "c1")
  abund <- matrix(
    c(
      4, 1, 0, 2, 3, 0, 1,
      0, 5, 0, 1, 0, 0, 2,
      2, 2, 0, 0, 6, 0, 0
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), ids)
  )
  set.seed(11)
  diss <- mat_to_table(rand_diss(length(ids)), ids)
  absent <- c("a3", "b3")
  with_absent <- diss$ID1 %in% absent | diss$ID2 %in% absent
  present <- setdiff(ids, absent)
  list(
    abund = abund,
    clust = c(a1 = "A", a2 = "A", a3 = "A", b1 = "B", b2 = "B", b3 = "B", c1 = "C"),
    diss = diss,
    # Without any pair of an absent subunit, and without some of them
    diss_none = diss[!with_absent, ],
    diss_some = diss[-which(with_absent)[c(1, 3)], ],
    present = present
  )
}

# Value and warnings of an expression
collect_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

absent_warning <- "zero abundance in every sample"

# Every index that reads the distances, as a function of abund, diss and clust
absent_indices <- list(
  raoQ = function(a, d, cl, ...) raoQuadratic(a, d, ...),
  FD_sigma = function(a, d, cl, ...) diversity.functional(a, d, sig = 0.7, ...),
  FD_q = function(a, d, cl, ...) diversity.functional.traditional(a, d, q = 2, ...),
  redundancy = function(a, d, cl, ...) redundancy(a, d, ...),
  M_sigma = function(a, d, cl, ...) multiplicity.distance(a, d, cl, sig = 0.8, ...),
  RM = function(a, d, cl, ...) relative.multiplicity(a, d, cl, ...),
  RM_max_ref = function(a, d, cl, ...) {
    relative.multiplicity(
      a, d, cl,
      assume_max_reference_distance = TRUE, assume_homogeneous_abundance = TRUE, ...
    )
  },
  RM_ref_div = function(a, d, cl, ...) {
    relative.multiplicity.ref_div(a, d, cl, ref_div = c(A = 2, B = 1.5, C = 1), ...)
  },
  AR = function(a, d, cl, ...) average.redundancy(a, d, cl, ...),
  divermeta = function(a, d, cl, ...) {
    divermeta(a, d, c("raoQ", "FD_sigma", "redundancy", "multiplicity_distance"), cl, sig = 0.8, ...)
  }
)


test_that("missing distances of absent subunits give a warning and change nothing", {
  st <- absent_study()
  for (name in names(absent_indices)) {
    f <- absent_indices[[name]]

    # Complete table: no warning
    expect_no_warning(expected <- f(st$abund, st$diss, st$clust))
    # Same as without the absent subunits, except for the homogeneous reference,
    # which counts every subunit of the unit
    if (name != "RM_max_ref") {
      expect_equal(expected, f(st$abund[, st$present], st$diss, st$clust[st$present]), info = name)
    }

    for (diss in list(st$diss_none, st$diss_some)) {
      # In memory
      expect_warning(res <- f(st$abund, diss, st$clust), absent_warning, info = name)
      expect_equal(res, expected, info = name)

      # In a file, unchecked and checked
      path <- write_diss(diss, "csv")
      out <- collect_warnings(f(st$abund, path, st$clust))
      expect_true(any(grepl("not checked", out$warnings)), info = name)
      expect_true(any(grepl(absent_warning, out$warnings)), info = name)
      expect_equal(out$value, expected, info = name)

      expect_warning(
        res <- f(st$abund, path, st$clust, check_large_distance_file = TRUE),
        absent_warning, info = name
      )
      expect_equal(res, expected, info = name)
    }

    # A complete file gives only the "not checked" warning
    out <- collect_warnings(f(st$abund, write_diss(st$diss, "csv"), st$clust))
    expect_false(any(grepl(absent_warning, out$warnings)), info = name)
  }
})


test_that("missing distances between present subunits are still errors", {
  st <- absent_study()
  # Drop a pair of present subunits of the same unit
  diss <- st$diss_none
  diss <- diss[-which(st$clust[diss$ID1] == st$clust[diss$ID2])[1], ]

  expect_error(raoQuadratic(st$abund, diss), "Missing distances")
  expect_error(multiplicity.distance(st$abund, diss, st$clust), "Missing distances")
  expect_error(relative.multiplicity(st$abund, diss, st$clust), "Missing distances")
  expect_error(unchecked(raoQuadratic(st$abund, write_diss(diss, "csv"))), "Missing distances")
})


test_that("the homogeneous reference of relative multiplicity still needs every pair", {
  st <- absent_study()
  expect_error(
    relative.multiplicity(st$abund, st$diss_none, st$clust, assume_homogeneous_abundance = TRUE),
    "Missing distances"
  )
  expect_no_warning(
    relative.multiplicity(st$abund, st$diss, st$clust, assume_homogeneous_abundance = TRUE)
  )
})


test_that("linkage methods use the listed distances of absent subunits", {
  st <- absent_study()
  # a3 is very close to b1: it sets the minimum distance between A and B
  diss <- st$diss
  close <- (diss$ID1 == "a3" & diss$ID2 == "b1") | (diss$ID1 == "b1" & diss$ID2 == "a3")
  diss$Distance[close] <- 0.001

  with_a3 <- multiplicity.distance(st$abund, diss, st$clust, method = "min")
  without_a3 <- multiplicity.distance(
    st$abund[, st$present], st$diss, st$clust[st$present], method = "min"
  )
  expect_false(isTRUE(all.equal(with_a3, without_a3)))
  expect_warning(
    res <- multiplicity.distance(st$abund, st$diss_none, st$clust, method = "min"),
    absent_warning
  )
  expect_equal(res, without_a3)

  # A unit whose subunits are all absent needs no distance to the other units
  clust <- c(st$clust, d1 = "D")
  abund <- cbind(st$abund, d1 = 0)
  expect_warning(
    res <- multiplicity.distance(abund, st$diss_none, clust, method = "average"),
    absent_warning
  )
  expect_equal(
    res,
    multiplicity.distance(st$abund[, st$present], st$diss, st$clust[st$present], method = "average")
  )
})


test_that("a file listing pairs of absent subunits too many times is an error", {
  st <- absent_study()
  with_absent <- which(st$diss$ID1 == "a3" | st$diss$ID2 == "a3")
  path <- write_diss(rbind(st$diss, st$diss[with_absent[1], ]), "csv")
  expect_error(unchecked(raoQuadratic(st$abund, path)), "more than once")
})
