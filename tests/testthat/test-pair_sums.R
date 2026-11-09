## Tests for the internal .pair_sums


# Former implementation, pair by pair, as the reference. With `group` (the
# group of every pair), the sums are kept apart for every group
pair_sums_reference <- function(X, i, j, w, Y = X, group = NULL, n_groups = 1L) {
  out <- matrix(0, nrow = nrow(X), ncol = if (is.null(group)) 1L else n_groups)
  if (length(w) == 0 || nrow(X) == 0) {
    return(out)
  }
  prod <- X[, i, drop = FALSE] * Y[, j, drop = FALSE]
  if (is.null(group)) {
    out[, 1] <- as.vector(as.matrix(prod %*% w))
  } else {
    W <- Matrix::sparseMatrix(i = seq_along(w), j = group, x = w, dims = c(length(w), n_groups))
    out <- as.matrix(prod %*% W)
  }
  out
}

# Sparse samples x subunits table, with subunits present in few samples
sparse_table <- function(n, p, fill = 0.15) {
  Matrix::rsparsematrix(n, p, density = fill, rand.x = function(k) runif(k, 0.1, 5))
}

# Random pairs of distinct subunits; with `unit_of`, only pairs inside units
random_pairs <- function(p, m, unit_of = NULL) {
  i <- sample.int(p, m, replace = TRUE)
  j <- sample.int(p, m, replace = TRUE)
  if (!is.null(unit_of)) {
    j <- vapply(i, function(k) {
      members <- which(unit_of == unit_of[k])
      members[sample.int(length(members), 1)]
    }, integer(1))
  }
  ok <- i != j
  list(i = i[ok], j = j[ok], w = runif(sum(ok)))
}

# All pairs inside every unit, sorted by the first subunit as in the distance
# files
sorted_unit_pairs <- function(unit_of) {
  idx <- which(outer(unit_of, unit_of, "==") & upper.tri(diag(length(unit_of))), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  list(i = idx[, 1], j = idx[, 2], w = runif(nrow(idx)))
}

expect_same_sums <- function(X, i, j, w, Y = X, unit_of = NULL, n_units = 1L) {
  if (is.null(unit_of)) {
    expect_equal(.pair_sums(X, i, j, w, Y = Y), pair_sums_reference(X, i, j, w, Y = Y), tolerance = 1e-12)
  } else {
    expect_equal(
      .pair_sums(X, i, j, w, Y = Y, col_group = unit_of, n_groups = n_units),
      pair_sums_reference(X, i, j, w, Y = Y, group = unit_of[i], n_groups = n_units),
      tolerance = 1e-12
    )
  }
}

make_case <- function() {
  set.seed(42)
  n <- 60
  p <- 300
  unit_of <- rep(1:15, each = 20)
  list(
    X = sparse_table(n, p), Y = sparse_table(n, p, 0.3), unit_of = unit_of, n_units = 15L,
    pairs = random_pairs(p, 3000), within = random_pairs(p, 3000, unit_of)
  )
}


test_that("matches the pair by pair sums, with and without groups", {
  cs <- make_case()
  pr <- cs$pairs
  wi <- cs$within
  expect_same_sums(cs$X, pr$i, pr$j, pr$w)
  expect_same_sums(cs$X, wi$i, wi$j, wi$w, unit_of = cs$unit_of, n_units = cs$n_units)
  # Y different from X, and the pairs swapped
  expect_same_sums(cs$X, pr$i, pr$j, pr$w, Y = cs$Y)
  expect_same_sums(cs$X, pr$j, pr$i, pr$w, Y = cs$Y)
  expect_same_sums(cs$X, wi$j, wi$i, wi$w, Y = cs$Y, unit_of = cs$unit_of, n_units = cs$n_units)
})


test_that("repeated pairs, zero and negative weights add up as in a sum", {
  cs <- make_case()
  pr <- cs$pairs
  rep_k <- c(seq_along(pr$w), sample(seq_along(pr$w), 500))
  expect_same_sums(cs$X, pr$i[rep_k], pr$j[rep_k], pr$w[rep_k])
  w <- pr$w
  w[sample(length(w), 1000)] <- 0
  expect_same_sums(cs$X, pr$i, pr$j, w)
  expect_same_sums(cs$X, pr$i, pr$j, pmin(pr$w, 0.5) - 0.5)
  expect_equal(.pair_sums(cs$X, pr$i, pr$j, 0 * pr$w), matrix(0, nrow(cs$X), 1))
})


test_that("no pairs or no samples give zeros", {
  cs <- make_case()
  pr <- cs$pairs
  expect_equal(.pair_sums(cs$X, integer(), integer(), numeric()), matrix(0, nrow(cs$X), 1))
  expect_equal(
    .pair_sums(cs$X, integer(), integer(), numeric(), col_group = cs$unit_of, n_groups = cs$n_units),
    matrix(0, nrow(cs$X), cs$n_units)
  )
  expect_equal(.pair_sums(cs$X[0, ], pr$i, pr$j, pr$w), matrix(0, 0, 1))
})


test_that("works with dense, indicator and full rows, empty samples and absent subunits", {
  cs <- make_case()
  pr <- cs$pairs
  wi <- cs$within
  X <- cs$X
  X[5, ] <- 0
  X[, 1:40] <- 0
  expect_same_sums(X, pr$i, pr$j, pr$w)
  expect_same_sums(as.matrix(X), pr$i, pr$j, pr$w)
  expect_same_sums((X > 0) * 1, pr$i, pr$j, pr$w)
  # The reference row of relative.multiplicity has no zeros
  expect_same_sums(rbind(X, 1), wi$i, wi$j, wi$w, unit_of = cs$unit_of, n_units = cs$n_units)
})


test_that("pairs sharing the first or the second subunit give the same sums", {
  cs <- make_case()
  sp <- sorted_unit_pairs(cs$unit_of)
  # Few distinct first subunits (as in the sorted files), then few distinct second ones
  first <- sp$i <= 25
  expect_lt(length(unique(sp$i[first])), length(unique(sp$j[first])))
  for (k in list(list(sp$i, sp$j), list(sp$j, sp$i))) {
    i <- k[[1]][first]
    j <- k[[2]][first]
    w <- sp$w[first]
    expect_same_sums(cs$X, i, j, w)
    expect_same_sums(cs$X, i, j, w, Y = cs$Y)
    expect_same_sums(cs$X, i, j, w, unit_of = cs$unit_of, n_units = cs$n_units)
  }
  expect_same_sums(cs$X, sp$i, sp$j, sp$w, unit_of = cs$unit_of, n_units = cs$n_units)
})


test_that("results do not depend on the size of the internal products", {
  cs <- make_case()
  pr <- cs$pairs
  wi <- cs$within
  sp <- sorted_unit_pairs(cs$unit_of)
  for (sizes in list(c(1, 1), c(13, 3), c(500, 1e6), c(1e9, 7))) {
    local_mocked_bindings(.chunk_work = sizes[1], .chunk_pairs = sizes[2])
    expect_same_sums(cs$X, pr$i, pr$j, pr$w, Y = cs$Y)
    expect_same_sums(cs$X, wi$i, wi$j, wi$w, unit_of = cs$unit_of, n_units = cs$n_units)
    expect_same_sums(cs$X, sp$i, sp$j, sp$w, unit_of = cs$unit_of, n_units = cs$n_units)
  }
})


test_that("grouped sums reject pairs of different groups", {
  cs <- make_case()
  pr <- cs$pairs
  expect_error(
    .pair_sums(cs$X, pr$i, pr$j, pr$w, col_group = cs$unit_of, n_groups = cs$n_units),
    "pairs of subunits of different groups"
  )
})
