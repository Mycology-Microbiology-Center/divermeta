# Legacy implementations
# ----------------------------
# The single sample, matrix based implementations of the indices, kept in
# tests/testthat/legacy/ as a reference for the tests. They are sourced into
# their own environment, so they call each other and not the package versions:
# use them as legacy$multiplicity.distance(...), legacy$raoQuadratic(...), etc.
legacy <- new.env(parent = globalenv())
for (f in sort(list.files(test_path("legacy"), pattern = "\\.R$", full.names = TRUE))) {
  sys.source(f, envir = legacy)
}


# Shared helpers
# ----------------------------

# Upper triangle of a square matrix as a three-column distance table
mat_to_table <- function(m, ids = NULL) {
  m <- as.matrix(m)
  if (is.null(ids)) {
    ids <- if (is.null(rownames(m))) as.character(seq_len(nrow(m))) else rownames(m)
  }
  idx <- which(upper.tri(m), arr.ind = TRUE)
  data.frame(
    ID1 = ids[idx[, 1]],
    ID2 = ids[idx[, 2]],
    Distance = m[idx],
    stringsAsFactors = FALSE
  )
}

# Single sample abundance table (one row, subunits in columns)
one_row <- function(ab, ids = NULL, sample = "S1") {
  if (is.null(ids)) {
    ids <- as.character(seq_along(ab))
  }
  matrix(ab, nrow = 1, dimnames = list(sample, ids))
}

# Random symmetric distance matrix with zero diagonal
rand_diss <- function(n, min = 0.1, max = 0.9) {
  m <- matrix(runif(n * n, min = min, max = max), nrow = n, ncol = n)
  m <- (m + t(m)) / 2
  diag(m) <- 0
  m
}

# Matrix with every pair at distance sig
max_diss <- function(n, sig = 1) {
  m <- matrix(sig, nrow = n, ncol = n)
  diag(m) <- 0
  m
}


# Block diagonal matrix: blocks on the diagonal, `off` everywhere else
block_matrix <- function(blocks, off) {
  sizes <- vapply(blocks, nrow, integer(1))
  diss <- matrix(off, nrow = sum(sizes), ncol = sum(sizes))
  offset <- 0
  for (M in blocks) {
    idx <- offset + seq_len(nrow(M))
    diss[idx, idx] <- M
    offset <- offset + nrow(M)
  }
  diag(diss) <- 0
  diss
}
