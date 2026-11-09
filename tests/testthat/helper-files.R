# Distance files
# ----------------------------

# Writes a three-column distance table to a temporary file, deleted when the
# calling test ends. `format` sets the extension and the layout:
#   "csv", "tsv", "csv.gz", "tsv.bz2", "tsv.xz": with a header
#   "txt": white space separated, without a header
#   "tab": tab separated, no extension, without a header
# (see has_header). Rows are shuffled and some orientations flipped
write_diss <- function(diss, format = "csv", shuffle = TRUE, env = parent.frame()) {
  diss <- as.data.frame(diss, stringsAsFactors = FALSE)
  names(diss)[1:3] <- c("ID1", "ID2", "Distance")
  if (shuffle && nrow(diss) > 1) {
    diss <- diss[sample(nrow(diss)), ]
    flip <- runif(nrow(diss)) < 0.5
    diss[flip, c("ID1", "ID2")] <- diss[flip, c("ID2", "ID1")]
  }

  ext <- switch(format, tab = "", paste0(".", format))
  path <- withr::local_tempfile(fileext = ext, .local_envir = env)
  sep <- switch(sub("\\..*$", "", format), csv = ",", txt = " ", "\t")
  header <- has_header(format)

  con <- if (grepl("\\.gz$", format)) {
    gzfile(path, "w")
  } else if (grepl("\\.bz2$", format)) {
    bzfile(path, "w")
  } else if (grepl("\\.xz$", format)) {
    xzfile(path, "w")
  } else {
    file(path, "w")
  }
  utils::write.table(
    diss, con,
    sep = sep, quote = sep == ",", row.names = FALSE, col.names = header
  )
  close(con)
  path
}

file_formats <- c("csv", "tsv", "csv.gz", "tsv.bz2", "tsv.xz", "txt", "tab")

# Whether the files written by write_diss() in `format` have a header
has_header <- function(format) !(format %in% c("txt", "tab"))


# Random study: samples x subunits abundances (with an empty sample, absent
# subunits and an absent unit in some samples), a clustering in four units and
# the complete distance table
make_file_study <- function(n = 24, n_samples = 6, seed = 1) {
  set.seed(seed)
  ids <- paste0("f", seq_len(n))
  clust <- stats::setNames(sample(rep(c("A", "B", "C", "D"), length.out = n)), ids)
  abund <- matrix(
    runif(n_samples * n, 1, 50) * (runif(n_samples * n) > 0.3),
    nrow = n_samples,
    dimnames = list(paste0("S", seq_len(n_samples)), ids)
  )
  abund[n_samples, ] <- 0
  abund[1, clust == "A"] <- 0
  list(
    abund = abund,
    clust = clust,
    diss = mat_to_table(rand_diss(n, 0.05, 1.2), ids)
  )
}


# Evaluates an index with a file expecting the "not checked" warning
unchecked <- function(expr) {
  expect_warning(res <- expr, "not checked")
  res
}
