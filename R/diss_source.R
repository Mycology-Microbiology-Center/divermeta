#' Distances from files
#'
#' Every index that uses distances accepts `diss` either as a three-column table in memory or as
#' the path to a file (or several files, read one after the other) with that table. A file is read
#' in chunks of `chunk_size` rows, so the table never has to fit in memory: memory depends on the
#' number of samples, subunits and units, and on `chunk_size`, not on the number of pairs.
#'
#' @section File format:
#' \itemize{
#'   \item Delimited text with at least three columns, taken by position: the identifiers of two
#'     subunits and the distance between them. Further columns are ignored.
#'   \item Files ending in `.csv` are comma separated and files ending in `.tsv` are tab
#'     separated. For other extensions the separator (tab, comma or white space) is detected
#'     from the first line.
#'   \item With `header = TRUE` (default) the first line of every file is a header and is
#'     skipped; if its third field is a number, the line looks like data and the index stops, so a
#'     row is never dropped silently. Files without a header need `header = FALSE`, and then a
#'     first row with a non-numeric distance is an error like any other row.
#'   \item Identifiers are read as text and matched as text against the column names of `abund`
#'     (so `100000` in a file does not match a column named `1e+05`).
#'   \item Files may be plain or compressed with gzip (`.gz`), bzip2 (`.bz2`) or xz (`.xz`).
#'     Files compressed with zstd (`.zst`) are read through the `zstd` command line tool, which
#'     must be installed.
#' }
#'
#' @section Every pair exactly once:
#' Pairs are not kept in memory, so they cannot be de-duplicated while reading: a file must list
#' every pair of subunits the index needs **exactly once** (in either orientation). Rows with
#' identifiers not in `abund`, self pairs and pairs the index does not use (e.g. pairs of subunits
#' of different units for indices that only use the pairs inside units) are ignored, also when
#' their distance is NA, infinite or negative. The exception is [metagenomic.alpha.index()],
#' which checks the distance of every row with two known identifiers: there, an NA, infinite or
#' negative distance is an error even in a pair it does not use.
#'
#' Pairs with a subunit that has zero abundance in every sample are not needed (except for the
#' homogeneous reference of [relative.multiplicity()], also when requested through
#' [divermeta()]): if some of them are missing, a warning is
#' given and the index is computed without them. When they are listed, they count as any other
#' row. When `diss` is a file:
#' \itemize{
#'   \item The file is always checked to exist, to be readable and to have three columns, and
#'     every chunk is checked for NA, infinite or negative distances in the pairs the index
#'     uses. Every distance must still be readable as a number.
#'   \item The number of rows of the needed pairs is compared with the number of pairs: too few
#'     or too many rows is an error.
#'   \item With `check_large_distance_file = FALSE` (default), nothing else is checked and a
#'     warning says so. A missing pair would count as distance 0 (as if the two subunits were
#'     identical) and a duplicated pair twice, and a missing pair and a duplicated one cancel out
#'     in the row count.
#'   \item With `check_large_distance_file = TRUE`, [check_distance_file()] is run first, which
#'     finds every missing, duplicated, conflicting or invalid row, and the index stops if there
#'     is any. The check reads the file once more and needs about 16 bytes of temporary disk
#'     space per pair.
#' }
#' In-memory tables are always fully checked, and there a pair listed twice with the same
#' distance is accepted.
#'
#' @seealso [check_distance_file()]
#'
#' @examples
#' abund <- matrix(c(2, 1, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("S1", "S2"), c("a", "b", "c")))
#' diss <- data.frame(ID1 = c("a", "a", "b"), ID2 = c("b", "c", "c"), Distance = c(0.4, 0.9, 0.6))
#'
#' diss_file <- tempfile(fileext = ".csv")
#' write.csv(diss, diss_file, row.names = FALSE)
#'
#' raoQuadratic(abund, diss_file, check_large_distance_file = TRUE)
#'
#' # A file without a header
#' write.table(diss, diss_file, sep = ",", row.names = FALSE, col.names = FALSE)
#' raoQuadratic(abund, diss_file, check_large_distance_file = TRUE, header = FALSE)
#'
#' @name distance-files
NULL


# Pair sources
# ----------------------------
# The indices read the distances through a pair source: a list with
#   - next_chunk(): the next pairs as list(i, j, d), with i < j indices into the
#     subunit identifiers and d their distances, or NULL when there are no more
#   - close(): releases the connections
# An in-memory table is a single chunk, already validated and de-duplicated,
# and is also kept as `pairs` for exact checks. A file is read chunk by chunk.


.is_diss_file <- function(diss) {
  is.character(diss) && is.null(dim(diss))
}


.check_chunk_size <- function(chunk_size) {
  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || !is.finite(chunk_size) ||
    chunk_size < 1) {
    stop("`chunk_size` must be a single positive number")
  }
}


# Pair source of `diss` (table or file paths) indexed by `ids`. `scope` is the
# set of pairs the caller needs ("all", "within" or "between" units, or "none"
# when the caller checks its pairs itself): rows outside it are dropped before
# their distances are checked, and a file is fully checked for it when `check`
# is TRUE, and a warning is given otherwise. `header` says whether every file
# starts with a header line. `need` (logical, one per subunit; default all TRUE)
# marks the subunits whose pairs must be listed: pairs with a subunit not needed
# (zero abundance in every sample) are used when listed, but may be missing
.diss_source <- function(
  diss,
  ids,
  chunk_size = 1e6,
  check = FALSE,
  scope = "all",
  cl = NULL,
  name = "diss",
  header = TRUE,
  need = NULL
) {
  .check_chunk_size(chunk_size)
  .check_flag(check, "check_large_distance_file")
  if (is.null(need)) {
    need <- rep(TRUE, length(ids))
  }
  needed <- .needed_subunits(ids, cl, need)

  if (!.is_diss_file(diss)) {
    pairs <- .parse_diss(diss, ids, name, scope, cl)
    done <- FALSE
    return(list(
      next_chunk = function() {
        if (done) {
          return(NULL)
        }
        done <<- TRUE
        pairs
      },
      close = function() invisible(NULL),
      pairs = pairs,
      need = need
    ))
  }

  .check_flag(header, "header")
  .check_diss_files(diss, name, header)
  if (scope != "none") {
    if (check) {
      report <- .check_distance_file(
        diss, needed$ids, needed$cl, scope, chunk_size, name = name, header = header
      )
      if (!report$ok) {
        stop(paste(c("The distance file did not pass the check:", format(report)), collapse = "\n"))
      }
    } else {
      warning(paste0(
        "The distance file was not checked (`check_large_distance_file = FALSE`): every needed ",
        "pair of subunits is assumed to be listed exactly once. Only the number of rows is ",
        "compared with the number of pairs; a missing pair would count as distance 0 and a ",
        "duplicated one twice. Use `check_large_distance_file = TRUE` or `check_distance_file()` ",
        "to verify the file."
      ), call. = FALSE)
    }
  }

  pairs_needed <- scope != "none" &&
    sum(.expected_partners(scope, length(needed$ids), needed$cl)) > 0
  src <- .diss_file_source(diss, ids, chunk_size, name, scope, cl, header, pairs_needed)
  src$need <- need
  src
}


# The subunits whose pairs must be listed (`need`), with the clustering
# restricted to them. Units are kept, also when none of their subunits is needed
.needed_subunits <- function(ids, cl, need) {
  cl_need <- if (is.null(cl)) NULL else list(unit_of = cl$unit_of[need], unit_ids = cl$unit_ids)
  list(ids = ids[need], cl = cl_need)
}


# Pair source reading the files one after the other, `chunk_size` rows at a time.
# Rows outside `scope` are dropped (see .index_pairs). With `pairs_needed`, a
# file with no row pairing two subunits is an error (its first columns are
# likely not subunit identifiers)
.diss_file_source <- function(
  paths,
  ids,
  chunk_size,
  name = "diss",
  scope = "all",
  cl = NULL,
  header = TRUE,
  pairs_needed = TRUE
) {
  k <- 0L
  con <- NULL
  sep <- NULL
  path <- NULL
  n_rows <- 0
  n_kept <- 0

  close_current <- function() {
    if (!is.null(con)) {
      close(con)
      con <<- NULL
    }
    invisible(NULL)
  }

  next_chunk <- function() {
    repeat {
      if (is.null(con)) {
        if (k >= length(paths)) {
          if (pairs_needed && n_rows > 0 && n_kept == 0) {
            stop(paste0(
              "No row of the distance file pairs two subunits of `abund`: the first two ",
              "columns must hold subunit identifiers (the column names of `abund`)"
            ))
          }
          return(NULL)
        }
        k <<- k + 1L
        path <<- paths[k]
        con <<- .open_diss_file(path)
        sep <<- .read_diss_header(con, path, header)$sep
      }

      x <- .scan_diss_chunk(con, sep, chunk_size, path)
      if (is.null(x)) {
        close_current()
        next
      }
      pairs <- .index_pairs(x[[1]], x[[2]], x[[3]], ids, name, scope, cl)
      n_rows <<- n_rows + length(x[[1]])
      n_kept <<- n_kept + pairs$n_known
      return(pairs)
    }
  }

  list(next_chunk = next_chunk, close = close_current)
}


# Opens a plain or compressed text file for reading
.open_diss_file <- function(path) {
  if (grepl("\\.zst$", path, ignore.case = TRUE)) {
    if (!nzchar(Sys.which("zstd"))) {
      stop(paste0(
        "Reading `.zst` files requires the `zstd` command line tool, which was not found: ",
        path
      ))
    }
    return(pipe(paste("zstd -dc", shQuote(path)), "r"))
  }
  # gzfile() reads plain, gzip, bzip2 and xz files
  gzfile(path, "r")
}


# Field separator: from the extension for .csv and .tsv files (compressed or
# not), otherwise from the first line
.diss_sep <- function(path, line) {
  base <- sub("\\.(gz|bz2|xz|zst)$", "", path, ignore.case = TRUE)
  if (grepl("\\.csv$", base, ignore.case = TRUE)) {
    return(",")
  }
  if (grepl("\\.tsv$", base, ignore.case = TRUE)) {
    return("\t")
  }
  if (grepl("\t", line, fixed = TRUE)) {
    return("\t")
  }
  if (grepl(",", line, fixed = TRUE)) {
    return(",")
  }
  ""
}


# Reads the first line of a file: finds the separator and checks there are
# three columns. With `header`, the line is skipped, and it must not look like
# data (a number as third field); otherwise it is put back. Returns the separator
.read_diss_header <- function(con, path, header = TRUE) {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (length(line) == 0) {
    return(list(sep = "\t"))
  }

  sep <- .diss_sep(path, line)
  fields <- if (sep == "") {
    strsplit(trimws(line), "[[:space:]]+")[[1]]
  } else {
    trimws(strsplit(line, sep, fixed = TRUE)[[1]])
  }
  fields <- gsub("^\"|\"$", "", fields)
  if (length(fields) < 3) {
    stop(paste0(
      "The distance file must have three columns (id_1, id_2, distance), separated by ",
      "tabs, commas or spaces: ", path
    ))
  }

  if (!header) {
    pushBack(line, con)
  } else if (!is.na(suppressWarnings(as.numeric(fields[3])))) {
    stop(paste0(
      "The first line of the distance file looks like data, not a header (its third field is ",
      "the number ", fields[3], "). Use `header = FALSE` for files without a header: ", path
    ))
  }
  list(sep = sep)
}


# Next `chunk_size` rows of a file (the first three columns), or NULL at the end.
# The distances are read as numbers, or as text with `distance = ""`
.scan_diss_chunk <- function(con, sep, chunk_size, path, distance = 0) {
  x <- tryCatch(
    scan(
      con,
      what = list("", "", distance),
      nmax = chunk_size,
      sep = sep,
      quote = "\"",
      strip.white = TRUE,
      flush = TRUE,
      multi.line = FALSE,
      na.strings = "NA",
      quiet = TRUE
    ),
    error = function(e) {
      stop(paste0(
        "Could not read the distance file ", path, ": ", conditionMessage(e),
        ". It must have three columns (id_1, id_2, distance) with numeric distances"
      ), call. = FALSE)
    }
  )
  if (length(x[[1]]) == 0) {
    return(NULL)
  }
  x
}


# Checks that the files exist, can be read and have three columns (and a
# header, if expected)
.check_diss_files <- function(paths, name = "diss", header = TRUE) {
  if (length(paths) == 0 || anyNA(paths)) {
    stop(paste0("`", name, "` must be a data frame or the path to a distance file"))
  }
  missing <- paths[!file.exists(paths) | dir.exists(paths)]
  if (length(missing) > 0) {
    stop(paste0("Distance file not found: ", paste(missing, collapse = ", ")))
  }
  unreadable <- paths[file.access(paths, 4) != 0]
  if (length(unreadable) > 0) {
    stop(paste0("Distance file cannot be read: ", paste(unreadable, collapse = ", ")))
  }
  for (path in paths) {
    con <- .open_diss_file(path)
    tryCatch(.read_diss_header(con, path, header), finally = close(con))
  }
  invisible(NULL)
}


# Reading the pairs
# ----------------------------

# Reads the pair source once and feeds every chunk to every accumulator, a list
# with `scope` (the pairs it needs: "all", "within", "between" or "none"),
# `update(pairs)` and `finish()`. Checks that the needed pairs are listed: exactly
# for in-memory tables, by the number of rows for files. Pairs with a subunit
# that is not needed (`src$need`) may be missing, with a warning. Returns the
# result of `finish()` of every accumulator
.consume_pairs <- function(src, accs, ids, cl = NULL, name = "diss") {
  on.exit(src$close(), add = TRUE)
  scopes <- intersect(c("all", "within", "between"), vapply(accs, `[[`, "", "scope"))
  need <- if (is.null(src$need)) rep(TRUE, length(ids)) else src$need
  in_memory <- !is.null(src$pairs)
  if (in_memory) {
    .check_pairs_in_memory(src$pairs, scopes, ids, cl, name, need)
  }

  n_units <- if (is.null(cl)) 0L else length(cl$unit_ids)
  tally <- list(all = .new_tally(n_units), need = .new_tally(n_units))
  repeat {
    pairs <- src$next_chunk()
    if (is.null(pairs)) {
      break
    }
    if (length(pairs$i) == 0) {
      next
    }
    if (!in_memory) {
      nd <- need[pairs$i] & need[pairs$j]
      tally$all <- .add_tally(tally$all, pairs$i, pairs$j, cl)
      tally$need <- .add_tally(tally$need, pairs$i[nd], pairs$j[nd], cl)
    }
    for (acc in accs) {
      acc$update(pairs)
    }
  }

  if (!in_memory) {
    needed <- .needed_subunits(ids, cl, need)
    .check_pair_counts(tally$need, scopes, needed$ids, needed$cl, name)
    .check_absent_pairs(tally$all, tally$need, scopes, ids, cl, need, name)
  }
  lapply(accs, function(acc) acc$finish())
}


# Number of pairs read: in total, inside every unit and between units
.new_tally <- function(n_units) {
  list(rows = 0, within = numeric(n_units), between = 0)
}


.add_tally <- function(tally, i, j, cl) {
  tally$rows <- tally$rows + length(i)
  if (!is.null(cl) && length(i) > 0) {
    u <- cl$unit_of[i]
    within <- u == cl$unit_of[j]
    tally$within <- tally$within + tabulate(u[within], length(cl$unit_ids))
    tally$between <- tally$between + sum(!within)
  }
  tally
}


# Number of pairs of `scope` listed in a tally
.tally_scope <- function(tally, scope) {
  switch(scope,
    all = tally$rows,
    within = sum(tally$within),
    between = tally$between
  )
}


# Compares the pairs with a subunit that is not needed (zero abundance in every
# sample) with their number: missing ones give a warning, and more rows than
# pairs (only possible in files) are an error. `listed_all` and `listed_need` are
# tallies of all pairs and of the pairs of needed subunits
.check_absent_pairs <- function(listed_all, listed_need, scopes, ids, cl, need, name = "diss") {
  if (all(need)) {
    return(invisible(NULL))
  }
  needed <- .needed_subunits(ids, cl, need)
  missing <- FALSE
  for (scope in scopes) {
    expected <- sum(.expected_partners(scope, length(ids), cl)) -
      sum(.expected_partners(scope, length(needed$ids), needed$cl))
    listed <- .tally_scope(listed_all, scope) - .tally_scope(listed_need, scope)
    if (listed > expected) {
      stop(paste0(
        "The distance file lists ", listed, " rows for the ", expected, " pairs with subunits ",
        "that have zero abundance in every sample: some pairs are listed more than once. ",
        "Use `check_distance_file()` to find the pairs"
      ))
    }
    missing <- missing || listed < expected
  }
  if (missing) {
    absent <- ids[!need]
    warning(paste0(
      "`", name, "` does not list every pair of the subunits with zero abundance in every sample (",
      paste(utils::head(absent, 5), collapse = ", "), if (length(absent) > 5) ", ..." else "",
      "). These pairs are not needed, and the missing ones are ignored."
    ), call. = FALSE)
  }
  invisible(NULL)
}


# Stops unless every needed pair of an in-memory table is listed. Only the pairs
# of needed subunits (`need`) are required; missing pairs with other subunits
# give a warning
.check_pairs_in_memory <- function(pairs, scopes, ids, cl, name = "diss", need = NULL) {
  if (is.null(need)) {
    need <- rep(TRUE, length(ids))
  }
  nd <- need[pairs$i] & need[pairs$j]
  needed <- .needed_subunits(ids, cl, need)
  pos <- which(need)
  pairs_need <- list(i = match(pairs$i[nd], pos), j = match(pairs$j[nd], pos), d = pairs$d[nd])

  .require_scope_pairs(pairs_need, scopes, needed$ids, needed$cl, name)

  if (!all(need)) {
    n_units <- if (is.null(cl)) 0L else length(cl$unit_ids)
    listed_all <- .add_tally(.new_tally(n_units), pairs$i, pairs$j, cl)
    listed_need <- .add_tally(.new_tally(n_units), pairs$i[nd], pairs$j[nd], cl)
    .check_absent_pairs(listed_all, listed_need, scopes, ids, cl, need, name)
  }
  invisible(NULL)
}


# Stops unless every pair of `scopes` is listed in the (de-duplicated) pairs
.require_scope_pairs <- function(pairs, scopes, ids, cl, name = "diss") {
  if ("all" %in% scopes) {
    .require_all_pairs(pairs, ids, name)
  }
  if ("within" %in% scopes) {
    .within_unit_pairs(pairs, cl, ids, name)
  }
  if ("between" %in% scopes) {
    expected <- .n_between_pairs(cl)
    listed <- sum(cl$unit_of[pairs$i] != cl$unit_of[pairs$j])
    if (listed < expected) {
      stop(paste0(
        "Missing distances in `", name, "` for ", expected - listed,
        " pair(s) of subunits of different units"
      ))
    }
  }
  invisible(NULL)
}


.n_between_pairs <- function(cl) {
  sizes <- tabulate(cl$unit_of, length(cl$unit_ids))
  (sum(sizes)^2 - sum(sizes^2)) / 2
}


# Stops unless the number of rows read from a file matches the number of
# needed pairs
.check_pair_counts <- function(tally, scopes, ids, cl, name = "diss") {
  hint <- paste0(
    "A distance file must list every pair exactly once: use `check_distance_file()` ",
    "to find the pairs"
  )
  count_error <- function(listed, expected, what) {
    if (listed < expected) {
      stop(paste0(
        "Missing distances in `", name, "` for ", expected - listed, " pair(s) of ", what,
        ": the distance file lists ", listed, " of the ", expected, " pairs. ", hint
      ))
    }
    if (listed > expected) {
      stop(paste0(
        "The distance file lists ", listed, " rows for the ", expected, " pairs of ", what,
        ": some pairs are listed more than once. ", hint
      ))
    }
  }

  if ("all" %in% scopes) {
    n <- length(ids)
    count_error(tally$rows, n * (n - 1) / 2, "subunits")
  }
  if ("within" %in% scopes) {
    sizes <- tabulate(cl$unit_of, length(cl$unit_ids))
    expected <- sizes * (sizes - 1) / 2
    wrong <- which(tally$within != expected)
    if (length(wrong) > 0) {
      u <- wrong[1]
      count_error(tally$within[u], expected[u], paste0("subunits in unit: ", cl$unit_ids[u]))
    }
  }
  if ("between" %in% scopes) {
    count_error(tally$between, .n_between_pairs(cl), "subunits of different units")
  }
  invisible(NULL)
}


# Number of partners j > i of every subunit i in the pairs of `scope`
.expected_partners <- function(scope, n, cl = NULL) {
  after <- n - seq_len(n)
  if (scope == "all") {
    return(after)
  }
  u <- cl$unit_of
  sizes <- tabulate(u, length(cl$unit_ids))
  o <- order(u)
  rank <- integer(n)
  rank[o] <- seq_len(n) - (cumsum(sizes) - sizes)[u[o]]
  same_after <- sizes[u] - rank
  if (scope == "within") same_after else after - same_after
}
