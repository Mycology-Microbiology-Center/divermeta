#' Check a large distance file
#'
#' Checks that a distance file lists every needed pair of subunits exactly once, with a valid
#' distance, without loading the file into memory. The indices run this check when called with
#' `check_large_distance_file = TRUE`; it can also be run on its own, e.g. once, before computing
#' several indices from the same file.
#'
#' @param file Path to the distance file, or several paths read one after the other. See
#'   [distance-files] for the format.
#' @param ids Character vector with the identifiers of the subunits, e.g. `colnames(abund)`. Rows
#'   with other identifiers are ignored.
#' @param clust Vector of cluster memberships of the subunits: named by subunit identifier, or
#'   following the order of `ids`. A data frame is also accepted, with two columns (subunit
#'   identifiers, units) or one column of units. Only needed for `scope` "within" or "between".
#' @param scope Which pairs of subunits must be listed: "all" (default) every pair, as needed by
#'   [raoQuadratic()], [diversity.functional()], [diversity.functional.traditional()],
#'   [redundancy()] and [multiplicity.distance()] with a linkage method; "within" the pairs of
#'   subunits of the same unit, as needed by [multiplicity.distance()] with `method = "sigma"`,
#'   [relative.multiplicity()], [relative.multiplicity.ref_div()] and [average.redundancy()];
#'   "between" the pairs of subunits of different units, as needed by [unit_distances()].
#'   [divermeta()] needs "all" when any of the requested indices does, and "within" otherwise.
#'   Rows of pairs outside the scope are ignored.
#' @param chunk_size Number of rows read at a time (default `1e6`).
#' @param max_pairs_in_memory Maximum number of pairs expected in memory at a time (default
#'   `1e7`, about 160 MB). Lower values use less memory and more temporary files.
#' @param tmp_dir Directory for the temporary files (default [tempdir()]).
#' @param header Logical (default `TRUE`). Whether every file starts with a header line. See
#'   [distance-files].
#'
#' @return A list of class `divermeta_distance_check`, printed as a summary, with:
#' \describe{
#'   \item{ok}{`TRUE` if every needed pair is listed exactly once with a valid distance.}
#'   \item{scope}{The scope checked.}
#'   \item{n_rows}{Number of rows read.}
#'   \item{n_expected}{Number of pairs needed.}
#'   \item{n_used}{Number of valid rows of pairs in the scope.}
#'   \item{n_ignored}{Number of rows ignored: unknown identifiers, self pairs and pairs outside
#'     the scope.}
#'   \item{n_invalid}{Number of rows with a non-numeric distance, and of rows in the scope with an
#'     NA, infinite or negative distance. NA, infinite and negative distances outside the scope are
#'     ignored, as the indices do.}
#'   \item{n_missing}{Number of needed pairs with no row at all (a row with an invalid distance
#'     lists its pair, and is counted in `n_invalid`).}
#'   \item{n_duplicated}{Number of rows with a valid distance listing a pair already listed with a
#'     valid distance (in either orientation).}
#'   \item{n_conflicting}{Number of those rows with a different distance.}
#'   \item{examples}{An example of every problem found.}
#' }
#'
#' @details
#' The file is read twice, never as a whole:
#' \enumerate{
#'   \item It is read in chunks of `chunk_size` rows. Every row of a needed pair is written to a
#'     temporary binary file on disk (a bucket), chosen by its first subunit so that every bucket
#'     holds about `max_pairs_in_memory` pairs.
#'   \item Every bucket is loaded in turn to find the duplicated and conflicting pairs and to
#'     count the pairs of every subunit, which gives the missing pairs.
#' }
#' The temporary files take about 16 bytes of disk per pair (e.g. about 80 GB for
#' \eqn{5 \cdot 10^9}{5e9} pairs) and are deleted at the end.
#'
#' @seealso [distance-files]
#'
#' @export
#'
#' @examples
#' diss <- data.frame(
#'   ID1 = c("a", "a", "b", "b"),
#'   ID2 = c("b", "c", "c", "a"),
#'   Distance = c(0.4, 0.9, 0.6, 0.5)
#' )
#' diss_file <- tempfile(fileext = ".tsv")
#' write.table(diss, diss_file, sep = "\t", quote = FALSE, row.names = FALSE)
#'
#' # The pair a-b is listed twice, with different distances
#' check_distance_file(diss_file, ids = c("a", "b", "c"))
#'
check_distance_file <- function(
  file,
  ids,
  clust = NULL,
  scope = "all",
  chunk_size = 1e6,
  max_pairs_in_memory = 1e7,
  tmp_dir = tempdir(),
  header = TRUE
) {
  if (!.is_diss_file(file)) {
    stop("`file` must be the path to a distance file")
  }
  if (length(ids) == 0 || anyNA(ids)) {
    stop("`ids` must be a non-empty vector of subunit identifiers")
  }
  ids <- .as_ids(ids)
  if (anyDuplicated(ids) > 0) {
    stop("`ids` (subunit identifiers) must be unique")
  }
  if (!is.character(scope) || length(scope) != 1 ||
    !(scope %in% c("all", "within", "between"))) {
    stop("`scope` must be one of \"all\", \"within\" or \"between\"")
  }
  cl <- NULL
  if (scope != "all") {
    if (is.null(clust)) {
      stop(paste0("`clust` is needed for scope \"", scope, "\""))
    }
    cl <- .parse_clust(clust, ids)
  }
  .check_chunk_size(chunk_size)
  if (!is.numeric(max_pairs_in_memory) || length(max_pairs_in_memory) != 1 ||
    !is.finite(max_pairs_in_memory) || max_pairs_in_memory < 1) {
    stop("`max_pairs_in_memory` must be a single positive number")
  }
  if (!is.character(tmp_dir) || length(tmp_dir) != 1 || !dir.exists(tmp_dir)) {
    stop("`tmp_dir` must be an existing directory")
  }
  .check_flag(header, "header")
  .check_diss_files(file, "file", header)

  .check_distance_file(file, ids, cl, scope, chunk_size, max_pairs_in_memory, tmp_dir, header = header)
}


.check_distance_file <- function(
  paths,
  ids,
  cl,
  scope,
  chunk_size,
  max_pairs_in_memory = 1e7,
  tmp_dir = tempdir(),
  name = "diss",
  header = TRUE
) {
  n <- length(ids)
  nn <- as.numeric(n)
  unit_of <- cl$unit_of

  # Buckets: ranges of the first subunit i, with about max_pairs_in_memory
  # expected pairs each
  expected <- .expected_partners(scope, n, cl)
  bucket_of <- pmax(1L, as.integer(ceiling(cumsum(expected) / max_pairs_in_memory)))
  n_buckets <- if (n > 0) max(bucket_of) else 0L

  dir <- tempfile("divermeta_check_", tmpdir = tmp_dir)
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  bucket_file <- function(b, what) file.path(dir, paste0("bucket_", b, "_", what, ".bin"))
  append_bin <- function(x, path) {
    con <- file(path, "ab")
    on.exit(close(con))
    writeBin(x, con)
  }

  n_rows <- 0
  n_used <- 0
  n_invalid <- 0
  ex_invalid <- NULL

  # Pass 1: read the file and write the rows of the needed pairs to the buckets
  for (path in paths) {
    con <- .open_diss_file(path)
    file_rows <- 0
    tryCatch(
      {
        sep <- .read_diss_header(con, path, header)$sep
        repeat {
          x <- .scan_diss_chunk(con, sep, chunk_size, path, distance = "")
          if (is.null(x)) {
            break
          }
          raw <- x[[3]]
          d <- suppressWarnings(as.numeric(raw))

          i1 <- match(.as_ids(x[[1]]), ids)
          i2 <- match(.as_ids(x[[2]]), ids)
          i <- pmin(i1, i2)
          j <- pmax(i1, i2)
          used <- !is.na(i) & !is.na(j) & i != j
          if (scope == "within") {
            used[used] <- unit_of[i[used]] == unit_of[j[used]]
          } else if (scope == "between") {
            used[used] <- unit_of[i[used]] != unit_of[j[used]]
          }

          non_numeric <- is.na(d) & !is.na(raw) & raw != ""
          invalid <- non_numeric | (used & !(is.finite(d) & d >= 0))
          if (any(invalid)) {
            k <- which(invalid)[1]
            if (is.null(ex_invalid)) {
              ex_invalid <- paste0(
                x[[1]][k], "-", x[[2]][k], " (distance: ", raw[k], ", row ", file_rows + k,
                " of ", basename(path), ")"
              )
            }
            n_invalid <- n_invalid + sum(invalid)
          }
          n_rows <- n_rows + length(d)
          file_rows <- file_rows + length(d)

          n_used <- n_used + sum(used & !invalid)
          # A row of a needed pair with an invalid distance still lists the pair
          # (its distance is kept as NA), so the pair is not also counted as missing
          if (!any(used)) {
            next
          }
          d[invalid] <- NA_real_
          i <- i[used]
          key <- (i - 1) * nn + (j[used] - 1)
          d <- d[used]
          rows <- split(seq_along(i), bucket_of[i])
          for (b in names(rows)) {
            r <- rows[[b]]
            append_bin(key[r], bucket_file(b, "key"))
            append_bin(d[r], bucket_file(b, "d"))
          }
        }
      },
      finally = close(con)
    )
  }

  # Pass 2: every bucket in turn
  n_missing <- 0
  n_duplicated <- 0
  n_conflicting <- 0
  ex_missing <- NULL
  ex_duplicated <- NULL
  ex_conflicting <- NULL
  pair_name <- function(key) {
    paste0(ids[key %/% nn + 1], "-", ids[key %% nn + 1])
  }
  read_bin <- function(path) {
    if (!file.exists(path)) {
      return(numeric(0))
    }
    readBin(path, "double", n = file.size(path) / 8)
  }

  for (b in seq_len(n_buckets)) {
    rows_b <- which(bucket_of == b)
    key <- read_bin(bucket_file(b, "key"))
    d <- read_bin(bucket_file(b, "d"))

    # Rows with an invalid distance (NA) list their pair, so it is not missing,
    # but they are counted in n_invalid only: duplicates are among valid rows
    valid <- !is.na(d)
    key_valid <- key[valid]
    d_valid <- d[valid]
    dup <- duplicated(key_valid)
    if (any(dup)) {
      n_duplicated <- n_duplicated + sum(dup)
      if (is.null(ex_duplicated)) {
        ex_duplicated <- pair_name(key_valid[which(dup)[1]])
      }
      first <- match(key_valid, key_valid)
      conflict <- dup & abs(d_valid - d_valid[first]) > sqrt(.Machine$double.eps)
      if (any(conflict)) {
        n_conflicting <- n_conflicting + sum(conflict)
        if (is.null(ex_conflicting)) {
          k <- which(conflict)[1]
          ex_conflicting <- paste0(
            pair_name(key_valid[k]), " (", d_valid[first[k]], " and ", d_valid[k], ")"
          )
        }
      }
    }
    key <- unique(key)

    # Partners listed for every subunit of the bucket
    i <- key %/% nn + 1
    counts <- tabulate(i - rows_b[1] + 1, length(rows_b))
    short <- counts < expected[rows_b]
    if (any(short)) {
      n_missing <- n_missing + sum(expected[rows_b] - counts)
      if (is.null(ex_missing)) {
        r <- rows_b[which(short)[1]]
        partners <- if (r < n) (r + 1):n else integer(0)
        if (scope == "within") {
          partners <- partners[unit_of[partners] == unit_of[r]]
        } else if (scope == "between") {
          partners <- partners[unit_of[partners] != unit_of[r]]
        }
        listed <- key[i == r] %% nn + 1
        ex_missing <- paste0(ids[r], "-", ids[setdiff(partners, listed)[1]])
      }
    }
  }

  structure(
    list(
      ok = n_missing == 0 && n_duplicated == 0 && n_invalid == 0,
      scope = scope,
      n_rows = n_rows,
      n_expected = sum(expected),
      n_used = n_used,
      n_ignored = n_rows - n_used - n_invalid,
      n_invalid = n_invalid,
      n_missing = n_missing,
      n_duplicated = n_duplicated,
      n_conflicting = n_conflicting,
      examples = list(
        invalid = ex_invalid,
        missing = ex_missing,
        duplicated = ex_duplicated,
        conflicting = ex_conflicting
      )
    ),
    class = "divermeta_distance_check"
  )
}


#' @export
format.divermeta_distance_check <- function(x, ...) {
  scope <- switch(x$scope,
    all = "every pair of subunits",
    within = "pairs of subunits of the same unit",
    between = "pairs of subunits of different units"
  )
  line <- function(label, n, example = NULL) {
    paste0(
      "  ", label, ": ", format(n, big.mark = ",", scientific = FALSE),
      if (!is.null(example)) paste0(" (e.g. ", example, ")") else ""
    )
  }
  c(
    paste0("Distance file check (", scope, "): ", if (x$ok) "OK" else "FAILED"),
    line("Rows read", x$n_rows),
    line("Pairs needed", x$n_expected),
    line("Rows used", x$n_used),
    line("Rows ignored", x$n_ignored),
    line("Invalid distances", x$n_invalid, x$examples$invalid),
    line("Missing pairs", x$n_missing, x$examples$missing),
    line("Duplicated pairs", x$n_duplicated, x$examples$duplicated),
    line("Conflicting duplicates", x$n_conflicting, x$examples$conflicting)
  )
}


#' @export
print.divermeta_distance_check <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
