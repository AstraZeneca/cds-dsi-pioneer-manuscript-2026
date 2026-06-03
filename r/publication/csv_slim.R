#' Project a CmdStan CSV to a narrow column subset
#'
#' Streams the input file line by line. Comment lines (starting with `#`) pass
#' through verbatim. The first non-comment line is the header; column indices
#' for `keep_col_names` are resolved once from it. Every subsequent non-comment
#' line (a draw row) is projected to those column indices.
#'
#' @param input_csv  Path to the source CmdStan CSV (can be very large).
#' @param output_csv Path for the projected output CSV (created/overwritten).
#' @param keep_col_names Character vector of column base-names to retain.
#'   These must match the header field names exactly (no index suffixes).
#'   Include the 7 sampler diagnostics (lp__, accept_stat__, ..., energy__)
#'   explicitly — they are not auto-retained.
#' @return `output_csv` invisibly.
project_cmdstan_csv <- function(input_csv, output_csv, keep_col_names) {
  # Validate keep_col_names before creating output file
  keep_idx <- .resolve_keep_idx(input_csv, keep_col_names)

  con_in  <- file(input_csv,  open = "r", encoding = "UTF-8")
  con_out <- file(output_csv, open = "w", encoding = "UTF-8")
  on.exit({ close(con_in); close(con_out) }, add = TRUE)

  repeat {
    line <- readLines(con_in, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (startsWith(line, "#")) {
      writeLines(line, con_out)
      next
    }
    writeLines(paste(strsplit(line, ",", fixed = TRUE)[[1]][keep_idx], collapse = ","), con_out)
  }

  invisible(output_csv)
}

.resolve_keep_idx <- function(input_csv, keep_col_names) {
  con <- file(input_csv, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("No header line found in ", input_csv)
    if (!startsWith(line, "#")) break
  }
  fields <- strsplit(line, ",", fixed = TRUE)[[1]]
  missing_cols <- setdiff(keep_col_names, fields)
  if (length(missing_cols) > 0L) {
    stop("keep_col_names columns not found in CSV header: ",
         paste(missing_cols, collapse = ", "))
  }
  match(keep_col_names, fields)
}
