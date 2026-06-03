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

# Column base-names kept in the CRC slim CSV — union of all four
# select_draws_single_chain() calls in publication_targets.R.
# Extend this vector if new draws targets are added that read new columns.
CRC_SLIM_KEEP_COLS <- c(
  # sampler diagnostics (always required by read_cmdstan_csv)
  "lp__", "accept_stat__", "stepsize__", "treedepth__",
  "n_leapfrog__", "divergent__", "energy__",
  # tumor_ssls_draws_pop_chain
  "measure_sd_sld",
  # tumor_ssls_draws_sld_recist_chain
  "rep_patient_log_sld", "forecast_patient_log_sld",
  "rep_recist", "forecast_obs_recist",
  # tumor_ssls_draws_endpoints_chain
  "recist_confusion_matrix"
)

#' Build the full keep-column list for a CRC fit CSV
#'
#' Reads the header of `input_csv`, matches all publication selection patterns,
#' and returns the unique set of column names to retain.
#'
#' @param input_csv Path to one CRC fit CSV (used only to read the header).
#' @return Character vector of column names (with index suffixes) to retain.
build_crc_keep_col_names <- function(input_csv) {
  con <- file(input_csv, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("No header line found in ", input_csv)
    if (!startsWith(line, "#")) break
  }
  all_cols <- strsplit(line, ",", fixed = TRUE)[[1]]

  # Strip trailing dotted numeric index suffix to get base names.
  base_names <- sub("\\.[0-9]+(\\.[0-9]+)*$", "", all_cols)

  # Patterns from publication_targets.R select_draws_single_chain() calls:
  patterns <- c(
    "_pop$",
    "^pop_",
    "_sd_level_",
    "^(time_invariant|time_varying)_coef",
    "^(frac|init|tr)_.+_patient",
    "patient_log_(growth|decrease)_rate",
    "(spop|sample)(_target|_ms)?_(((quant_)?(pfs|os))|km_est|right_censored|(pfs|os)_n)",
    "(spop|sample)_target_(((un)?confirmed_response)|orr)",
    "(spop|sample)_(os|pfs)_(quant|km_est|n)",
    "(spop|sample)_os(_censored)?",
    "(spop|sample)_(os|pfs)_quant_exceeds_max"
  )

  matched <- unique(c(
    all_cols[base_names %in% CRC_SLIM_KEEP_COLS],
    unlist(lapply(patterns, function(p) all_cols[grepl(p, base_names)]))
  ))
  matched
}

#' Project one CmdStan CSV to a slim copy, for use in a targets pipeline
#'
#' The output is written to `<output_dir>/<original_filename>`.
#'
#' @param input_csv    Path to source CmdStan CSV (one chain).
#' @param output_dir   Directory for the projected output (created if needed).
#' @param keep_col_names Character vector of column names from `build_crc_keep_col_names()`.
#' @return Path to the slim output CSV.
slim_cmdstan_csv_files <- function(input_csv, output_dir, keep_col_names) {
  fs::dir_create(output_dir, recurse = TRUE)
  output_csv <- file.path(output_dir, basename(input_csv))
  project_cmdstan_csv(input_csv, output_csv, keep_col_names)
  output_csv
}
