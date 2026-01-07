# Enhanced parquet_draws using Arrow streaming + targets branching
# nolint start: object_usage_linter

#' Convert Stan CSV files to parquet using DuckDB's optimized CSV reader
#'
#' This is much faster than Arrow approach because:
#' - DuckDB's CSV parser natively handles comments
#' - Can read multiple CSV files in parallel (one per chain)
#' - Column selection happens during read (not after)
#' - Direct CSV-to-Parquet streaming (no intermediate steps)
#' - No memory overhead from loading unwanted columns
#'
#' IMPORTANT: Stan CSVs have comments (#) interspersed throughout:
#' - Header comments at start
#' - "# Adaptation terminated" etc. AFTER the column names
#' - Footer comments at end (timing info)
#' DuckDB handles this natively with the comment='#' parameter.
#'
#' @param csv_files Vector of paths to CmdStan CSV outputs (one per chain)
#' @param output_parquet Path for output parquet file
#' @param ... Optional dplyr-style selection expressions to filter variables
#' @return Path to created parquet file
csv_to_parquet_duckdb <- function(csv_files, output_parquet, ..., max_threads = Inf, max_memory = "50GB") {
  # Get column selection expressions
  dots <- rlang::enquos(...)
  n_chains <- length(csv_files)

  # Stan diagnostic columns to always preserve
  diagnostic_cols <- c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")

  # Connect to DuckDB with explicit config for memory management
  # CRITICAL: Config must be passed to duckdb() driver, not dbConnect()
  drv <- duckdb::duckdb(
    config = list(
      "memory_limit" = max_memory,
      "temp_directory" = "/tmp/duckdb_temp",
      "threads" = if (is.finite(max_threads)) as.character(max_threads)
    ) |>
      purrr::compact()
  )
  con <- withr::local_db_connection(DBI::dbConnect(drv))

  # Also set preserve_insertion_order after connection
  DBI::dbExecute(con, "SET preserve_insertion_order = false")
  DBI::dbExecute(con, "SET enable_progress_bar = true")
  DBI::dbExecute(con, "SET enable_progress_bar_print = true")

  # Enable profiling
  # DBI::dbExecute(con, "SET enable_profiling = 'JSON'")
  # DBI::dbExecute(con, "SET profiling_mode = 'DETAILED'")

  if (length(dots) == 0) {
    select_clause <- "*"
    types_clause <- "types = {'*': 'DOUBLE'},"
  } else {
    # Read just the header to get column names for tidyselect
    header <- readr::read_lines(csv_files[1], n_max = 100) |>
      stringr::str_subset("^#", negate = TRUE) |>
      head(1) |>
      stringr::str_split_1(",")

    # Create simplified names for tidyselect matching (removes array indices)
    # This allows using bare symbols like `states` instead of matches("^states")
    # e.g., "states.1.1" -> "states", "beta.2.3" -> "beta"
    simplified_names <- stringr::str_replace(header, "\\.[0-9.]+$", "")

    # Apply tidyselect to get selected column indices
    # Use simplified names for matching, but preserve original column names
    selected_idx <- tidyselect::eval_select(
      rlang::expr(c(dplyr::any_of(!!diagnostic_cols), !!!dots)),
      data = rlang::set_names(header, simplified_names)
    )

    # Get the original column names (with array indices) for the selected columns
    selected_cols <- header[selected_idx]
    names(selected_cols) <- selected_cols # Set names for compatibility

    # Quote column names to handle dots and other special characters
    # DuckDB requires double quotes for identifiers with special chars
    quoted_cols <- stringr::str_c('"', names(selected_cols), '"')
    select_clause <- paste(quoted_cols, collapse = ", ")

    # Build explicit types map for selected columns
    types_list <- paste0("'", names(selected_cols), "': 'DOUBLE'", collapse = ", ")
    types_clause <- stringr::str_glue("types = {{{types_list}}},")
  }

  # Build file list for DuckDB (needs to be SQL array syntax)
  file_list <- paste(shQuote(csv_files, type = "sh"), collapse = ", ")

  # Profile CSV reading operation
  # DBI::dbExecute(con, "SET profiling_output = 'profile_read.json'")

  query_read <- stringr::str_glue(
    "
    CREATE OR REPLACE TABLE temp_data AS 
    SELECT 
      CAST(regexp_extract(filename, '-(\\d+)-[a-f0-9]+\\.csv', 1) AS INTEGER) as chain_id,
      {select_clause}
    FROM read_csv([{file_list}], 
                  filename = true,
                  comment = '#',
                  header = true,
                  max_line_size = 104857600,
                  auto_detect = false,
                  all_varchar = false,
                  {types_clause}
                  union_by_name = true,
                  parallel = true)
    "
  )

  # browser()

  DBI::dbExecute(con, query_read)

  # Profile Parquet writing operation
  # DBI::dbExecute(con, "SET profiling_output = 'profile_write.json'")

  query_write <- stringr::str_glue(
    "
    COPY temp_data TO '{output_parquet}' 
    (FORMAT PARQUET, COMPRESSION UNCOMPRESSED)
    "
  )

  DBI::dbExecute(con, query_write)

  # Clean up temporary table
  DBI::dbExecute(con, "DROP TABLE temp_data")

  return(output_parquet)
}

#' Targets format for parquet draws with proper hashing
#'
#' This format handles directories containing parquet files, computing a hash
#' of all files to detect changes. On read, it returns an Arrow dataset.
#'
#' @return A tar_format object for use with targets
#' @export
parquet_draws_format <- function() {
  targets::tar_format(
    read = function(path) {
      if (inherits(path, "parquet_draws")) {
        arrow::open_dataset(path$path)
      } else {
        arrow::open_dataset(path)
      }
    },
    write = function(object, path) {
      # Hash all files if object is a directory
      if (fs::is_dir(object)) {
        files <- fs::dir_ls(object, recurse = TRUE, type = "file")
        hash_val <- digest::digest(purrr::map(files, \(f) digest::digest(f, file = TRUE)))
      } else {
        hash_val <- digest::digest(object, file = TRUE)
      }

      structure(
        list(
          path = object,
          hash = hash_val
        ),
        class = "parquet_draws"
      )
    }
  )
}

# S3 methods for Arrow datasets to work with tidybayes/posterior ================

#' Convert Arrow Dataset to posterior draws_df
#'
#' This allows tidybayes::spread_rvars and gather_rvars to work
#' transparently with Arrow-backed draws.
#'
#' @param x An Arrow Dataset containing MCMC draws
#' @param ... Additional arguments (unused)
#' @return A draws_df object
#' @export
as_draws_df.Dataset <- function(x, ...) {
  df <- dplyr::collect(x)

  # Convert Stan CSV column names (dot notation) to bracket notation
  # e.g., "beta.1.2" -> "beta[1,2]"
  names(df) <- cmdstanr:::repair_variable_names(names(df))

  # Extract chain_id if present, otherwise assume single chain
  if ("chain_id" %in% names(df)) {
    chain_ids <- df$chain_id
    df$chain_id <- NULL
  } else {
    chain_ids <- rep(1L, nrow(df))
  }

  # Create iteration indices within each chain
  df$.chain <- chain_ids
  df$.iteration <- ave(seq_len(nrow(df)), chain_ids, FUN = seq_along)
  df$.draw <- seq_len(nrow(df))

  # Reorder to have metadata columns first
  df <- dplyr::select(df, .chain, .iteration, .draw, dplyr::everything())

  posterior::as_draws_df(df)
}

#' Convert Arrow Dataset to posterior draws_rvars
#'
#' @param x An Arrow Dataset containing MCMC draws
#' @param ... Additional arguments (unused)
#' @return A draws_rvars object
#' @export
as_draws_rvars.Dataset <- function(x, ...) {
  posterior::as_draws_rvars(as_draws_df.Dataset(x, ...))
}

#' Convert Arrow Dataset to posterior draws
#'
#' @param x An Arrow Dataset containing MCMC draws
#' @param ... Additional arguments (unused)
#' @return A draws_df object
#' @export
as_draws.Dataset <- function(x, ...) {
  as_draws_df.Dataset(x, ...)
}

#' Register S3 methods for Arrow datasets
#'
#' Call this function to register the as_draws methods for Arrow datasets.
#' This is called automatically when the package is loaded.
#' @export
register_arrow_draws_methods <- function() {
  if (requireNamespace("posterior", quietly = TRUE)) {
    registerS3method("as_draws", "Dataset", as_draws.Dataset, envir = asNamespace("posterior"))
    registerS3method("as_draws_df", "Dataset", as_draws_df.Dataset, envir = asNamespace("posterior"))
    registerS3method("as_draws_rvars", "Dataset", as_draws_rvars.Dataset, envir = asNamespace("posterior"))
  }
}

# Auto-register methods when this file is sourced
register_arrow_draws_methods()

# Column-selective spread_rvars and gather_rvars for Arrow Datasets ==========
#
# These methods solve the memory problem where as_draws_df.Dataset would load
# the entire parquet file into memory. Instead, we:
# 1. Parse the variable specifications to determine which columns are needed
# 2. Use Arrow's column projection to select only those columns
# 3. Collect the filtered data (much smaller)
# 4. Then call the standard tidybayes methods

#' Extract base variable names from tidybayes variable specifications
#'
#' Given quosures like `spop_target_pfs[i]` or `states[n, p]`, extracts
#' the base variable name (e.g., "spop_target_pfs", "states").
#'
#' @param quos List of quosures from enquos(...)
#' @return Character vector of base variable names
#' @keywords internal
extract_variable_names <- function(quos) {

  purrr::map_chr(quos, function(q) {
    spec <- tidybayes:::parse_variable_spec(q)
    spec[[1]] # First element is always the variable name
  })
}

#' Build column selection regex for Arrow dataset
#'
#' Given base variable names, builds regex patterns that match the Stan CSV
#' column format (e.g., "var.1.2" for arrays).
#'
#' @param var_names Character vector of base variable names
#' @return Character regex pattern
#' @keywords internal
build_column_regex <- function(var_names) {

  # Match either exact name or name followed by dot and indices

  # e.g., "states" matches "states" and "states.1.2"
  patterns <- paste0("^", var_names, "(\\.[0-9.]+)?$")
  paste(patterns, collapse = "|")
}

#' Select columns from Arrow Dataset matching variable specs
#'
#' @param dataset An Arrow Dataset
#' @param var_names Character vector of base variable names
#' @return Arrow Dataset with only matching columns (plus metadata columns)
#' @keywords internal
select_matching_columns <- function(dataset, var_names) {
  all_cols <- names(dataset)


  # Always keep chain_id and diagnostic columns
  keep_cols <- c("chain_id", "lp__", "accept_stat__", "stepsize__",
                 "treedepth__", "n_leapfrog__", "divergent__", "energy__")
  keep_cols <- intersect(keep_cols, all_cols)


  # Build regex and find matching columns

  pattern <- build_column_regex(var_names)
  matching_cols <- all_cols[grepl(pattern, all_cols)]

  # Combine and select
  cols_to_select <- unique(c(keep_cols, matching_cols))
  dplyr::select(dataset, dplyr::all_of(cols_to_select))
}

#' spread_rvars method for Arrow Datasets with lazy column selection
#'
#' This method intercepts spread_rvars calls on Arrow Datasets and only loads
#' the columns that are actually needed, dramatically reducing memory usage.
#'
#' @param model An Arrow Dataset containing MCMC draws
#' @param ... Variable specifications (e.g., `beta[i]`, `sigma`)
#' @param ndraws Optional number of draws to subsample
#' @return A tibble with rvars
#' @export
spread_rvars.Dataset <- function(model, ..., ndraws = NULL) {
  quos <- rlang::enquos(...)

  # Extract variable names from specs

  var_names <- extract_variable_names(quos)

  # Select only needed columns from Arrow dataset
  filtered_dataset <- select_matching_columns(model, var_names)

  # Collect the filtered data (now much smaller)
  df <- dplyr::collect(filtered_dataset)

  # Convert column names from Stan CSV format to bracket notation
  names(df) <- cmdstanr:::repair_variable_names(names(df))

  # Build draws_df
  if ("chain_id" %in% names(df)) {
    chain_ids <- df$chain_id
    df$chain_id <- NULL
  } else {
    chain_ids <- rep(1L, nrow(df))
  }

  df$.chain <- chain_ids
  df$.iteration <- ave(seq_len(nrow(df)), chain_ids, FUN = seq_along)
  df$.draw <- seq_len(nrow(df))
  df <- dplyr::select(df, .chain, .iteration, .draw, dplyr::everything())

  draws <- posterior::as_draws_df(df)

  # Preserve tidybayes_constructors if present (for recover_types)
  constructors <- attr(model, "tidybayes_constructors")
  if (!is.null(constructors)) {
    attr(draws, "tidybayes_constructors") <- constructors
  }

  # Call tidybayes spread_rvars on the draws_df

  tidybayes::spread_rvars(draws, !!!quos, ndraws = ndraws)
}

#' gather_rvars method for Arrow Datasets with lazy column selection
#'
#' This method intercepts gather_rvars calls on Arrow Datasets and only loads
#' the columns that are actually needed, dramatically reducing memory usage.
#'
#' @param model An Arrow Dataset containing MCMC draws
#' @param ... Variable specifications (e.g., `beta[i]`, `sigma`)
#' @param ndraws Optional number of draws to subsample
#' @return A tibble with rvars in long format
#' @export
gather_rvars.Dataset <- function(model, ..., ndraws = NULL) {
  quos <- rlang::enquos(...)

  # Extract variable names from specs
  var_names <- extract_variable_names(quos)

  # Select only needed columns from Arrow dataset
  filtered_dataset <- select_matching_columns(model, var_names)

  # Collect the filtered data (now much smaller)
  df <- dplyr::collect(filtered_dataset)

  # Convert column names from Stan CSV format to bracket notation
  names(df) <- cmdstanr:::repair_variable_names(names(df))

  # Build draws_df
  if ("chain_id" %in% names(df)) {
    chain_ids <- df$chain_id
    df$chain_id <- NULL
  } else {
    chain_ids <- rep(1L, nrow(df))
  }

  df$.chain <- chain_ids
  df$.iteration <- ave(seq_len(nrow(df)), chain_ids, FUN = seq_along)
  df$.draw <- seq_len(nrow(df))
  df <- dplyr::select(df, .chain, .iteration, .draw, dplyr::everything())

  draws <- posterior::as_draws_df(df)

  # Preserve tidybayes_constructors if present (for recover_types)
  constructors <- attr(model, "tidybayes_constructors")
  if (!is.null(constructors)) {
    attr(draws, "tidybayes_constructors") <- constructors
  }

  # Call tidybayes gather_rvars on the draws_df
  tidybayes::gather_rvars(draws, !!!quos, ndraws = ndraws)
}

#' Decorate Arrow Dataset with type recovery information
#'
#' This allows using recover_types() with Arrow datasets so that
#' spread_rvars and gather_rvars can convert indices back to factors.
#'
#' @param model An Arrow Dataset containing MCMC draws
#' @param ... Data frames or lists providing type prototypes
#' @return The dataset with tidybayes_constructors attribute
#' @export
recover_types.Dataset <- function(model, ...) {
  # Use tidybayes's make_constructors to build the constructor list
  constructors <- tidybayes::make_constructors(...)
  attr(model, "tidybayes_constructors") <- constructors
  model
}

# Register recover_types methods
register_recover_types_methods <- function() {
  if (requireNamespace("tidybayes", quietly = TRUE)) {
    registerS3method("recover_types", "Dataset", recover_types.Dataset, envir = asNamespace("tidybayes"))
  }
}

register_recover_types_methods()

# Register spread_rvars and gather_rvars methods for Arrow Datasets
register_rvars_methods <- function() {
  if (requireNamespace("tidybayes", quietly = TRUE)) {
    registerS3method("spread_rvars", "Dataset", spread_rvars.Dataset, envir = asNamespace("tidybayes"))
    registerS3method("gather_rvars", "Dataset", gather_rvars.Dataset, envir = asNamespace("tidybayes"))
  }
}

register_rvars_methods()

# nolint end: object_usage_linter
