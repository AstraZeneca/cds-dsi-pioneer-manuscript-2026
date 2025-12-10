# Enhanced parquet_draws using Arrow streaming + targets branching
# nolint start: object_usage_linter

#' Convert single CSV to parquet using Arrow's fast streaming reader
#'
#' This is much faster than manual readLines() approach because:
#' - Arrow's C++ CSV parser is 10-100x faster
#' - Handles type inference automatically
#' - Uses optimized memory layouts
#' - Can filter columns before materializing
#'
#' IMPORTANT: Stan CSVs have comments (#) interspersed throughout:
#' - Header comments at start
#' - "# Adaptation terminated" etc. AFTER the column names
#' - Footer comments at end (timing info)
#' Arrow's read_csv_arrow() can't handle interspersed comments, so we
#' filter them out first.
#'
#' @param csv_file Path to CmdStan CSV output
#' @param output_parquet Path for output parquet file
#' @param chain_id Chain identifier (integer)
#' @param variables_regex Optional vector of regex patterns to filter variables
#' @return Path to created parquet file
csv_to_parquet_arrow_single <- function(csv_file, output_parquet, chain_id, variables_regex = NULL) {
  # Stan CSVs have comments INTERSPERSED in data section!
  # Arrow can't handle this, so filter out ALL comment lines first
  lines <- readLines(csv_file, warn = FALSE)
  data_lines <- lines[!startsWith(lines, "#")]

  # Write clean CSV to temp file for Arrow to parse
  temp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_csv), add = TRUE)
  writeLines(data_lines, temp_csv)

  # Arrow's streaming CSV reader - this is the key optimization!
  # It reads in optimized chunks using fast C++ code
  tbl <- arrow::read_csv_arrow(
    temp_csv, # Clean CSV without any comments
    as_data_frame = FALSE, # Keep as Arrow Table for zero-copy operations
    read_options = arrow::CsvReadOptions(
      use_threads = TRUE, # Parallel CSV parsing
      block_size = 1048576 # 1MB chunks
    )
  )

  # Filter columns if patterns provided - happens in Arrow (fast!)
  if (!is.null(variables_regex)) {
    pattern <- stringr::str_c("^(", stringr::str_c(variables_regex, collapse = "|"), ")")

    all_cols <- names(tbl)

    # Separate parameter columns from Stan diagnostics
    diagnostic_cols <- c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")
    param_cols <- setdiff(all_cols, diagnostic_cols)

    # Select matching parameters plus key diagnostics
    matched_params <- param_cols[stringr::str_detect(param_cols, pattern)]

    if (length(matched_params) == 0) {
      stop("No variables matched pattern: ", pattern, " in file: ", csv_file)
    }

    # Always keep lp__ and accept_stat__ for diagnostics
    cols_to_keep <- unique(c(matched_params, "lp__", "accept_stat__"))

    # Arrow select() is zero-copy - just creates a view
    tbl <- tbl |> dplyr::select(dplyr::any_of(cols_to_keep))
  }

  # Add chain metadata - still in Arrow land (fast!)
  tbl <- tbl |>
    dplyr::mutate(
      .chain = chain_id,
      .iteration = dplyr::row_number(),
      .draw = .iteration # Posterior package convention
    )

  # Write to parquet - Arrow handles compression efficiently
  arrow::write_parquet(
    tbl,
    output_parquet,
    compression = "zstd", # Fast compression, good ratio
    compression_level = 3 # Balance speed vs size
  )

  invisible(output_parquet)
}

#' Create targets for parallel CSV-to-Parquet conversion using branching
#'
#' This creates separate targets for each chain conversion, allowing targets
#' to parallelize them automatically using its built-in crew workers.
#'
#' Benefits over future::future_map:
#' - Targets tracks each chain separately (better caching)
#' - Failed chains don't invalidate successful ones
#' - Built-in retry logic
#' - Better integration with targets dependency graph
#' - Uses existing crew workers (no separate future plan needed)
#'
#' @param name_prefix Prefix for target names (e.g., "my_fit_parquet")
#' @param fit_target Name of the target containing CmdStanFit object
#' @param variables_regex Optional vector of regex patterns to filter variables
#' @param ... Additional arguments passed to tar_target
#' @return List of targets: one per chain + final combiner
#'
#' @examples
#' # In your _targets.R:
#' list(
#'   tar_target(my_fit, model$sample(...)),
#'   tar_parquet_draws_branching("my_draws", my_fit)
#' )
tar_parquet_draws_branching <- function(name_prefix, fit_target, variables_regex = NULL, ...) {
  name_prefix <- as.character(substitute(name_prefix))
  fit_sym <- substitute(fit_target)

  list(
    # Step 1: Extract CSV file paths from fit
    targets::tar_target_raw(
      name = paste0(name_prefix, "_csv_info"),
      command = rlang::expr({
        fit <- !!fit_sym
        csv_files <- fit$output_files()

        tibble::tibble(
          chain_id = seq_along(csv_files),
          csv_file = csv_files,
          parquet_file = file.path(
            targets::tar_path_target(!!name_prefix),
            stringr::str_c("chain_", chain_id, ".parquet")
          )
        )
      })
    ),

    # Step 2: Create output directory
    targets::tar_target_raw(
      name = paste0(name_prefix, "_dir"),
      command = rlang::expr({
        output_dir <- targets::tar_path_target(!!name_prefix)
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        output_dir
      }),
      format = "file"
    ),

    # Step 3: Branch over chains - THIS IS WHERE PARALLEL MAGIC HAPPENS!
    # Each chain becomes a separate target that can run in parallel
    targets::tar_target_raw(
      name = paste0(name_prefix, "_chains"),
      command = rlang::expr({
        csv_info <- !!as.name(paste0(name_prefix, "_csv_info"))
        output_dir <- !!as.name(paste0(name_prefix, "_dir"))

        # This will be called once per chain in parallel!
        csv_to_parquet_arrow_single(
          csv_file = csv_info$csv_file,
          output_parquet = csv_info$parquet_file,
          chain_id = csv_info$chain_id,
          variables_regex = !!variables_regex
        )
      }),
      pattern = rlang::expr(map(!!as.name(paste0(name_prefix, "_csv_info")))),
      format = "file",
      ...
    ),

    # Step 4: Combine into parquet_draws object
    targets::tar_target_raw(
      name = name_prefix,
      command = rlang::expr({
        # Wait for all chains to complete
        chain_files <- !!as.name(paste0(name_prefix, "_chains"))
        output_dir <- !!as.name(paste0(name_prefix, "_dir"))

        # Return parquet_draws object pointing to directory
        new_parquet_draws(output_dir)
      }),
      format = tar_format(
        read = function(path) new_parquet_draws(path),
        write = function(object, path) object$path
      )
    )
  )
}

#' Simpler version that returns a single combined target
#'
#' This is more convenient when you just want the draws object and don't
#' need to see intermediate chain targets.
#'
#' @param name Target name for the final parquet_draws object
#' @param fit Expression that evaluates to CmdStanFit
#' @param variables_regex Optional variable filter
#' @param ... Passed to tar_target
#' @return A targets target
tar_parquet_draws_simple <- function(name, fit, variables_regex = NULL, ...) {
  name <- targets::tar_deparse_language(substitute(name))
  fit_expr <- substitute(fit)

  targets::tar_target_raw(
    name = name,
    command = rlang::expr({
      fit <- !!fit_expr
      csv_files <- fit$output_files()
      final_path <- targets::tar_path_target()

      dir.create(final_path, showWarnings = FALSE, recursive = TRUE)

      # Use Arrow streaming for each chain sequentially
      # (Targets will parallelize multiple tar_parquet_draws_simple calls)
      purrr::walk(seq_along(csv_files), \(i) {
        parquet_path <- file.path(final_path, stringr::str_c("chain_", i, ".parquet"))
        csv_to_parquet_arrow_single(
          csv_files[i],
          parquet_path,
          i,
          !!variables_regex
        )
        cat(sprintf("Chain %d: %s -> %s (%.1f MB)\n", i, basename(csv_files[i]), basename(parquet_path), file.size(parquet_path) / 1024^2))
      })

      new_parquet_draws(final_path)
    }),
    format = tar_format(
      read = function(path) new_parquet_draws(path),
      write = function(object, path) object$path
    ),
    ...
  )
}

# nolint end: object_usage_linter
