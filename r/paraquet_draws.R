# nolint start: object_usage_linter

# Define the parquet_draws class
new_parquet_draws <- function(path) {
  structure(
    list(path = path),
    class = c("parquet_draws", "draws")
  )
}

# {targets} format and functions ####

# Custom format for parquet_draws - use with regular draws objects
tar_format_parquet_draws <- tar_format(
  read = function(path) {
    new_parquet_draws(path)
  },

  write = function(object, path) {
    arrow::write_parquet(object, path)
    path
  }
)

# Memory-efficient CSV-to-Parquet converter for CmdStan fits
tar_parquet_draws <- function(name, fit, ..., variables_regex = NULL) {
  name <- targets::tar_deparse_language(substitute(name))
  fit_expr <- substitute(fit)

  targets::tar_target_raw(
    name = name,
    command = rlang::expr({
      fit <- !!fit_expr
      csv_files <- fit$output_files()
      metadata <- fit$metadata()
      iter_per_chain <- metadata[["iter_sampling"]]
      final_path <- targets::tar_path_target()

      # Create output directory for partitioned parquet files
      dir.create(final_path, showWarnings = FALSE, recursive = TRUE)

      # Function to stream CSV directly to parquet file
      csv_to_parquet_file <- function(input_csv, output_parquet, chain_id, patterns = NULL, batch_size = 10) {
        con_in <- withr::local_connection(file(input_csv, "r"))

        # Read and filter header (first non-comment line)
        header_line <- NULL
        repeat {
          line <- readLines(con_in, n = 1)
          if (length(line) == 0) {
            break
          } # EOF
          if (!stringr::str_starts(line, "#")) {
            header_line <- line
            break
          }
        }

        if (is.null(header_line)) {
          stop("No header found in ", input_csv)
        }

        # Parse column names and find matching columns
        col_names <- stringr::str_split_1(header_line, ",")

        if (!is.null(patterns)) {
          pattern <- stringr::str_c("^(", stringr::str_c(patterns, collapse = "|"), ")_")
          col_indices <- stringr::str_which(col_names, pattern)
          if (length(col_indices) == 0) {
            stop("No columns matched pattern: ", pattern)
          }
        } else {
          col_indices <- seq_along(col_names)
        }

        selected_col_names <- col_names[col_indices]

        # Collect all data rows in batches
        all_batches <- list()
        batch_rows <- list()

        repeat {
          line <- readLines(con_in, n = 1)
          if (length(line) == 0) {
            break
          } # EOF

          if (!stringr::str_starts(line, "#") && nchar(line) > 0) {
            fields <- stringr::str_split_1(line, ",")
            batch_rows[[length(batch_rows) + 1]] <- as.numeric(fields[col_indices])

            # When batch is full, convert to data frame and store
            if (length(batch_rows) >= batch_size) {
              df <- as.data.frame(do.call(rbind, batch_rows))
              names(df) <- selected_col_names
              all_batches[[length(all_batches) + 1]] <- df
              batch_rows <- list()
            }
          }
        }

        # Process remaining rows
        if (length(batch_rows) > 0) {
          df <- as.data.frame(do.call(rbind, batch_rows))
          names(df) <- selected_col_names
          all_batches[[length(all_batches) + 1]] <- df
        }

        # Combine all batches and add metadata
        if (length(all_batches) > 0) {
          full_df <- do.call(rbind, all_batches)
          full_df$.chain <- chain_id
          full_df$.draw <- seq_len(nrow(full_df))

          # Write once to parquet
          arrow::write_parquet(full_df, output_parquet)
        }

        output_parquet
      }

      # Process each CSV to separate parquet file
      purrr::walk(seq_along(csv_files), \(i) {
        parquet_path <- file.path(final_path, stringr::str_c("chain_", i, ".parquet"))
        csv_to_parquet_file(csv_files[i], parquet_path, i, !!variables_regex)

        cat("Converted ", csv_files[i], " to ", parquet_path, "\n")
      })

      # Return parquet_draws object
      new_parquet_draws(final_path)
    }),
    format = tar_format(
      read = function(path) {
        new_parquet_draws(path)
      },
      write = function(object, path) {
        # Files already written by command to object$path
        # Just return the path for targets to track
        object$path
      }
    ),
    ...
  )
}

# Utility and overridden functions

.get_draws <- function(data, ...) {
  dots <- rlang::enquos(...)
  ds <- arrow::open_dataset(data$path)

  if (length(dots) > 0) {
    all_cols <- names(ds)
    param_names <- get_param_names_from_dots(dots, all_cols)

    ds |>
      dplyr::select(dplyr::all_of(param_names), .chain, .iteration, .draw) |>
      dplyr::collect() |>
      posterior::as_draws_df()
  } else {
    posterior::as_draws_df(data)
  }
}

# Add method for spread_rvars (already a generic in tidybayes)
spread_rvars.parquet_draws <- function(data, ...) {
  .get_draws(data, ...) |>
    tidybayes::spread_rvars(...)
}

# Add method for gather_rvars (already a generic in tidybayes)
gather_rvars.parquet_draws <- function(data, ...) {
  .get_draws(data, ...) |>
    tidybayes::gather_rvars(...)
}

# posterior package methods ####

# Subset draws - this is where performance wins happen
subset_draws.parquet_draws <- function(x, variable = NULL, draw = NULL, chain = NULL, iteration = NULL, regex = FALSE, ...) {
  ds <- arrow::open_dataset(x$path)

  # Start with metadata columns
  cols_to_read <- c(".chain", ".iteration", ".draw")

  # Add variable selection
  if (!is.null(variable)) {
    all_cols <- names(ds)
    param_cols <- setdiff(all_cols, c(".chain", ".iteration", ".draw"))

    if (regex) {
      # Treat variable as regex patterns
      matching <- param_cols[grepl(paste(variable, collapse = "|"), param_cols)]
    } else {
      # Handle Stan bracket notation: "beta" should match "beta[1]", "beta[2]", etc.
      patterns <- paste0("^", variable, "(\\[|$)")
      matching <- param_cols[grepl(paste(patterns, collapse = "|"), param_cols)]
    }

    cols_to_read <- c(cols_to_read, matching)
  }

  # Read selected columns
  result <- ds |>
    dplyr::select(dplyr::all_of(cols_to_read)) |>
    dplyr::collect()

  # Apply filters if specified
  if (!is.null(draw)) {
    result <- result |> dplyr::filter(.draw %in% draw)
  }
  if (!is.null(chain)) {
    result <- result |> dplyr::filter(.chain %in% chain)
  }
  if (!is.null(iteration)) {
    result <- result |> dplyr::filter(.iteration %in% iteration)
  }

  posterior::as_draws_df(result)
}

# Metadata methods - avoid reading full dataset
variables.parquet_draws <- function(x, reserved = FALSE, ...) {
  ds <- arrow::open_dataset(x$path)
  all_vars <- names(ds)

  if (reserved) {
    all_vars
  } else {
    setdiff(all_vars, c(".chain", ".iteration", ".draw"))
  }
}

nchains.parquet_draws <- function(x, ...) {
  arrow::open_dataset(x$path) |>
    dplyr::select(.chain) |>
    dplyr::summarise(n = max(.chain)) |>
    dplyr::collect() |>
    dplyr::pull(n)
}

niterations.parquet_draws <- function(x, ...) {
  arrow::open_dataset(x$path) |>
    dplyr::select(.iteration) |>
    dplyr::summarise(n = max(.iteration)) |>
    dplyr::collect() |>
    dplyr::pull(n)
}

ndraws.parquet_draws <- function(x, ...) {
  arrow::open_dataset(x$path) |>
    dplyr::count() |>
    dplyr::collect() |>
    dplyr::pull(n)
}

# Print method for better REPL experience
print.parquet_draws <- function(x, ..., n = 10) {
  params <- posterior::variables(x)
  n_params <- length(params)
  n_chains <- posterior::nchains(x)
  n_iters <- posterior::niterations(x)
  total_draws <- posterior::ndraws(x)

  cat("# A parquet_draws object: ", total_draws, " draws from ", n_chains, " chains\n", sep = "")
  cat("# Path: ", x$path, "\n", sep = "")
  cat("# ", n_params, " variables: ", sep = "")

  if (n_params <= 20) {
    cat(paste(params, collapse = ", "))
  } else {
    cat(paste(head(params, 20), collapse = ", "), ", ...")
  }
  cat("\n\n")

  # Show preview of first n rows
  preview <- arrow::open_dataset(x$path) |>
    dplyr::slice_head(n = n) |>
    dplyr::collect() |>
    posterior::as_draws_df()

  print(preview)

  if (total_draws > n) {
    cat("# ... with ", total_draws - n, " more draws\n", sep = "")
  }

  invisible(x)
}

# Convert back to regular draws_df - materializes entire dataset
as_draws_df.parquet_draws <- function(x, ...) {
  arrow::read_parquet(x$path) |>
    posterior::as_draws_df()
}

# Convert to draws_array
as_draws_array.parquet_draws <- function(x, ...) {
  x |>
    posterior::as_draws_df() |>
    posterior::as_draws_array()
}

# # Subset method for direct parameter access
# `[.parquet_draws` <- function(x, i, ...) {
#   if (!missing(i)) {
#     if (is.character(i)) {
#       draws_subset <- arrow::open_dataset(x$path) |>
#         select(all_of(i), .chain, .iteration, .draw) |>
#         collect() |>
#         as_draws_df()
#       return(draws_subset)
#     }
#   }
#   stop("Subsetting parquet_draws requires character vector of parameter names")
# }

# summarise_draws support - just materialize and delegate
summarise_draws.parquet_draws <- function(x, ...) {
  posterior::as_draws_df(x) |>
    posterior::summarise_draws(...)
}

# nolint end: object_usage_linter
