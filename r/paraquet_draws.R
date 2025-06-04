
# Define the parquet_draws class
new_parquet_draws <- function(path) {
  structure(
    list(path = path),
    class = c("parquet_draws", "draws_source")
  )
}

# {targets} format and functions ####

# Custom format for targets
tar_format_parquet_draws <- tar_format(
  read = function(path) {
    new_parquet_draws(path)
  },
  
  write = function(object, path) {
    arrow::write_parquet(object, path)
    path
  }
)

tar_parquet_draws <- function(name, fit, ...) {
  tar_target(name, fit$draws(format = "draws_df"), format = tar_format_parquet_draws, ...)
}

# Utility and overridden functions

.get_draws <- function(data, ...) {
  
  # Parse the ... to figure out which columns we need
  dots <- rlang::enquos(...)
  
  # Get all column names from the dataset
  ds <- arrow::open_dataset(data$path)
 
  if (length(dots) > 0) { 
    all_cols <- names(ds)
    param_names <- get_param_names_from_dots(dots, all_cols)
  
    # Read only the needed parameters
    ds |>  
      select(all_of(param_names), .chain, .iteration, .draw) |>  
      collect() |>  
      as_draws_df() 
  } else {
    as_draws_df(data)
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
    tidybayes::gather_rvars(draws_subset, ...)
}

# Print method for better REPL experience
print.parquet_draws <- function(x, ...) {
  ds <- arrow::open_dataset(x$path)
  n_draws <- ds |> count() |> collect() |> pull(n)
  params <- setdiff(names(ds), c(".chain", ".iteration", ".draw"))
  
  cat("Parquet-backed draws object\n")
  cat("Path:", x$path, "\n")
  cat("Parameters:", length(params), "\n")
  cat("Total draws:", n_draws, "\n")
  cat("\nAvailable parameters: ")
  if(length(params) <= 20) {
    cat(paste(params, collapse = ", "))
  } else {
    cat(paste(head(params, 20), collapse = ", "), ", ...")
  }
  cat("\n")
  invisible(x)
}

# Convert back to regular draws_df if needed
as_draws_df.parquet_draws <- function(x, ...) {
  arrow::open_dataset(x$path) |>
    collect() |>
    as_draws_df()
}

# Convert to draws_array
as_draws_array.parquet_draws <- function(x, ...) {
  x |>
    as_draws_df() |>
    as_draws_array()
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

# summarise_draws support
summarise_draws.parquet_draws <- function(x, ...) {
  # Convert to draws_df and use the default method
  as_draws_df(x) |>
    posterior::summarise_draws(...)
}

# Add methods for posterior diagnostics
rhat.parquet_draws <- function(x, ...) {
  .get_draws(x, ...) |> 
    posterior::rhat(...)
}

ess_bulk.parquet_draws <- function(x, ...) {
  .get_draws(x, ...) |> 
    posterior::ess_bulk(...)
}

ess_tail.parquet_draws <- function(x, ...) {
  .get_draws(x, ...) |> 
    posterior::ess_tail()
}