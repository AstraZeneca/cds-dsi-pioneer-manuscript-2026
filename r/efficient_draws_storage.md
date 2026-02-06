# Efficient Stan Draws Storage Solutions

## The Core Problem
Stan outputs CSV files, and every time you access the fit object with `cmdstanr_format`, it re-reads those CSVs. When multiple targets need the draws, this becomes very slow.

## Solutions (Ranked by Efficiency)

### 1. **Use CmdStan's output_dir with persistence** (Fastest, no conversion)
Keep Stan CSVs but use `targets` to avoid redundant reads:

```r
# Store fit once
tar_target(my_fit, 
  model$sample(data = stan_data, ...),
  format = cmdstanr_format
)

# Extract draws ONCE to an efficient format
tar_target(my_draws,
  my_fit$draws(format = "draws_rvars"),  # rvar format is memory efficient
  format = "qs2"  # qs2 is fast compressed format
)

# All downstream targets use my_draws (loaded once by targets)
tar_target(analysis1, do_analysis(my_draws))
tar_target(analysis2, do_other(my_draws))
```

**Pros**: Simple, no conversion overhead
**Cons**: First load still reads all CSVs, but targets caches it

### 2. **Parquet with Partitioned Storage** (Best for huge models)
Use partitioned parquet files - one per chain:

```r
# In paraquet_draws.R - simpler approach
convert_fit_to_partitioned_parquet <- function(fit, output_dir) {
  csv_files <- fit$output_files()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get metadata line count
  comment_lines <- sum(startsWith(readLines(csv_files[1], n = 100), "#"))
  
  # Convert each CSV to its own parquet (parallel if desired)
  purrr::walk2(csv_files, seq_along(csv_files), function(csv, i) {
    arrow::read_csv_arrow(csv, skip = comment_lines, as_data_frame = FALSE) |>
      arrow::write_parquet(file.path(output_dir, sprintf("chain_%d.parquet", i)))
  })
  
  output_dir
}

# Then open as dataset (lazy evaluation)
my_draws <- arrow::open_dataset(partition_dir, format = "parquet")
```

**Pros**: Lazily loads only what you need, can filter before loading
**Cons**: Still needs to read CSVs once for initial conversion

### 3. **CmdStanR Direct Parquet Output** (Future solution)
Tell CmdStan to output parquet directly (requires Stan 2.35+):

```r
fit <- model$sample(
  data = stan_data,
  output_format = "parquet"  # Not yet available!
)
```

### 4. **Use $save_output_files() + Memoization** (Current best practice)
```r
# Save CSVs to permanent location
sample_and_save <- function(model, data, output_dir, ...) {
  fit <- model$sample(data = data, ...)
  fit$save_output_files(dir = output_dir, random = FALSE)
  fit
}

# Use posterior::read_cmdstan_csv with caching
library(memoise)
read_draws_cached <- memoise(function(csv_files) {
  posterior::read_cmdstan_csv(csv_files, format = "draws_df")
})

# In targets:
tar_target(my_fit, sample_and_save(model, data, "outputs/fits/myfit"))
tar_target(my_draws, read_draws_cached(my_fit$output_files()), format = "qs2")
```

## Recommendation for Your Project

Since you have `qs2` format already and are using `targets`:

```r
# 1. Keep cmdstanr_format for fit objects (for diagnostics)
tar_target(tumor_fit,
  model$sample(...),
  format = cmdstanr_format
)

# 2. Extract draws ONCE using qs2 (which is already very fast)
tar_target(tumor_draws,
  tumor_fit$draws(format = "draws_df"),
  format = "qs2",  # Your existing fast format!
  resources = tar_resources(qs = tar_resources_qs(preset = "fast"))
)

# 3. Use draws in downstream targets
tar_target(posterior_summary,
  tumor_draws |> spread_rvars(beta[i])
)
```

**This is actually faster than parquet for most use cases!** `qs2` with preset="fast" is:
- Faster than parquet for R native objects
- Handles `draws_df` natively
- Already integrated with targets
- No conversion overhead

## When to Use Parquet

Use parquet when:
1. You have extremely large models (>1GB of draws)
2. You need to share data with non-R tools
3. You want selective parameter loading (though qs2 can do this too with some work)
4. You're working with distributed computing

For most Bayesian workflows with targets, **qs2 format is the winner**.
