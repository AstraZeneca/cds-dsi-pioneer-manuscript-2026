# All-in-one testthat file for Stan pos.stan utilities
# Follows the robust pattern in README.md

library(here)
library(testthat)
library(tidybayes)
library(posterior)
source(here("tests/testthat/helper-stan.R"))

# Define test cases for position array utilities
cases <- list(
  list(
    group_sizes = c(2, 2, 3),
    flat_data = c(10, 20, 30, 40, 50, 60, 70),
    n_groups = 3,
    n_flat = 7,
    expect_pos = c(1, 3, 5, 8),
    expect_pos_size = c(2, 2, 3),
    expect_max = c(20, 40, 70),
    expect_min = c(10, 30, 50),
    expect_max_idx = c(11, 11, 21),
    expect_min_pos = c(10, 30, 50),
    expect_max_pos = c(20, 40, 70),
    expect_last_int = c(20, 40, 70)
  )
)


# Pad only group_sizes and flat_data, never the expect_* fields
max_groups <- max(sapply(cases, function(x) length(x$group_sizes)))
max_size <- max(sapply(cases, function(x) length(x$flat_data)))

for (i in seq_along(cases)) {
  cases[[i]]$group_sizes <- c(cases[[i]]$group_sizes, rep(0, max_groups - length(cases[[i]]$group_sizes)))
  cases[[i]]$flat_data <- c(cases[[i]]$flat_data, rep(0, max_size - length(cases[[i]]$flat_data)))
}

data_list <- list(
  N_CASES = length(cases),
  MAX_GROUPS = max_groups,
  MAX_SIZE = max_size,
  group_sizes = do.call(rbind, lapply(cases, function(x) x$group_sizes)),
  flat_data = do.call(rbind, lapply(cases, function(x) x$flat_data)),
  n_groups = sapply(cases, function(x) x$n_groups),
  n_flat = sapply(cases, function(x) x$n_flat)
)

stan_file <- "tests/testthat/stan/test_pos_all.stan"


fit <- test_stan_function(
  stan_file,
  data = data_list
)


# Use tidybayes to extract and reshape Stan outputs robustly
res_df <- as_draws_df(fit$draws())

# Helper to extract first draw for a variable (matrix output)
extract_first_draw_mat <- function(df, var, idx1, idx2) {
  out <- spread_draws(df, !!as.name(var)[[idx1]], !!as.name(var)[[idx2]])
  # Get first draw for each [idx1, idx2] pair
  mat <- with(out, tapply(.data[[var]], list(.data[[idx1]], .data[[idx2]]), function(x) x[1]))
  return(mat)
}

# Helper to extract first draw for a variable (vector output)
extract_first_draw_vec <- function(df, var, idx1) {
  out <- spread_draws(df, !!as.name(var)[[idx1]])
  vec <- tapply(out[[var]], out[[idx1]], function(x) x[1])
  return(as.numeric(vec))
}

# Extract all outputs as matrices/vectors for easy indexing
res <- list(
  pos_out = spread_draws(res_df, pos_out[case, group]) |> with(tapply(pos_out, list(case, group), function(x) x[1])),
  pos_size_out = spread_draws(res_df, pos_size_out[case, group]) |> with(tapply(pos_size_out, list(case, group), function(x) x[1])),
  max_out = spread_draws(res_df, max_out[case, group]) |> with(tapply(max_out, list(case, group), function(x) x[1])),
  min_out = spread_draws(res_df, min_out[case, group]) |> with(tapply(min_out, list(case, group), function(x) x[1])),
  max_idx_out = spread_draws(res_df, max_idx_out[case, group]) |> with(tapply(max_idx_out, list(case, group), function(x) x[1])),
  min_pos_out = spread_draws(res_df, min_pos_out[case, group]) |> with(tapply(min_pos_out, list(case, group), function(x) x[1])),
  max_pos_out = spread_draws(res_df, max_pos_out[case, group]) |> with(tapply(max_pos_out, list(case, group), function(x) x[1])),
  last_int_out = spread_draws(res_df, last_int_out[case, group]) |> with(tapply(last_int_out, list(case, group), function(x) x[1]))
)

# Check outputs for each case
for (i in seq_along(cases)) {
  case <- cases[[i]]
  # cat(sprintf("\n--- CASE %d ---\n", i))
  # cat("full pos_out:    actual=", toString(res$pos_out[i,]), " expected=", toString(case$expect_pos), "\n")
  # cat("pos_out:    actual=", toString(res$pos_out[i, 1:(case$n_groups+1)]), " expected=", toString(case$expect_pos), "\n")
  # cat("pos_size:   actual=", toString(res$pos_size_out[i, 1:case$n_groups]), " expected=", toString(case$expect_pos_size), "\n")
  # cat("max:        actual=", toString(res$max_out[i, 1:case$n_groups]), " expected=", toString(case$expect_max), "\n")
  # cat("min:        actual=", toString(res$min_out[i, 1:case$n_groups]), " expected=", toString(case$expect_min), "\n")
  # cat("max_idx:    actual=", toString(res$max_idx_out[i, 1:case$n_groups]), " expected=", toString(case$expect_max_idx), "\n")
  # cat("min_pos:    actual=", toString(res$min_pos_out[i, 1:case$n_groups]), " expected=", toString(case$expect_min_pos), "\n")
  # cat("max_pos:    actual=", toString(res$max_pos_out[i, 1:case$n_groups]), " expected=", toString(case$expect_max_pos), "\n")
  # cat("last_int:   actual=", toString(res$last_int_out[i, 1:case$n_groups]), " expected=", toString(case$expect_last_int), "\n")
  test_that(sprintf("pos.stan utilities: case %d", i), {
    expect_equal(unname(res$pos_out[i, 1:(case$n_groups+1)]), case$expect_pos)
    expect_equal(unname(res$pos_size_out[i, 1:case$n_groups]), case$expect_pos_size)
    expect_equal(unname(res$max_out[i, 1:case$n_groups]), case$expect_max)
    expect_equal(unname(res$min_out[i, 1:case$n_groups]), case$expect_min)
    expect_equal(unname(res$max_idx_out[i, 1:case$n_groups]), case$expect_max_idx)
    expect_equal(unname(res$min_pos_out[i, 1:case$n_groups]), case$expect_min_pos)
    expect_equal(unname(res$max_pos_out[i, 1:case$n_groups]), case$expect_max_pos)
    expect_equal(unname(res$last_int_out[i, 1:case$n_groups]), case$expect_last_int)
  })
}
