# All-in-one testthat file for Stan util.stan (and pos.stan) utilities
# Follows the robust pattern in README.md

library(here)
library(testthat)
library(tidybayes)
library(posterior)
source(here("tests/testthat/helper-stan.R"))

# Define test cases for a subset of util.stan functions
cases <- list(
  # get_max_t: 2 patients, 1 tumor each, 2 measures each
  list(
    t_measure = c(1, 2, 3, 4),
    n_measures = c(2, 2),
    n_patient_tumors = c(1, 1),
    x = c(1, 2, 2, 3, 3, 1),
    n_x = 6,
    all = c(1, 2, 3, 2, 2, 4),
    what = c(2, 3),
    n_all = 6,
    n_what = 2,
    n_succ = 3,
    expect_max_t = c(2, 4),
    expect_num_unique = 3,
    expect_unique = c(1, 2, 3),
    expect_find_first = 2
  ),
  # num_unique/unique: 5 elements, some repeated
  list(
    t_measure = c(1, 2, 0, 0),
    n_measures = c(1, 1),
    n_patient_tumors = c(1, 1),
    x = c(5, 5, 5, 5, 5),
    n_x = 5,
    all = c(1, 1, 1, 1, 1, 1),
    what = c(1),
    n_all = 6,
    n_what = 1,
    n_succ = 2,
    expect_max_t = c(1, 2),
    expect_num_unique = 1,
    expect_unique = c(5),
    expect_find_first = 1
  )
)

N_CASES <- length(cases)
MAX_LEN <- max(sapply(cases, function(x) max(length(x$t_measure), length(x$n_measures), length(x$n_patient_tumors), length(x$x), length(x$all), length(x$what))))

# Pad all arrays to MAX_LEN
for (i in seq_along(cases)) {
  for (nm in c("t_measure", "n_measures", "n_patient_tumors", "x", "all", "what")) {
    cases[[i]][[nm]] <- c(cases[[i]][[nm]], rep(0, MAX_LEN - length(cases[[i]][[nm]])))
  }
}

data_list <- list(
  N_CASES = N_CASES,
  MAX_LEN = MAX_LEN,
  t_measure = do.call(rbind, lapply(cases, `[[`, "t_measure")),
  n_measures = do.call(rbind, lapply(cases, `[[`, "n_measures")),
  n_patient_tumors = do.call(rbind, lapply(cases, `[[`, "n_patient_tumors")),
  x = do.call(rbind, lapply(cases, `[[`, "x")),
  n_x = sapply(cases, `[[`, "n_x"),
  all = do.call(rbind, lapply(cases, `[[`, "all")),
  what = do.call(rbind, lapply(cases, `[[`, "what")),
  n_all = sapply(cases, `[[`, "n_all"),
  n_what = sapply(cases, `[[`, "n_what"),
  n_succ = sapply(cases, `[[`, "n_succ")
)


fit <- test_stan_function(
  stan_file = here("tests/testthat/stan/test_util_all.stan"),
  data = data_list
)

# Use posterior::as_draws_df for robust extraction
library(posterior)
draws_df <- as_draws_df(fit$draws())

# Helper to extract array elements from draws_df
get_array <- function(df, prefix, dims) {
  # dims: vector of dimension sizes, e.g. c(N_CASES, MAX_LEN)
  arr <- array(NA_integer_, dims)
  idx <- 1
  for (i in seq_len(dims[1])) {
    if (length(dims) == 1) {
      name <- sprintf("%s[%d]", prefix, i)
      arr[i] <- as.integer(df[[name]])
    } else {
      for (j in seq_len(dims[2])) {
        name <- sprintf("%s[%d,%d]", prefix, i, j)
        arr[i, j] <- as.integer(df[[name]])
      }
    }
  }
  arr
}


max_t_out <- get_array(draws_df, "max_t_out", c(N_CASES, MAX_LEN))
num_unique_out <- get_array(draws_df, "num_unique_out", c(N_CASES))
unique_out <- get_array(draws_df, "unique_out", c(N_CASES, MAX_LEN))
find_first_out <- get_array(draws_df, "find_first_out", c(N_CASES))

# Extract and check outputs for each case
for (i in seq_along(cases)) {
  expect_equal(max_t_out[i, 1:length(cases[[i]]$expect_max_t)], cases[[i]]$expect_max_t)
  expect_equal(num_unique_out[i], cases[[i]]$expect_num_unique)
  expect_equal(unique_out[i, 1:length(cases[[i]]$expect_unique)], cases[[i]]$expect_unique)
  expect_equal(find_first_out[i], cases[[i]]$expect_find_first)
}
