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
    x = c(1, 2, 3),
    n_x = 3,
    all = c(1, 2, 3, 2, 2, 4),
    what = c(2, 3),
    n_all = 6,
    n_what = 2,
    n_succ = 3,
    expect_max_t = c(2, 4),
    expect_num_unique = 3,
    expect_unique = c(1, 2, 3),
    expect_find_first = 0
  ),
  # num_unique/unique: 2 elements, both unique
  list(
    t_measure = c(1, 2),
    n_measures = c(1, 1),
    n_patient_tumors = c(1, 1),
    x = c(5, 6),
    n_x = 2,
    all = c(1, 1),
    what = c(1),
    n_all = 2,
    n_what = 1,
    n_succ = 2,
    expect_max_t = c(1, 2),
    expect_num_unique = 2,
    expect_unique = c(5, 6),
    expect_find_first = 1
  )
)

N_CASES <- length(cases)
MAX_LEN <- max(sapply(cases, function(x) {
  max(
    length(x$t_measure),
    length(x$n_measures),
    length(x$n_patient_tumors),
    length(x$x),
    length(x$all),
    length(x$what)
  )
}))

# Pad all arrays to MAX_LEN
for (i in seq_along(cases)) {
  for (nm in c(
    "t_measure",
    "n_measures",
    "n_patient_tumors",
    "x",
    "all",
    "what"
  )) {
    cases[[i]][[nm]] <- c(
      cases[[i]][[nm]],
      rep(0, MAX_LEN - length(cases[[i]][[nm]]))
    )
  }
}

data_list <- list(
  N_CASES = N_CASES,
  MAX_LEN = MAX_LEN,
  t_measure = do.call(rbind, lapply(cases, `[[`, "t_measure")),
  n_measures = do.call(rbind, lapply(cases, `[[`, "n_measures")),
  n_patient_tumors = do.call(rbind, lapply(cases, `[[`, "n_patient_tumors")),
  x = do.call(rbind, lapply(cases, `[[`, "x")),
  n_x = sapply(cases, function(x) length(x$expect_unique)),
  all = do.call(rbind, lapply(cases, `[[`, "all")),
  what = do.call(rbind, lapply(cases, `[[`, "what")),
  n_all = sapply(cases, function(x) length(x$expect_max_t)),
  n_what = sapply(cases, function(x) length(x$expect_unique)),
  n_succ = sapply(cases, `[[`, "n_succ"),
  # Explicit per-case dimensions for Stan test harness
  n_patients_case = sapply(cases, function(x) length(x$expect_max_t)),
  n_tumors_case = sapply(cases, function(x) length(x$n_measures)),
  n_meas_sum_case = sapply(cases, function(x) sum(x$n_measures))
)


fit <- test_stan_function(
  stan_file = here("tests/testthat/stan/test_util_all.stan"),
  data = data_list
)

# Direct extraction from fit$draws()
library(posterior)
draws_array <- as_draws_array(fit$draws())

# Helper to extract a vector for a given variable and case
extract_vec <- function(var, case, max_len = MAX_LEN) {
  as.numeric(draws_array[1, 1, str_c(var, "[", case, ",", 1:max_len, "]")])
}

# Helper to extract a scalar for a given variable and case
extract_scalar <- function(var, case) {
  as.numeric(draws_array[1, 1, str_c(var, "[", case, "]")])
}

# Helper to extract a vector for a given variable (no case)
extract_vec_nocase <- function(var, len) {
  as.numeric(draws_array[1, 1, str_c(var, "[", 1:len, "]")])
}

# Helper to extract a scalar variable (no case)
extract_scalar_nocase <- function(var) {
  as.numeric(draws_array[1, 1, var])
}

first_draw <- list(
  max_t_out = lapply(1:N_CASES, function(i) extract_vec("max_t_out", i)),
  num_unique_out = sapply(1:N_CASES, function(i) {
    extract_scalar("num_unique_out", i)
  }),
  unique_out = lapply(1:N_CASES, function(i) extract_vec("unique_out", i)),
  find_first_out = sapply(1:N_CASES, function(i) {
    extract_scalar("find_first_out", i)
  }),
  min_eig = extract_scalar_nocase("min_eig"),
  max_eig = extract_scalar_nocase("max_eig"),
  cond_num = extract_scalar_nocase("cond_num"),
  n_missing_measures = extract_vec_nocase("n_missing_measures", 2),
  idx0 = extract_vec_nocase("idx0", 2),
  idx1 = extract_vec_nocase("idx1", 3),
  uniq_vals = extract_vec_nocase("uniq_vals", 4),
  uniq_pos = extract_vec_nocase("uniq_pos", 3),
  idx_dict = extract_vec_nocase("idx_dict", 7),
  mean_tumor = extract_scalar_nocase("mean_tumor"),
  sd_tumor = extract_scalar_nocase("sd_tumor"),
  std_vals = extract_vec_nocase("std_vals", 3)
)
test_that("util: max_t, unique, and find_first utilities for all cases", {
  for (i in seq_along(cases)) {
    actual_max_t <- as.numeric(first_draw$max_t_out[[i]])
    actual_max_t <- actual_max_t[actual_max_t != 0]
    expect_equal(actual_max_t, cases[[i]]$expect_max_t,
                 label = sprintf("max_t case %d", i))

    expect_equal(first_draw$num_unique_out[[i]], cases[[i]]$expect_num_unique,
                 label = sprintf("num_unique case %d", i))

    actual_unique <- as.numeric(first_draw$unique_out[[i]])
    actual_unique <- actual_unique[actual_unique != 0]
    expect_equal(actual_unique, cases[[i]]$expect_unique,
                 label = sprintf("unique values case %d", i))

    expect_equal(first_draw$find_first_out[[i]], cases[[i]]$expect_find_first,
                 label = sprintf("find_first case %d", i))
  }
})

# Expanded util function checks
test_that("summarize_matrix_eigenvalues returns correct values", {
  expect_equal(first_draw$min_eig, 2)
  expect_equal(first_draw$max_eig, 8)
  expect_equal(first_draw$cond_num, 4)
})

test_that("calculate_n_missing_measures returns correct values", {
  expect_equal(first_draw$n_missing_measures, c(0, 0))
})

test_that("get_mask_idx returns correct indices", {
  expect_equal(first_draw$idx0, c(1, 3))
  expect_equal(first_draw$idx1, c(2, 4, 5))
})

test_that("unique_by_pos returns correct values and positions", {
  expect_equal(first_draw$uniq_vals, c(1, 2, 3, 4))
  expect_equal(first_draw$uniq_pos, c(1, 3, 5))
})

test_that("get_idx_dict returns correct mapping", {
  expect_equal(first_draw$idx_dict, c(0, 1, 0, 2, 0, 0, 4))
})

test_that("standardize_tumor_sizes returns correct mean, sd, and standardized values", {
  expect_equal(first_draw$mean_tumor, 2)
  expect_equal(first_draw$sd_tumor, 1)
  expect_equal(first_draw$std_vals, c(-1, 0, 1))
})
