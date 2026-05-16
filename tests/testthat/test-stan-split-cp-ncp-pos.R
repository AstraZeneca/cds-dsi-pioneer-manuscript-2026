library(testthat)
library(here)
library(posterior)

# R oracle: route a group to _raw_ if mode in {1 FE, 2 RE, 3 RE_GP}; route to
# _cp_ if mode == 4 RE_CP. Skip levels with mode == 0 NONE.
r_split_cp_ncp_pos <- function(mode, n_groups) {
  n_levels <- length(mode)
  raw_counts <- integer(n_levels)
  cp_counts  <- integer(n_levels)
  for (lv in seq_len(n_levels)) {
    if (mode[lv] == 4L) {
      cp_counts[lv] <- n_groups[lv]
    } else if (mode[lv] %in% c(1L, 2L, 3L)) {
      raw_counts[lv] <- n_groups[lv]
    }
  }
  list(
    n_raw   = sum(raw_counts),
    n_cp    = sum(cp_counts),
    raw_pos = as.integer(c(1L, cumsum(raw_counts) + 1L)),
    cp_pos  = as.integer(c(1L, cumsum(cp_counts)  + 1L))
  )
}

run_split_test <- function(mode, n_groups) {
  n_levels <- length(mode)
  stan_data <- list(
    n_levels = n_levels,
    mode     = as.array(as.integer(mode)),
    n_groups = as.array(as.integer(n_groups))
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_split_cp_ncp_pos.stan"),
    data = stan_data
  )
  d <- as_draws_df(fit$draws())
  exp <- r_split_cp_ncp_pos(mode, n_groups)
  expect_equal(as.integer(d$n_raw[1]), exp$n_raw)
  expect_equal(as.integer(d$n_cp[1]),  exp$n_cp)
  for (i in seq_len(n_levels + 1)) {
    expect_equal(as.integer(d[[paste0("raw_pos[", i, "]")]][1]), exp$raw_pos[i])
    expect_equal(as.integer(d[[paste0("cp_pos[",  i, "]")]][1]), exp$cp_pos[i])
  }
}

test_that("split: all NONE", {
  run_split_test(c(0L, 0L, 0L), c(5L, 10L, 100L))
})

test_that("split: all RE (legacy)", {
  run_split_test(c(2L, 2L, 2L), c(5L, 10L, 100L))
})

test_that("split: all RE_CP", {
  run_split_test(c(4L, 4L, 4L), c(5L, 10L, 100L))
})

test_that("split: mixed NONE/FE/RE/RE_CP 4-level", {
  run_split_test(c(0L, 1L, 2L, 4L), c(3L, 5L, 10L, 100L))
})

test_that("split: FE goes to raw bucket", {
  run_split_test(c(1L, 4L), c(5L, 50L))
})

test_that("split: RE_GP goes to raw bucket", {
  run_split_test(c(3L, 4L), c(5L, 50L))
})
