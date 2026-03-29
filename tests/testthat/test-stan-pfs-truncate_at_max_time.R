library(testthat)
library(cmdstanr)

test_that("truncate_at_max_time: all cases produce correct truncated pfs and right_censored", {
  cases <- list(
    list(
      case_name     = "no_truncation_needed",
      pfs           = c(5L, 10L, 15L),
      right_censored = c(0L, 0L, 0L),
      max_time      = c(20L, 20L, 20L)
    ),
    list(
      case_name     = "all_truncated",
      pfs           = c(25L, 30L, 40L),
      right_censored = c(0L, 0L, 0L),
      max_time      = c(20L, 20L, 20L)
    ),
    list(
      case_name     = "mixed_truncation_and_censoring",
      pfs           = c(5L, 15L, 25L, 35L),
      right_censored = c(0L, 1L, 0L, 1L),
      max_time      = c(10L, 20L, 20L, 30L)
    ),
    list(
      case_name     = "exact_boundary_not_truncated",
      pfs           = c(10L, 20L, 30L),
      right_censored = c(0L, 0L, 0L),
      max_time      = c(10L, 20L, 30L)
    ),
    list(
      case_name     = "already_censored_no_truncation",
      pfs           = c(8L, 12L),
      right_censored = c(1L, 1L),
      max_time      = c(20L, 20L)
    ),
    list(
      case_name     = "already_censored_truncation_needed",
      pfs           = c(25L, 30L),
      right_censored = c(1L, 1L),
      max_time      = c(20L, 20L)
    )
  )

  N_CASES <- length(cases)
  N       <- max(sapply(cases, function(x) length(x$pfs)))

  # Pad all cases to N patients using dummy values that won't be truncated
  pfs_mat <- matrix(1L, N_CASES, N)
  rc_mat  <- matrix(0L, N_CASES, N)
  mt_mat  <- matrix(100L, N_CASES, N)

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    n  <- length(cc$pfs)
    pfs_mat[i, seq_len(n)] <- cc$pfs
    rc_mat[i,  seq_len(n)] <- cc$right_censored
    mt_mat[i,  seq_len(n)] <- cc$max_time
  }

  stan_data <- list(
    N_CASES        = N_CASES,
    N              = N,
    pfs            = pfs_mat,
    right_censored = rc_mat,
    max_time       = mt_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_pfs_functions_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  for (i in seq_along(cases)) {
    cc       <- cases[[i]]
    expected <- r_truncate_at_max_time(cc$pfs, cc$right_censored, cc$max_time)
    n        <- length(cc$pfs)
    for (j in seq_len(n)) {
      expect_equal(
        get_val("out_pfs", i, j),
        expected$pfs[j],
        label = str_c("case ", i, " (", cc$case_name, ") out_pfs[", j, "]")
      )
      expect_equal(
        get_val("out_right_censored", i, j),
        expected$right_censored[j],
        label = str_c("case ", i, " (", cc$case_name, ") out_right_censored[", j, "]")
      )
    }
  }
})
