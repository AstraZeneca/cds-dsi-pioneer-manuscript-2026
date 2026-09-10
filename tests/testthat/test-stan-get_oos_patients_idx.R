library(testthat)
library(cmdstanr)
library(posterior)

test_that("get_oos_patients_idx: all cases produce correct OOS start indices", {
  cases <- list(
    list(
      case_name        = "single_patient_after_cutoff",
      n_patients       = 1L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(20L),
      cutoffs          = c(15L),
      expected         = c(1L)   # patient 1: last=20 > 15 → OOS from idx 1
    ),
    list(
      case_name        = "single_patient_before_cutoff",
      n_patients       = 1L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(10L),
      cutoffs          = c(15L),
      expected         = c(0L)   # patient 1: last=10 <= 15 → no OOS patients
    ),
    list(
      case_name        = "two_patients_split",
      n_patients       = 2L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(10L, 20L),
      cutoffs          = c(15L),
      expected         = c(2L)   # patient 1 skipped (10 <= 15), patient 2 is first OOS
    ),
    list(
      case_name        = "two_patients_both_oos",
      n_patients       = 2L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(20L, 30L),
      cutoffs          = c(15L),
      expected         = c(1L)   # patient 1: last=20 > 15 → OOS from idx 1
    ),
    list(
      case_name        = "exact_boundary_not_oos",
      n_patients       = 3L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(10L, 15L, 20L),
      cutoffs          = c(15L),
      expected         = c(3L)   # patients 1,2 have last <= 15; patient 3 (last=20) is first OOS
    ),
    list(
      case_name        = "two_cutoffs_three_patients",
      n_patients       = 3L,
      n_cutoffs        = 2L,
      sorted_last_visit = c(10L, 20L, 30L),
      cutoffs          = c(15L, 25L),
      expected         = c(2L, 3L)  # cutoff=15 → idx 2; cutoff=25 → idx 3 (continues from curr_patient_idx=2)
    )
  )

  N_CASES      <- length(cases)
  MAX_PATIENTS <- max(sapply(cases, function(x) x$n_patients))
  MAX_CUTOFFS  <- max(sapply(cases, function(x) x$n_cutoffs))

  sorted_last_visit_mat <- matrix(0L, N_CASES, MAX_PATIENTS)
  cutoffs_mat           <- matrix(0L, N_CASES, MAX_CUTOFFS)
  n_patients_vec        <- integer(N_CASES)
  n_cutoffs_vec         <- integer(N_CASES)

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    n_patients_vec[i] <- cc$n_patients
    n_cutoffs_vec[i]  <- cc$n_cutoffs
    sorted_last_visit_mat[i, seq_len(cc$n_patients)] <- cc$sorted_last_visit
    cutoffs_mat[i, seq_len(cc$n_cutoffs)]            <- cc$cutoffs
  }

  stan_data <- list(
    N_CASES    = N_CASES,
    MAX_PATIENTS = MAX_PATIENTS,
    MAX_CUTOFFS  = MAX_CUTOFFS,
    n_patients = n_patients_vec,
    n_cutoffs  = n_cutoffs_vec,
    sorted_last_visit_calendar_day = sorted_last_visit_mat,
    cutoff_calendar_day            = cutoffs_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_get_oos_patients_idx_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    for (c in seq_len(cc$n_cutoffs)) {
      expect_equal(
        get_val("oos_idx_out", i, c),
        cc$expected[c],
        label = str_c("case ", i, " (", cc$case_name, ") cutoff ", c)
      )
    }
  }
})
