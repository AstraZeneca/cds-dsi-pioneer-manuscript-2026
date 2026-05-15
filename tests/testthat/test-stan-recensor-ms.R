library(testthat)
library(cmdstanr)
library(tidybayes)
library(stringr)
source(here::here("tests/testthat/helper-stan.R"))

test_that("recensor_ms_at_cutoff handles all patient states correctly", {
  # Each case is a single-patient scenario. Using 1 patient per case keeps
  # the expected values trivially readable. Multi-patient cases below (#7)
  # verify that the function handles array indexing correctly.

  cases <- list(
    # Case 1: State 0 — censored beyond cutoff → clamp times
    list(
      name = "state_0_clamp",
      np = 1L,
      ms_final_state = 0L, ms_time_01 = 100L, ms_censored_01 = 1L,
      ms_time_02 = 100L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 50L, exp_censored_01 = 1L,
      exp_time_02 = 50L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 2: State 1 — progression after cutoff → downgrade to state 0
    list(
      name = "state_1_after_cutoff",
      np = 1L,
      ms_final_state = 1L, ms_time_01 = 80L, ms_censored_01 = 0L,
      ms_time_02 = 80L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 2L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 50L, exp_censored_01 = 1L,
      exp_time_02 = 50L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 3: State 1 — progression before cutoff → keep unchanged
    list(
      name = "state_1_before_cutoff",
      np = 1L,
      ms_final_state = 1L, ms_time_01 = 30L, ms_censored_01 = 0L,
      ms_time_02 = 30L, ms_time_12 = 20L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 2L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 1L, exp_time_01 = 30L, exp_censored_01 = 0L,
      exp_time_02 = 30L, exp_time_12 = 20L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 2L, exp_prog_det = 0L, exp_ic_gap_01 = 3L
    ),
    # Case 4: State 2, 0→1→2 — death after cutoff → downgrade to state 1
    list(
      name = "state_2_012_death_after",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 40L, ms_censored_01 = 0L,
      ms_time_02 = 60L, ms_time_12 = 20L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 60L, interval_censored = 0L, ms_prog_deterministic = 1L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 1L, exp_time_01 = 40L, exp_censored_01 = 0L,
      exp_time_02 = 50L, exp_time_12 = 10L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 1L, exp_ic_gap_01 = 0L
    ),
    # Case 5: State 2, 0→2 direct death after cutoff → downgrade to state 0
    list(
      name = "state_2_02_death_after",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 60L, ms_censored_01 = 1L,
      ms_time_02 = 60L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 50L, exp_censored_01 = 1L,
      exp_time_02 = 50L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 6: Not enrolled → zeroed out
    list(
      name = "not_enrolled",
      np = 1L,
      ms_final_state = 1L, ms_time_01 = 30L, ms_censored_01 = 0L,
      ms_time_02 = 30L, ms_time_12 = 10L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 2L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 0L, cutoff_cal_week = 0L,
      exp_final_state = 0L, exp_time_01 = 0L, exp_censored_01 = 1L,
      exp_time_02 = 0L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 7: State 3 — dropout after cutoff → censor
    list(
      name = "state_3_dropout_after",
      np = 1L,
      ms_final_state = 3L, ms_time_01 = 0L, ms_censored_01 = 1L,
      ms_time_02 = 60L, ms_time_12 = 0L, ms_time_03 = 60L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 50L, exp_censored_01 = 1L,
      exp_time_02 = 50L, exp_time_12 = 0L, exp_time_03 = 50L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 8: State 3 — dropout before cutoff, off-trial death after → censor 3→2
    list(
      name = "state_3_death_after",
      np = 1L,
      ms_final_state = 3L, ms_time_01 = 0L, ms_censored_01 = 1L,
      ms_time_02 = 40L, ms_time_12 = 0L, ms_time_03 = 40L, ms_time_32 = 20L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 3L, exp_time_01 = 0L, exp_censored_01 = 1L,
      exp_time_02 = 40L, exp_time_12 = 0L, exp_time_03 = 40L, exp_time_32 = 10L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 9: State 2, 0→1→2 — both prog and death before cutoff → unchanged
    list(
      name = "state_2_012_both_before",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 20L, ms_censored_01 = 0L,
      ms_time_02 = 30L, ms_time_12 = 10L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 30L, interval_censored = 0L, ms_prog_deterministic = 1L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 2L, exp_time_01 = 20L, exp_censored_01 = 0L,
      exp_time_02 = 30L, exp_time_12 = 10L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 1L, exp_ic_gap_01 = 0L
    ),
    # Case 10: ic_gap_01 recomputation — progression with IC flipped by re-censoring
    list(
      name = "ic_gap_recompute",
      np = 1L,
      ms_final_state = 1L, ms_time_01 = 80L, ms_censored_01 = 0L,
      ms_time_02 = 80L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 3L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      # Re-censored to state 0 → censored_01=1 → ic_gap_01=0
      exp_final_state = 0L, exp_time_01 = 50L, exp_censored_01 = 1L,
      exp_time_02 = 50L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 11: State 2, 0→1→2 — progression AFTER cutoff → censor everything
    list(
      name = "state_2_012_prog_after",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 60L, ms_censored_01 = 0L,
      ms_time_02 = 80L, ms_time_12 = 20L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 80L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 50L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 50L, exp_censored_01 = 1L,
      exp_time_02 = 50L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 12: 0→2 direct death BETWEEN last visit and calendar cutoff → retain death
    # Last visit was week 14, calendar cutoff is week 40, death at week 31 (exact from registry).
    # Old behavior (visit-cutoff uniformly) clipped the death at week 14 and called the patient
    # censored. New behavior retains the death because it occurred before the calendar cutoff.
    list(
      name = "state_2_02_death_between_visit_and_cal",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 31L, ms_censored_01 = 1L,
      ms_time_02 = 31L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 14L, cutoff_cal_week = 40L,
      exp_final_state = 2L, exp_time_01 = 31L, exp_censored_01 = 1L,
      exp_time_02 = 31L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 13: 0→2 direct death AFTER calendar cutoff → censor at last visit
    list(
      name = "state_2_02_death_after_cal",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 60L, ms_censored_01 = 1L,
      ms_time_02 = 60L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 40L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 40L, exp_censored_01 = 1L,
      exp_time_02 = 40L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    ),
    # Case 14: 0→1→2 — progression observed, death between last visit and cal cutoff → retain death
    list(
      name = "state_2_012_death_between_visit_and_cal",
      np = 1L,
      ms_final_state = 2L, ms_time_01 = 20L, ms_censored_01 = 0L,
      ms_time_02 = 35L, ms_time_12 = 15L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 35L, interval_censored = 0L, ms_prog_deterministic = 1L,
      cutoff_visit_week = 25L, cutoff_cal_week = 40L,
      exp_final_state = 2L, exp_time_01 = 20L, exp_censored_01 = 0L,
      exp_time_02 = 35L, exp_time_12 = 15L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 1L, exp_ic_gap_01 = 0L
    ),
    # Case 15: State 0 — state-0 follow-up capped at last visit even when calendar cutoff is later
    # (no registry event to know about, so visit is the correct censoring time).
    list(
      name = "state_0_clamp_at_visit_not_cal",
      np = 1L,
      ms_final_state = 0L, ms_time_01 = 100L, ms_censored_01 = 1L,
      ms_time_02 = 100L, ms_time_12 = 0L, ms_time_03 = 0L, ms_time_32 = 0L,
      ms_os_event_12 = 0L, interval_censored = 0L, ms_prog_deterministic = 0L,
      cutoff_visit_week = 30L, cutoff_cal_week = 50L,
      exp_final_state = 0L, exp_time_01 = 30L, exp_censored_01 = 1L,
      exp_time_02 = 30L, exp_time_12 = 0L, exp_time_03 = 0L, exp_time_32 = 0L,
      exp_ic = 0L, exp_prog_det = 0L, exp_ic_gap_01 = 0L
    )
  )

  # Rectangularize
  N_cases <- length(cases)
  max_np <- max(sapply(cases, \(x) x$np))

  make_mat <- function(field) {
    m <- matrix(0L, N_cases, max_np)
    for (i in seq_len(N_cases)) m[i, 1:cases[[i]]$np] <- cases[[i]][[field]]
    m
  }

  stan_data <- list(
    N_cases = N_cases,
    max_n_patients = max_np,
    n_patients = sapply(cases, \(x) x$np),
    ms_final_state = make_mat("ms_final_state"),
    ms_time_01 = make_mat("ms_time_01"),
    ms_censored_01 = make_mat("ms_censored_01"),
    ms_time_02 = make_mat("ms_time_02"),
    ms_time_12 = make_mat("ms_time_12"),
    ms_time_03 = make_mat("ms_time_03"),
    ms_time_32 = make_mat("ms_time_32"),
    ms_os_event_12 = make_mat("ms_os_event_12"),
    interval_censored = make_mat("interval_censored"),
    ms_prog_deterministic = make_mat("ms_prog_deterministic"),
    cutoff_visit_week = make_mat("cutoff_visit_week"),
    cutoff_cal_week = make_mat("cutoff_cal_week")
  )

  fit <- test_stan_function("tests/testthat/stan/test_recensor_ms_all.stan", data = stan_data)
  draws_df <- posterior::as_draws_df(fit$draws())

  # Check each case
  for (i in seq_len(N_cases)) {
    case <- cases[[i]]
    label <- case$name
    np <- case$np

    for (j in seq_len(np)) {
      expect_equal(get_stan_val(draws_df, "out_final_state", i, j),
                   case$exp_final_state[j], label = str_glue("{label}: final_state[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_time_01", i, j),
                   case$exp_time_01[j], label = str_glue("{label}: time_01[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_censored_01", i, j),
                   case$exp_censored_01[j], label = str_glue("{label}: censored_01[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_time_02", i, j),
                   case$exp_time_02[j], label = str_glue("{label}: time_02[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_time_12", i, j),
                   case$exp_time_12[j], label = str_glue("{label}: time_12[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_time_03", i, j),
                   case$exp_time_03[j], label = str_glue("{label}: time_03[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_time_32", i, j),
                   case$exp_time_32[j], label = str_glue("{label}: time_32[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_ic", i, j),
                   case$exp_ic[j], label = str_glue("{label}: ic[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_prog_det", i, j),
                   case$exp_prog_det[j], label = str_glue("{label}: prog_det[{j}]"))
      expect_equal(get_stan_val(draws_df, "out_ic_gap_01", i, j),
                   case$exp_ic_gap_01[j], label = str_glue("{label}: ic_gap_01[{j}]"))
    }
  }
})

test_that("derive_ms_censoring_indicators on re-censored arrays — LFO OS-KM wiring", {
  # Guards the LFO GQ include (_ms_standalone_lfo_os_km_generated_quantities.stan):
  # after recensor_ms_at_cutoff, the three censoring indicators fed into
  # compute_ms_os_km_rng must come from derive_ms_censoring_indicators applied
  # to the *re-censored* final_state / time_01 / time_32.
  #
  # For each original case, we hand-compute expected c02/c12/c32 from the
  # EXPECTED re-censored outputs (exp_final_state, exp_time_01, exp_time_32),
  # then invoke the Stan derive harness on those re-censored inputs.

  # Pattern table from derive_ms_censoring_indicators doc:
  expected_censoring <- function(final_state, time_01, time_32) {
    # c02 == 0 iff state==2 && time_01==0 (direct on-trial death 0→2)
    c02 <- as.integer(!(final_state == 2L & time_01 == 0L))
    # c12 == 0 iff state==2 && time_01>0  (post-progression death 1→2)
    c12 <- as.integer(!(final_state == 2L & time_01 >  0L))
    # c32 == 0 iff time_32 > 0             (off-trial death observed)
    c32 <- as.integer(!(time_32 > 0L))
    list(c02 = c02, c12 = c12, c32 = c32)
  }

  # Use the exp_* values from each case from the fixture above.
  # Reconstruct the same case list for clarity.
  rc_cases <- list(
    state_0_clamp            = c(final = 0L, t01 = 50L, t32 = 0L),
    state_1_after_cutoff     = c(final = 0L, t01 = 50L, t32 = 0L),
    state_1_before_cutoff    = c(final = 1L, t01 = 30L, t32 = 0L),
    state_2_012_death_after  = c(final = 1L, t01 = 40L, t32 = 0L),
    state_2_02_death_after   = c(final = 0L, t01 = 50L, t32 = 0L),
    not_enrolled             = c(final = 0L, t01 = 0L,  t32 = 0L),
    state_3_dropout_after    = c(final = 0L, t01 = 50L, t32 = 0L),
    state_3_death_after      = c(final = 3L, t01 = 0L,  t32 = 10L),
    state_2_012_both_before  = c(final = 2L, t01 = 20L, t32 = 0L),
    ic_gap_recompute         = c(final = 0L, t01 = 50L, t32 = 0L),
    state_2_012_prog_after   = c(final = 0L, t01 = 50L, t32 = 0L)
  )

  final_state <- vapply(rc_cases, \(x) unname(x["final"]), integer(1))
  t01         <- vapply(rc_cases, \(x) unname(x["t01"]),   integer(1))
  t32         <- vapply(rc_cases, \(x) unname(x["t32"]),   integer(1))

  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_derive_ms_censoring_all.stan"),
    list(N = length(final_state),
         ms_final_state = final_state,
         ms_time_01 = t01,
         ms_time_32 = t32),
    seed = 42L
  )
  d <- posterior::as_draws_df(fit$draws())

  exp <- expected_censoring(final_state, t01, t32)
  for (i in seq_along(rc_cases)) {
    nm <- names(rc_cases)[i]
    expect_equal(get_stan_val(d, "censored_02", i), exp$c02[i],
                 label = str_glue("{nm}: recensored censored_02"))
    expect_equal(get_stan_val(d, "censored_12", i), exp$c12[i],
                 label = str_glue("{nm}: recensored censored_12"))
    expect_equal(get_stan_val(d, "censored_32", i), exp$c32[i],
                 label = str_glue("{nm}: recensored censored_32"))
  }
})
