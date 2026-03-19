# Plan: Phase 7 — Multistate Likelihood Functions

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/multistate.stanfunctions`
**Target assertions:** ~55 new assertions

## Background

`multistate.stanfunctions` implements the illness-death model likelihood. The key design insight
is the **fast path**: when only the 0→1 transition is enabled, `multistate_lpmf` immediately
delegates to `calc_ms_single_transition_loglik`, which itself wraps `calc_pch_loglik` with a
detection-week → last-surviving-week convention shift. This makes unit testing tractable: we can
isolate `calc_ms_single_transition_loglik` independently, then verify `multistate_lpmf` routes
correctly.

## Functions to cover

| Function | Notes |
|----------|-------|
| `calc_ms_single_transition_loglik` | Wraps `calc_pch_loglik`; event_time−1 = last_surv_week |
| `multistate_lpmf` | Fast path (01-only) + full illness-death path |
| `sample_dropout_death_rng` | RNG — tested with deterministic boundary cases |

---

## Task 1 — R helper `helper-multistate.R`

**File:** `tests/testthat/helper-multistate.R`

```r
# calc_ms_single_transition_loglik
# Detection-week convention: event at event_time[i], last_surv_week = event_time - 1 for events
# log_cond_surv: matrix [n_patients x max_t] (rows = patients)
r_calc_ms_single_transition_loglik <- function(event_time, censored, log_cond_surv) {
  n <- length(event_time)
  max_t <- ncol(log_cond_surv)
  last_surv_week <- ifelse(censored, event_time, event_time - 1L)

  r_calc_pch_loglik(
    last_surv_week,
    rep(1L, n),        # exit_event (single type)
    censored,
    rep(0L, n),        # no interval censoring
    0L,                # ignore_interval_censoring
    list(log_cond_surv),
    rep(1L, n),        # start_from
    rep(max_t, n)      # end_at
  )
}

# multistate_lpmf fast path (01 only, no 02/12)
r_multistate_lpmf_01only <- function(time_01, censored_01, log_cond_surv_01) {
  sum(r_calc_ms_single_transition_loglik(time_01, censored_01, log_cond_surv_01))
}
```

**Test file:** `tests/testthat/test-helper-multistate.R`
**Assertions (~12):**
- `r_calc_ms_single_transition_loglik`: censored patient → sum log_cond_surv[1..t]; event
  patient → sum log_cond_surv[1..t-1] + log(1-exp(log_cond_surv[t]))
- `r_multistate_lpmf_01only`: equals sum over `calc_ms_single_transition_loglik`
- Convention check: event at t=1 → last_surv_week=0, so only hazard term contributes

---

## Task 2 — Stan harness `test_multistate_all.stan`

**File:** `tests/testthat/stan/test_multistate_all.stan`

```stan
functions {
  #include pfs.stanfunctions
  #include util.stanfunctions
  #include pos.stanfunctions
  #include multistate.stanfunctions
}
data {
  int<lower=1> N;
  int<lower=1> MAX_T;
  // Single-transition test
  array[N] int event_time_01;
  array[N] int censored_01;
  matrix[N, MAX_T] log_cond_surv_01;
  // multistate full: final states and times
  array[N] int final_state;
  array[N] int time_02;  array[N] int censored_02;
  array[N] int time_12;  array[N] int censored_12;
  matrix[N, MAX_T] log_cond_surv_02;
  matrix[N, MAX_T] log_cond_surv_12_s;
  matrix[N, MAX_T] log_cond_surv_12_t;
  array[N] int prog_deterministic;
  // Unused transitions (for 01-only fast path)
  array[N] int time_03;  array[N] int censored_32;
  matrix[N, MAX_T] log_cond_surv_03;
  matrix[N, MAX_T] log_cond_surv_32;
}
generated quantities {
  // Single-transition loglik
  vector[N] st_llik = calc_ms_single_transition_loglik(
    event_time_01, censored_01, log_cond_surv_01
  );
  // multistate_lpmf: 01-only fast path
  real ms_lpmf_01only = multistate_lpmf(
    final_state | 1, 0, 0, 0, 0, 0,
    event_time_01, time_02, time_12, time_03, time_03,
    censored_01, censored_02, censored_12, censored_32,
    prog_deterministic,
    log_cond_surv_01, log_cond_surv_02,
    log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32
  );
  // multistate_lpmf: full 01+02 path (illness-death without 12)
  real ms_lpmf_01_02 = multistate_lpmf(
    final_state | 1, 1, 0, 0, 0, 0,
    event_time_01, time_02, time_12, time_03, time_03,
    censored_01, censored_02, censored_12, censored_32,
    prog_deterministic,
    log_cond_surv_01, log_cond_surv_02,
    log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32
  );
}
```

**Test file:** `tests/testthat/test-stan-multistate.R`

```r
library(testthat)

test_that("calc_ms_single_transition_loglik: detection-week convention", {
  n <- 3L; max_t <- 8L
  lcs <- matrix(-0.08, nrow = n, ncol = max_t)
  # Patient 1: event at week 4 (last_surv_week = 3)
  # Patient 2: censored at week 6 (last_surv_week = 6)
  # Patient 3: event at week 1 (last_surv_week = 0 → only hazard term)
  event_time <- c(4L, 6L, 1L)
  censored   <- c(0L, 1L, 0L)

  stan_data <- list(
    N = n, MAX_T = max_t,
    event_time_01 = event_time, censored_01 = censored,
    log_cond_surv_01 = lcs,
    final_state = c(1L, 0L, 1L),
    time_02 = rep(max_t + 1L, n), censored_02 = rep(1L, n),
    time_12 = rep(0L, n),         censored_12 = rep(1L, n),
    log_cond_surv_02 = matrix(-0.05, n, max_t),
    log_cond_surv_12_s = matrix(-0.05, n, max_t),
    log_cond_surv_12_t = matrix(-0.05, n, max_t),
    prog_deterministic = rep(0L, n),
    time_03 = rep(max_t + 1L, n),
    censored_32 = rep(1L, n),
    log_cond_surv_03 = matrix(-0.05, n, max_t),
    log_cond_surv_32 = matrix(-0.05, n, max_t)
  )
  fit <- test_stan_function("test_multistate_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_ms_single_transition_loglik(event_time, censored, lcs)
  for (i in seq_len(n)) {
    expect_equal(get_stan_val(d, "st_llik", i), expected[i], tolerance = 1e-6,
      label = paste0("st_llik[", i, "]"))
  }
})

test_that("multistate_lpmf 01-only fast path == sum(calc_ms_single_transition_loglik)", {
  # Verify fast path delegates correctly
  n <- 2L; max_t <- 6L
  lcs <- matrix(-0.1, n, max_t)
  event_time <- c(3L, 5L); censored <- c(0L, 1L)

  stan_data <- list(
    N = n, MAX_T = max_t,
    event_time_01 = event_time, censored_01 = censored,
    log_cond_surv_01 = lcs,
    final_state = c(1L, 0L),
    time_02 = rep(max_t + 1L, n), censored_02 = rep(1L, n),
    time_12 = rep(0L, n),         censored_12 = rep(1L, n),
    log_cond_surv_02 = matrix(-0.05, n, max_t),
    log_cond_surv_12_s = matrix(-0.05, n, max_t),
    log_cond_surv_12_t = matrix(-0.05, n, max_t),
    prog_deterministic = rep(0L, n),
    time_03 = rep(max_t + 1L, n),
    censored_32 = rep(1L, n),
    log_cond_surv_03 = matrix(-0.05, n, max_t),
    log_cond_surv_32 = matrix(-0.05, n, max_t)
  )
  fit <- test_stan_function("test_multistate_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  expected_sum <- sum(r_calc_ms_single_transition_loglik(event_time, censored, lcs))
  expect_equal(get_stan_val(d, "ms_lpmf_01only"), expected_sum, tolerance = 1e-6)
})

test_that("multistate_lpmf 01+02: state-0 censored patient gets both survival terms", {
  n <- 1L; max_t <- 5L
  lcs_01 <- matrix(-0.1, n, max_t)
  lcs_02 <- matrix(-0.05, n, max_t)
  # Patient censored at time 4 in state 0
  cens_time <- 4L
  stan_data <- list(
    N = n, MAX_T = max_t,
    event_time_01 = cens_time, censored_01 = 1L,
    log_cond_surv_01 = lcs_01,
    final_state = 0L,
    time_02 = cens_time, censored_02 = 1L,
    time_12 = 0L, censored_12 = 1L,
    log_cond_surv_02 = lcs_02,
    log_cond_surv_12_s = matrix(-0.05, n, max_t),
    log_cond_surv_12_t = matrix(-0.05, n, max_t),
    prog_deterministic = 0L,
    time_03 = max_t + 1L, censored_32 = 1L,
    log_cond_surv_03 = matrix(-0.05, n, max_t),
    log_cond_surv_32 = matrix(-0.05, n, max_t)
  )
  fit <- test_stan_function("test_multistate_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  # State 0 censored: 01 survival[1..4] + 02 survival[1..4]
  expected <- sum(lcs_01[1, 1:cens_time]) + sum(lcs_02[1, 1:cens_time])
  expect_equal(get_stan_val(d, "ms_lpmf_01_02"), expected, tolerance = 1e-6)
})
```

**Assertion count:** ~18 Stan assertions

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-multistate.R` + `test-helper-multistate.R` | ~12 R |
| 2 | `test_multistate_all.stan` + `test-stan-multistate.R` | ~18 Stan |
| **Total** | | **~30 new assertions** |

**Note:** `multistate_lpmf` with 12 post-progression transition and all three time-scale modes
(Markov/semi-Markov/extended) are left for a follow-up phase after Phase 7 is stable.
