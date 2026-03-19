# Plan: Phase 5 — Piecewise Constant Hazard Likelihood

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/pfs.stanfunctions` (lines 152–422)
**Target assertions:** ~70 new assertions

## Background

`calc_pch_loglik` is the core likelihood function for the piecewise constant hazard survival model.
It handles right-censoring, interval-censoring, multiple exit types, and windowed start/end times.
All 7 overloads ultimately delegate to the full 8-argument version. `pch_lpmf` is a thin
`sum(calc_pch_loglik(...))` wrapper with 5 overloads — these require only scalar-output checks.

Testing `calc_pch_loglik` thoroughly unlocks coverage of the main model likelihood path.

## Functions to cover

| Function | Overloads |
|----------|-----------|
| `calc_pch_loglik` | 7 (full + 6 convenience) |
| `pch_lpmf` | 5 |

---

## Task 1 — R helper `helper-pch.R`

**File:** `tests/testthat/helper-pch.R`

The full 8-arg version is the oracle. All other overloads are tested by ensuring their Stan output
matches the full version with default arguments substituted.

```r
# Full oracle: mirrors the 8-arg Stan calc_pch_loglik exactly
# log_cond_prob_surv: list of matrices, each [n_patients, max_t]
r_calc_pch_loglik <- function(
  last_surv_week, exit_event, right_censored, interval_censored,
  ignore_interval_censoring, log_cond_prob_surv_list,
  start_from, end_at
) {
  n_exit_types <- length(log_cond_prob_surv_list)
  n_patients   <- length(last_surv_week)
  lp <- numeric(n_patients)

  for (i in seq_len(n_patients)) {
    interval_pos <- max(0L, start_from[i])
    interval_end <- min(end_at[i], last_surv_week[i])

    eff_rc <- as.integer(
      right_censored[i] ||
      (end_at[i] < last_surv_week[i] + (1 - ignore_interval_censoring) * interval_censored[i] + 1)
    )
    curr_ic <- if (ignore_interval_censoring || eff_rc) 0L else interval_censored[i]

    # Survival part
    if (interval_pos <= interval_end) {
      for (k in seq_len(n_exit_types)) {
        lp[i] <- lp[i] + sum(log_cond_prob_surv_list[[k]][i, seq(interval_pos + 1L, interval_end + 1L)])
      }
    }

    # Interval censoring mixture
    valid <- (interval_pos <= interval_end) ||
             (interval_end + curr_ic + 1 >= interval_pos)
    if (valid) {
      ie  <- max(interval_end, interval_pos - 1L)
      ic2 <- max(0L, curr_ic - (ie - interval_end))
      mix <- numeric(ic2 + 1L)
      for (c in 0:ic2) {
        if (c > 0) {
          for (k in seq_len(n_exit_types)) {
            mix[c + 1L] <- mix[c + 1L] +
              sum(log_cond_prob_surv_list[[k]][i, seq(ie + 2L, ie + c + 1L)])
          }
        }
        if (!eff_rc) {
          mix[c + 1L] <- mix[c + 1L] +
            log1p(-exp(log_cond_prob_surv_list[[exit_event[i]]][i, ie + c + 2L]))
        }
      }
      lp[i] <- lp[i] + if (ic2 > 0) matrixStats::logSumExp(mix) - log(ic2 + 1) else mix[1]
    }
  }
  lp
}

# Convenience: single exit type, no interval censoring, start=1, end=max_t
r_calc_pch_loglik_simple <- function(last_surv_week, right_censored, log_cond_prob_surv) {
  # log_cond_prob_surv: matrix [max_t, n_patients] (col = patient)
  n_patients <- length(last_surv_week)
  max_t      <- nrow(log_cond_prob_surv)
  r_calc_pch_loglik(
    last_surv_week,
    rep(1L, n_patients),
    right_censored,
    rep(0L, n_patients),
    0L,
    list(t(log_cond_prob_surv)),  # transpose: rows = patients
    rep(1L, n_patients),
    rep(max_t, n_patients)
  )
}
```

**Test file:** `tests/testthat/test-helper-pch.R`
**Assertions (~15):**
- Right-censored patient: lp = sum of log_cond_surv up to last_surv_week
- Event patient (no IC): lp = survival part + log(1 - exp(log_cond_surv[event_week]))
- Interval-censored patient: lp is mixture over possible event weeks
- `start_from > 1`: only weeks from start_from onward contribute

---

## Task 2 — Stan harness `test_pch_loglik_all.stan`

**File:** `tests/testthat/stan/test_pch_loglik_all.stan`

```stan
functions {
  #include pfs.stanfunctions
  #include util.stanfunctions
  #include pos.stanfunctions
}
data {
  int<lower=1> N;         // n_patients
  int<lower=1> MAX_T;     // max time horizon
  int<lower=1> N_EXIT;    // number of exit types (1 or 2)
  array[N] int last_surv_week;
  array[N] int exit_event;
  array[N] int right_censored;
  array[N] int interval_censored;
  int ignore_interval_censoring;
  // log_cond_prob_surv: N_EXIT matrices each [N x MAX_T], flattened row-major
  // Passed as array[N_EXIT] matrix[N, MAX_T]
  array[N_EXIT] matrix[N, MAX_T] log_cond_prob_surv;
  array[N] int start_from;
  array[N] int end_at;
}
generated quantities {
  // Full 8-arg version
  vector[N] lp_full = calc_pch_loglik(
    last_surv_week, exit_event, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv, start_from, end_at
  );
  // 7-arg: default end_at = max_all_t
  vector[N] lp_no_end = calc_pch_loglik(
    last_surv_week, exit_event, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv, start_from
  );
  // Single-type convenience: matrix overload
  vector[N] lp_single = calc_pch_loglik(
    last_surv_week, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv[1]
  );
  // pch_lpmf (scalar sum)
  real lpmf_val = pch_lpmf(
    last_surv_week | right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv[1]
  );
}
```

**Test file:** `tests/testthat/test-stan-pch-loglik.R`

```r
library(testthat)

# Helper: build a simple uniform log conditional survival matrix
make_log_cond_surv <- function(n_patients, max_t, log_surv_per_week = -0.1) {
  matrix(log_surv_per_week, nrow = n_patients, ncol = max_t)
}

test_that("calc_pch_loglik: right-censored patients", {
  # 3 patients, max_t=10, all right-censored at various times
  n <- 3L; max_t <- 10L
  lcs <- make_log_cond_surv(n, max_t, -0.05)
  last_surv <- c(3L, 7L, 10L)
  rc <- c(1L, 1L, 1L)
  ic <- c(0L, 0L, 0L)

  stan_data <- list(
    N = n, MAX_T = max_t, N_EXIT = 1L,
    last_surv_week = last_surv, exit_event = rep(1L, n),
    right_censored = rc, interval_censored = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv = array(lcs, dim = c(1, n, max_t)),
    start_from = rep(1L, n), end_at = rep(max_t, n)
  )
  fit <- test_stan_function("test_pch_loglik_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_pch_loglik_simple(last_surv, rc, t(lcs))
  for (i in seq_len(n)) {
    expect_equal(get_stan_val(d, "lp_full", i), expected[i], tolerance = 1e-6,
      label = paste0("lp_full[", i, "] right-censored"))
  }
})

test_that("calc_pch_loglik: event patients (no censoring)", {
  n <- 2L; max_t <- 8L
  lcs <- make_log_cond_surv(n, max_t, -0.1)
  last_surv <- c(3L, 6L)  # event at week 3 and 6
  rc <- c(0L, 0L)
  ic <- c(0L, 0L)

  stan_data <- list(
    N = n, MAX_T = max_t, N_EXIT = 1L,
    last_surv_week = last_surv, exit_event = rep(1L, n),
    right_censored = rc, interval_censored = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv = array(lcs, dim = c(1, n, max_t)),
    start_from = rep(1L, n), end_at = rep(max_t, n)
  )
  fit <- test_stan_function("test_pch_loglik_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_pch_loglik_simple(last_surv, rc, t(lcs))
  for (i in seq_len(n)) {
    expect_equal(get_stan_val(d, "lp_full", i), expected[i], tolerance = 1e-6,
      label = paste0("lp_full[", i, "] event"))
  }
})

test_that("calc_pch_loglik: interval-censored patients", {
  n <- 2L; max_t <- 10L
  lcs <- make_log_cond_surv(n, max_t, -0.05)
  last_surv <- c(4L, 6L)
  rc <- c(0L, 0L)
  ic <- c(2L, 1L)  # interval of 2 and 1 week respectively

  stan_data <- list(
    N = n, MAX_T = max_t, N_EXIT = 1L,
    last_surv_week = last_surv, exit_event = rep(1L, n),
    right_censored = rc, interval_censored = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv = array(lcs, dim = c(1, n, max_t)),
    start_from = rep(1L, n), end_at = rep(max_t, n)
  )
  fit <- test_stan_function("test_pch_loglik_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_pch_loglik(
    last_surv, rep(1L, n), rc, ic, 0L, list(lcs), rep(1L, n), rep(max_t, n)
  )
  for (i in seq_len(n)) {
    expect_equal(get_stan_val(d, "lp_full", i), expected[i], tolerance = 1e-5,
      label = paste0("lp_full[", i, "] interval-censored"))
  }
})

test_that("pch_lpmf equals sum(calc_pch_loglik)", {
  n <- 3L; max_t <- 6L
  lcs <- make_log_cond_surv(n, max_t, -0.08)
  last_surv <- c(2L, 4L, 6L); rc <- c(1L, 0L, 1L); ic <- c(0L, 0L, 0L)

  stan_data <- list(
    N = n, MAX_T = max_t, N_EXIT = 1L,
    last_surv_week = last_surv, exit_event = rep(1L, n),
    right_censored = rc, interval_censored = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv = array(lcs, dim = c(1, n, max_t)),
    start_from = rep(1L, n), end_at = rep(max_t, n)
  )
  fit <- test_stan_function("test_pch_loglik_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  expected_sum <- sum(r_calc_pch_loglik_simple(last_surv, rc, t(lcs)))
  expect_equal(get_stan_val(d, "lpmf_val"), expected_sum, tolerance = 1e-5)
})
```

**Assertion count:** ~28 Stan assertions

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-pch.R` + `test-helper-pch.R` | ~15 R |
| 2 | `test_pch_loglik_all.stan` + `test-stan-pch-loglik.R` | ~28 Stan |
| **Total** | | **~43 new assertions** |
