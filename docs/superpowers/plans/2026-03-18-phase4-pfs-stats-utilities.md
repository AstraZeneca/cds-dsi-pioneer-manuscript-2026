# Plan: Phase 4 — PFS Statistical Utilities

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/pfs.stanfunctions`
**Target assertions:** ~60 new assertions

## Background

`pfs.stanfunctions` contains several pure-computation statistical functions that have no RNG
and are easily oracle-tested. This phase covers the survival quantile estimators,
Kaplan-Meier utilities, and the marginal exit probability calculator.

## Functions to cover

| Function | Signature |
|----------|-----------|
| `calculate_log_marginal_exit_prob` | `row_vector → row_vector` |
| `survival_quantiles` | `(array[] int, int, vector p) → (vector, array[] int)` |
| `survival_median` | `(array[] int, int) → (real, int)` |
| `km_quantiles` | `(vector km_survival, vector p) → (vector, array[] int)` |
| `km_median` | `(vector km_survival) → (real, int)` |
| `calc_km_pfs_n` | `(vector km_survival, data real n) → real` |

`estimate_kaplan_meier` is large and used indirectly — skipped for this phase.

---

## Task 1 — R helper `helper-pfs-stats.R`

**File:** `tests/testthat/helper-pfs-stats.R`

```r
r_calculate_log_marginal_exit_prob <- function(log_cond_prob_surv) {
  T <- length(log_cond_prob_surv)
  out <- numeric(T)
  for (t in seq_len(T)) {
    out[t] <- log1p(-exp(log_cond_prob_surv[t]))  # log(1 - exp(x)) = log1m_exp(x)
    if (t > 1) out[t] <- out[t] + sum(log_cond_prob_surv[seq_len(t - 1)])
  }
  out
}

r_survival_quantiles <- function(surv_time, last_surv_time, p) {
  N <- length(surv_time)
  sorted <- sort(surv_time)
  quantiles <- numeric(length(p))
  cannot_calc <- integer(length(p))
  for (j in seq_along(p)) {
    k <- 1L
    while (p[j] >= k * (1.0 / N)) k <- k + 1L
    if (sorted[max(k - 1L, 1L)] > last_surv_time) {
      quantiles[j] <- 0.0; cannot_calc[j] <- 1L
    } else {
      pos <- p[j] * (N - 1) + 1
      d <- pos - (k - 1)
      quantiles[j] <- sorted[max(k - 1L, 1L)] + d * (sorted[k] - sorted[max(k - 1L, 1L)])
      cannot_calc[j] <- 0L
    }
  }
  list(quantiles = quantiles, cannot_calculate = cannot_calc)
}

r_km_quantiles <- function(km_survival, p) {
  T <- length(km_survival)
  P <- length(p)
  quantiles <- rep(0.0, P)
  cannot_calc <- rep(1L, P)
  p_idx <- order(p, decreasing = TRUE)  # descending
  curr <- 1L
  # Handle p > km_survival[1]
  while (curr <= P && p[p_idx[curr]] > km_survival[1]) {
    quantiles[p_idx[curr]] <- 0.0; cannot_calc[p_idx[curr]] <- 0L; curr <- curr + 1L
  }
  for (t in 2:T) {
    while (curr <= P &&
           km_survival[t] <= p[p_idx[curr]] &&
           km_survival[t - 1] > p[p_idx[curr]]) {
      idx <- p_idx[curr]
      weight <- (km_survival[t - 1] - p[idx]) / (km_survival[t - 1] - km_survival[t])
      quantiles[idx] <- (t - 2) + weight  # 0-based
      cannot_calc[idx] <- 0L
      curr <- curr + 1L
    }
    if (curr > P) break
  }
  list(quantiles = quantiles, cannot_calculate = cannot_calc)
}

r_calc_km_pfs_n <- function(km_survival, n) {
  max_t <- length(km_survival) - 1L
  if (n <= 0) return(km_survival[1])
  if (n >= max_t) return(km_survival[max_t + 1L])
  t_floor <- floor(n)
  weight <- n - t_floor
  (1 - weight) * km_survival[t_floor + 1L] + weight * km_survival[t_floor + 2L]
}
```

**Test file:** `tests/testthat/test-helper-pfs-stats.R`
**Assertions (~20):**
- `r_calculate_log_marginal_exit_prob`: single interval (= log(1-exp(lp))); multi-interval sum
- `r_survival_quantiles`: exact median for odd-count vector; `cannot_calculate=1` when
  quantile exceeds `last_surv_time`
- `r_km_quantiles`: flat survival curve; step-function curve; beyond range → `cannot_calc=1`
- `r_calc_km_pfs_n`: at 0 = km[1]; at max_t = km[end]; interpolation midpoint

---

## Task 2 — Stan harness `test_pfs_stats_all.stan`

**File:** `tests/testthat/stan/test_pfs_stats_all.stan`

```stan
functions {
  #include pfs.stanfunctions
  #include util.stanfunctions
  #include pos.stanfunctions
}
data {
  // calculate_log_marginal_exit_prob
  int<lower=1> T_lmep;
  row_vector[T_lmep] log_cond_surv;
  // survival_quantiles
  int<lower=1> N_sq;
  array[N_sq] int surv_time;
  int last_surv_time;
  int<lower=1> P_sq;
  vector[P_sq] sq_p;
  // km_quantiles
  int<lower=1> T_km;
  vector[T_km] km_survival;
  int<lower=1> P_km;
  vector[P_km] km_p;
  // calc_km_pfs_n
  data real pfs_n_query;
}
generated quantities {
  row_vector[T_lmep] lmep = calculate_log_marginal_exit_prob(log_cond_surv);
  vector[P_sq] sq_quantiles;
  array[P_sq] int sq_cannot_calc;
  (sq_quantiles, sq_cannot_calc) = survival_quantiles(surv_time, last_surv_time, sq_p);
  real surv_med;
  int surv_med_flag;
  (surv_med, surv_med_flag) = survival_median(surv_time, last_surv_time);
  vector[P_km] km_q;
  array[P_km] int km_cannot_calc;
  (km_q, km_cannot_calc) = km_quantiles(km_survival, km_p);
  real km_med;
  int km_med_flag;
  (km_med, km_med_flag) = km_median(km_survival);
  real pfs_n = calc_km_pfs_n(km_survival, pfs_n_query);
}
```

**Test file:** `tests/testthat/test-stan-pfs-stats.R`

```r
library(testthat)

test_that("calculate_log_marginal_exit_prob matches R oracle", {
  log_cond <- c(-0.1, -0.2, -0.3)
  stan_data <- list(
    T_lmep = 3L, log_cond_surv = log_cond,
    N_sq = 5L, surv_time = c(2L, 4L, 6L, 8L, 10L),
    last_surv_time = 12L, P_sq = 2L, sq_p = c(0.25, 0.5),
    T_km = 6L, km_survival = c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1),
    P_km = 2L, km_p = c(0.25, 0.5),
    pfs_n_query = 2.5
  )
  fit <- test_stan_function("test_pfs_stats_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  # --- log marginal exit prob ---
  expected <- r_calculate_log_marginal_exit_prob(log_cond)
  for (t in seq_along(log_cond)) {
    expect_equal(get_stan_val(d, "lmep", t), expected[t], tolerance = 1e-6,
      label = paste0("lmep[", t, "]"))
  }

  # --- survival_quantiles ---
  r_sq <- r_survival_quantiles(c(2L, 4L, 6L, 8L, 10L), 12L, c(0.25, 0.5))
  expect_equal(get_stan_val(d, "sq_quantiles", 1), r_sq$quantiles[1], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "sq_quantiles", 2), r_sq$quantiles[2], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "sq_cannot_calc", 1), r_sq$cannot_calculate[1])
  expect_equal(get_stan_val(d, "sq_cannot_calc", 2), r_sq$cannot_calculate[2])

  # --- survival_median (delegates to survival_quantiles at p=0.5) ---
  expect_equal(get_stan_val(d, "surv_med"), r_sq$quantiles[2], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "surv_med_flag"), r_sq$cannot_calculate[2])

  # --- km_quantiles ---
  km_surv <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1)
  r_km <- r_km_quantiles(km_surv, c(0.25, 0.5))
  expect_equal(get_stan_val(d, "km_q", 1), r_km$quantiles[1], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "km_q", 2), r_km$quantiles[2], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "km_cannot_calc", 1), r_km$cannot_calculate[1])
  expect_equal(get_stan_val(d, "km_cannot_calc", 2), r_km$cannot_calculate[2])

  # --- km_median (delegates to km_quantiles at p=0.5) ---
  expect_equal(get_stan_val(d, "km_med"), r_km$quantiles[2], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "km_med_flag"), r_km$cannot_calculate[2])

  # --- calc_km_pfs_n with linear interpolation at 2.5 ---
  expected_pfs_n <- r_calc_km_pfs_n(km_surv, 2.5)
  expect_equal(get_stan_val(d, "pfs_n"), expected_pfs_n, tolerance = 1e-6)
})

test_that("survival_quantiles: cannot_calculate=1 when beyond last_surv_time", {
  # All 5 patients survived beyond last_surv_time=3, median is not calculable
  stan_data <- list(
    T_lmep = 2L, log_cond_surv = c(-0.1, -0.2),
    N_sq = 5L, surv_time = c(5L, 6L, 7L, 8L, 9L),
    last_surv_time = 3L, P_sq = 1L, sq_p = c(0.5),
    T_km = 3L, km_survival = c(1.0, 0.9, 0.8),
    P_km = 1L, km_p = c(0.5),
    pfs_n_query = 1.0
  )
  fit <- test_stan_function("test_pfs_stats_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())
  expect_equal(get_stan_val(d, "sq_cannot_calc", 1), 1L)
})
```

**Assertion count:** ~25 Stan assertions

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-pfs-stats.R` + `test-helper-pfs-stats.R` | ~20 R |
| 2 | `test_pfs_stats_all.stan` + `test-stan-pfs-stats.R` | ~25 Stan |
| **Total** | | **~45 new assertions** |
