library(testthat)

test_that("r_calc_pch_loglik_simple: right-censored patient accumulates survival", {
  # Patient survives weeks 1..3, censored at week 3
  # start_from=1, end_at=8 (max_t), interval_pos=1, interval_end=min(8,3)=3
  # eff_rc=1 (right_censored=1), known survival sum = lcs[1, 1:3] = 3 * (-0.1)
  # mix[1]: eff_rc=1 so no hazard term, mix[1]=0 -> lp += 0
  lcs    <- matrix(-0.1, nrow = 1, ncol = 8)
  result <- r_calc_pch_loglik_simple(3L, 1L, lcs)
  expect_equal(result[1], 3 * (-0.1), tolerance = 1e-9)
})

test_that("r_calc_pch_loglik_simple: event at week t accumulates survival[1..t-1] + hazard[t]", {
  # Event at week 4: last_surv_week=4, right_censored=0, interval_censored=0
  # interval_pos=1, interval_end=min(8,4)=4
  # eff_rc=0 (end_at=8 >= 4+0+1=5)
  # Known survival sum: lcs[1, 1:4] = 4 * (-0.1)  (cols 1,2,3,4)
  # After adjustment: ie=max(4, 0)=4, ic2=0
  # mix[1]: c=0, col_ev = 4+0+1 = 5, log1p(-exp(-0.1))
  # But wait: known survival already included col 4, so lp = 4*(-0.1) + log1p(-exp(-0.1))
  # However the canonical interpretation is survival weeks 1..3 then hazard at 4.
  # Let's verify: Stan sums cols interval_pos:interval_end = 1:4 (four terms).
  # The hazard is at col interval_end+1 = 5. So lp = 4*(-0.1) + log1p(-exp(-0.1)).
  lcs      <- matrix(-0.1, nrow = 1, ncol = 8)
  result   <- r_calc_pch_loglik_simple(4L, 0L, lcs)
  expected <- 4 * (-0.1) + log1p(-exp(-0.1))
  expect_equal(result[1], expected, tolerance = 1e-9)
})

test_that("r_calc_pch_loglik_simple: event at week 1", {
  # last_surv_week=1, right_censored=0, interval_censored=0, start_from=1, end_at=5
  # interval_pos=1, interval_end=min(5,1)=1
  # eff_rc=0 (end_at=5 >= 1+0+1=2)
  # Known survival sum: lcs[1, 1:1] = -0.2
  # After adjustment: ie=max(1, 0)=1, ic2=0
  # mix[1]: c=0, col_ev=1+0+1=2, log1p(-exp(-0.2))
  # lp = -0.2 + log1p(-exp(-0.2))
  lcs      <- matrix(-0.2, nrow = 1, ncol = 5)
  result   <- r_calc_pch_loglik_simple(1L, 0L, lcs)
  expected <- -0.2 + log1p(-exp(-0.2))
  expect_equal(result[1], expected, tolerance = 1e-9)
})

test_that("r_calc_pch_loglik_simple: multi-patient vectorized", {
  lcs       <- matrix(-0.1, nrow = 3, ncol = 6)
  last_surv <- c(2L, 4L, 6L)
  rc        <- c(1L, 0L, 1L)
  result    <- r_calc_pch_loglik_simple(last_surv, rc, lcs)
  expect_length(result, 3L)
  # Patient 1: right-censored at 2, known survival sum = lcs[1, 1:2] = 2*(-0.1)
  # mix[1]: eff_rc=1, no hazard term, 0. lp = 2*(-0.1)
  expect_equal(result[1], 2 * (-0.1), tolerance = 1e-9)
  # Patient 2: event at 4, known survival sum = lcs[2, 1:4] = 4*(-0.1)
  # hazard at col 5: log1p(-exp(-0.1)). lp = 4*(-0.1) + log1p(-exp(-0.1))
  expect_equal(result[2], 4 * (-0.1) + log1p(-exp(-0.1)), tolerance = 1e-9)
  # Patient 3: right-censored at 6 (max_t=6), known survival sum = lcs[3, 1:6] = 6*(-0.1)
  # mix[1]: eff_rc=1, no hazard. lp = 6*(-0.1)
  expect_equal(result[3], 6 * (-0.1), tolerance = 1e-9)
})

test_that("r_calc_pch_loglik: two exit types, event in type 2", {
  # With two exit types and an event patient, both matrices contribute to survival sum
  # but only the exit_event matrix contributes to the hazard term
  n      <- 1L
  max_t  <- 5L
  lcs1   <- matrix(-0.1, nrow = n, ncol = max_t)  # type 1: progression
  lcs2   <- matrix(-0.2, nrow = n, ncol = max_t)  # type 2: death
  # Patient: last_surv_week=3, exit_event=2 (death), not censored
  result <- r_calc_pch_loglik(
    last_surv_week          = 3L,
    exit_event              = 2L,
    right_censored          = 0L,
    interval_censored       = 0L,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv_list = list(lcs1, lcs2),
    start_from              = 1L,
    end_at                  = max_t
  )
  # Known survival sum over weeks 1:3 for both types:
  # lcs1[1, 1:3] + lcs2[1, 1:3] = 3*(-0.1) + 3*(-0.2) = -0.9
  # After adjustment: ie=3, ic2=0
  # mix[1]: col_ev=4, exit_event=2: log1p(-exp(lcs2[1,4])) = log1p(-exp(-0.2))
  expected <- (3 * (-0.1) + 3 * (-0.2)) + log1p(-exp(-0.2))
  expect_equal(result[1], expected, tolerance = 1e-9)
})

test_that("r_calc_pch_loglik: ignore_interval_censoring=1 zeroes out IC", {
  n     <- 1L
  max_t <- 8L
  lcs   <- matrix(-0.1, nrow = n, ncol = max_t)
  # Patient with interval_censored=3 but ignore_interval_censoring=1
  # Should behave identically to interval_censored=0
  result_ic  <- r_calc_pch_loglik(
    last_surv_week          = 5L,
    exit_event              = 1L,
    right_censored          = 0L,
    interval_censored       = 3L,
    ignore_interval_censoring = 1L,
    log_cond_prob_surv_list = list(lcs),
    start_from              = 1L,
    end_at                  = max_t
  )
  result_no_ic <- r_calc_pch_loglik(
    last_surv_week          = 5L,
    exit_event              = 1L,
    right_censored          = 0L,
    interval_censored       = 0L,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv_list = list(lcs),
    start_from              = 1L,
    end_at                  = max_t
  )
  expect_equal(result_ic[1], result_no_ic[1], tolerance = 1e-12)
})
