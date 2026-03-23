# Test for visit_only_survival_time_rng in Stan
#
# Exercises both the unconditional and conditional (forecast-from-obs_time) overloads.
# Verifies:
#   - Events can only occur at visit weeks (never between visits)
#   - Near-zero survival → event at first (eligible) visit
#   - Near-one survival → censored at max_t
#   - Conditional overload with right_censored=0 passes through observed event
#   - Conditional overload skips visits <= obs_time

library(testthat)
library(cmdstanr)

test_that("visit_only_survival_time_rng: events only at visit weeks, correct censoring", {
  stan_file_path <- here::here(
    "tests", "testthat", "stan",
    "test_visit_only_survival_time_rng.stan"
  )

  N_DRAWS <- 500L

  mod <- cmdstan_model(stan_file_path, include_paths = here::here("stan"),
                       force_recompile = TRUE, quiet = FALSE)
  fit <- mod$sample(
    data = list(N_DRAWS = N_DRAWS),
    seed = 12345,
    chains = 1,
    iter_sampling = 1,
    iter_warmup = 0,
    fixed_param = TRUE,
    show_messages = TRUE,
    refresh = 0
  )

  draws_df <- posterior::as_draws_df(fit$draws())

  # Helper to extract array[N_DRAWS] from draws
  get_vec <- function(var) {
    vapply(seq_len(N_DRAWS), function(j) {
      as.integer(draws_df[[sprintf("%s[%d]", var, j)]][1])
    }, integer(1))
  }

  # ── Case 2: Near-zero survival → event at first visit (week 10) ──────
  case2_time <- get_vec("case2_time")
  case2_cens <- get_vec("case2_cens")
  expect_true(all(case2_time == 10L),
              label = "Case 2: near-zero survival → event at first visit week 10")
  expect_true(all(case2_cens == 0L),
              label = "Case 2: all events (not censored)")

  # ── Case 3: Near-one survival → censored at max_t ────────────────────
  case3_time <- get_vec("case3_time")
  case3_cens <- get_vec("case3_cens")
  expect_true(all(case3_time == 100L),
              label = "Case 3: near-one survival → censored at max_t=100")
  expect_true(all(case3_cens == 1L),
              label = "Case 3: all censored")

  # ── Case 4: Conditional, observed event passed through ────────────────
  case4_time <- get_vec("case4_time")
  case4_cens <- get_vec("case4_cens")
  expect_true(all(case4_time == 15L),
              label = "Case 4: observed event at week 15 passed through")
  expect_true(all(case4_cens == 0L),
              label = "Case 4: all uncensored (observed event)")

  # ── Case 5: Conditional, forecast from obs_time=12 ───────────────────
  # Visits at {6,12,18,24,30}, obs_time=12 → skip 6,12 → first eligible is 18
  # Near-zero survival → event at week 18
  case5_time <- get_vec("case5_time")
  case5_cens <- get_vec("case5_cens")
  expect_true(all(case5_time == 18L),
              label = "Case 5: conditional near-zero → event at first future visit (week 18)")
  expect_true(all(case5_cens == 0L),
              label = "Case 5: all events")

  # ── Case 6: Conditional, high survival → censored ────────────────────
  case6_time <- get_vec("case6_time")
  case6_cens <- get_vec("case6_cens")
  expect_true(all(case6_time == 100L),
              label = "Case 6: conditional near-one → censored at max_t")
  expect_true(all(case6_cens == 1L),
              label = "Case 6: all censored")

  # ── Case 7: Events only at visit weeks (moderate hazard) ─────────────
  # Visits at {8, 16, 24, 32}. All event times must be in this set or max_t.
  case7_time <- get_vec("case7_time")
  case7_cens <- get_vec("case7_cens")
  valid_weeks <- c(8L, 16L, 24L, 32L, 100L)
  expect_true(all(case7_time %in% valid_weeks),
              label = "Case 7: all event times are visit weeks or max_t")
  # Events should be uncensored, max_t should be censored
  event_idx <- case7_cens == 0L
  cens_idx  <- case7_cens == 1L
  if (any(event_idx)) {
    expect_true(all(case7_time[event_idx] %in% c(8L, 16L, 24L, 32L)),
                label = "Case 7: events only at visit weeks")
  }
  if (any(cens_idx)) {
    expect_true(all(case7_time[cens_idx] == 100L),
                label = "Case 7: censored draws at max_t")
  }
  # With ~30% per-visit hazard and 4 visits, most draws should have events
  expect_gt(sum(event_idx), N_DRAWS * 0.5,
            label = "Case 7: majority of draws have events (moderate hazard)")
})
