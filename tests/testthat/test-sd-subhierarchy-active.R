# tests/testthat/test-sd-subhierarchy-active.R
# Smoke test: the SD sub-hierarchy active code path executes end-to-end on a
# minimal synthetic input. Asserts only structural facts:
#   - New parameters allocate (non-zero size)
#   - Sampling completes without crashing (5 warmup + 5 sampling iterations)
#   - tr_sd_intercept_pergroup varies across groups at the activated level
# NOT a convergence test. NOT a diagnostic test.

library(testthat)

test_that("active sub-hierarchy allocates parameters and samples without crashing", {
  fixture <- testthat::test_path("fixtures", "sd-hierarchy", "minimal-active-stan-data.rds")
  skip_if(
    !file.exists(fixture),
    "Minimal active-path fixture not yet captured — see README.md"
  )
  skip_on_ci()  # fit takes a few seconds; keep CI lightweight

  library(cmdstanr)

  stan_data <- readRDS(fixture)

  # Precondition: mode matrix should have at least one active entry at (patient, arm).
  # Level positions depend on the fixture; we check by name via patient_level_groups.
  expect_true(
    sum(stan_data$enable_sd_level_intercept_mode_tr) > 0,
    info = "fixture mode matrix should have at least one non-zero entry"
  )

  model <- cmdstan_model(
    testthat::test_path("..", "..", "stan", "psa", "pioneer.stan"),
    include_paths = c(
      testthat::test_path("..", "..", "stan"),
      testthat::test_path("..", "..", "stan", "psa")
    )
  )

  source(testthat::test_path("..", "..", "r", "pioneer", "initializers.R"))

  fit <- model$sample(
    data = stan_data,
    seed = 1,
    iter_warmup = 5,
    iter_sampling = 5,
    chains = 1,
    parallel_chains = 1,
    init = create_pioneer_initializer(stan_data),
    refresh = 0,
    show_messages = FALSE
  )

  # Structural assertions
  draws <- posterior::as_draws_df(fit$draws(
    variables = c(
      "tr_log_sd_level_intercept_pop",
      "tr_cp_log_sd_level_intercept",
      "tr_sd_intercept_pergroup"
    )
  ))

  # 1. Population log-SD is present (size > 0 because at least one level is active).
  pop_cols <- grep("^tr_log_sd_level_intercept_pop", names(draws), value = TRUE)
  expect_gt(length(pop_cols), 0,
            label = "tr_log_sd_level_intercept_pop should be non-empty when sub-hierarchy is active")

  # 2. CP bucket is present if any RE_CP mode is configured.
  # (If only FE or RE modes are active, cp bucket is 0 — we relax the test to cover both cases.)
  has_re_cp <- any(stan_data$enable_sd_level_intercept_mode_tr == 4L)
  if (has_re_cp) {
    cp_cols <- grep("^tr_cp_log_sd_level_intercept", names(draws), value = TRUE)
    expect_gt(length(cp_cols), 0, label = "tr_cp_log_sd_level_intercept should be non-empty")
  }

  # 3. Per-group SD at the active level varies across groups within any draw.
  # The pergroup param is indexed [lv, g]. For each active level, check variation.
  # Find the active level (first level L with has_sd_subhierarchy > 0).
  active_levels <- which(
    rowSums(stan_data$enable_sd_level_intercept_mode_tr) > 0
  )
  expect_gt(length(active_levels), 0)

  # Take the first active level; extract pergroup columns for that level.
  lv <- active_levels[1]
  pg_pattern <- sprintf("^tr_sd_intercept_pergroup\\[%d,", lv)
  pg_cols <- grep(pg_pattern, names(draws), value = TRUE)
  expect_gt(length(pg_cols), 1L, label = "pergroup should have at least 2 groups to vary over")

  # Within a single draw, the per-group SDs should differ.
  row1 <- as.numeric(draws[1, pg_cols])
  expect_gt(
    stats::sd(row1, na.rm = TRUE), 0,
    label = "per-group SD should vary across groups when sub-hierarchy is active"
  )
})
