test_that("add_lfo_fields adds required LFO fields", {
  mock_stan_data <- list(
    n_patients = 10L,
    calendar_day = seq(100L, 190L, by = 10L),
    t_patient_visits_day = rep(c(7L, 42L, 84L), 10L),
    n_patient_visits = rep(3L, 10L)
  )

  cutoffs <- c(150L, 180L, 210L)
  result <- add_lfo_fields(mock_stan_data, cutoffs)

  expect_equal(result$n_cutoffs, 3L)
  expect_equal(result$cutoff_calendar_day, cutoffs)
  expect_equal(result$max_n_rows, 1L)
  expect_equal(result$max_forecast_horizon, 2L)

  # Original data preserved

  expect_equal(result$n_patients, 10L)
  expect_equal(result$calendar_day, mock_stan_data$calendar_day)
})

test_that("add_lfo_fields respects PSIS mode dimensions", {
  mock_stan_data <- list(
    n_patients = 5L,
    calendar_day = 1:5,
    t_patient_visits_day = 1:15,
    n_patient_visits = rep(3L, 5L)
  )

  cutoffs <- c(100L, 200L, 300L)
  result <- add_lfo_fields(mock_stan_data, cutoffs,
                           max_n_rows = 3L,
                           max_forecast_horizon = 3L)

  expect_equal(result$max_n_rows, 3L)
  expect_equal(result$max_forecast_horizon, 3L)
})

test_that("add_lfo_fields errors without calendar_day", {
  mock_stan_data <- list(n_patients = 5L)
  expect_error(add_lfo_fields(mock_stan_data, c(100L)))
})
