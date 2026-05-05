library(testthat)
library(cmdstanr)
library(posterior)

test_that("find_first_week and find_first_forecast_week: Stan matches R reference", {

  # ---- find_first_week cases ----
  ffw_cases <- list(
    list(
      name           = "found_in_obs_visits",
      arr            = c(3L, 3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L, 20L),
      max_all_t      = 52L
    ),
    list(
      name           = "found_in_forecast",
      arr            = c(3L, 3L, 3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L, 20L),
      max_all_t      = 52L
    ),
    list(
      name           = "not_found",
      arr            = c(3L, 3L, 3L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L, 20L),
      max_all_t      = 52L
    ),
    list(
      name           = "min_run_length_2",
      arr            = c(3L, 4L, 4L, 3L),
      values         = c(4L),
      min_run_length = 2L,
      curr_visits    = c(4L, 8L, 12L, 16L),
      forecast_time  = c(20L),
      max_all_t      = 52L
    ),
    list(
      name           = "found_at_first_position",
      arr            = c(4L, 3L, 3L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L),
      max_all_t      = 52L
    ),
    list(
      name           = "multiple_values_in_what",
      arr            = c(3L, 1L, 3L),
      values         = c(1L, 2L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L),
      max_all_t      = 52L
    ),
    list(
      name           = "empty_curr_visits_match_in_forecast",
      arr            = c(3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = integer(0),
      forecast_time  = c(8L, 16L),
      max_all_t      = 52L
    )
  )

  # ---- find_first_forecast_week cases ----
  fffc_cases <- list(
    list(
      name           = "found_at_position_2",
      arr            = c(3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      forecast_time  = c(8L, 16L),
      max_all_t      = 52L
    ),
    list(
      name           = "found_at_first_position",
      arr            = c(4L),
      values         = c(4L),
      min_run_length = 1L,
      forecast_time  = c(8L),
      max_all_t      = 52L
    ),
    list(
      name           = "not_found",
      arr            = c(3L, 3L),
      values         = c(4L),
      min_run_length = 1L,
      forecast_time  = c(8L, 16L),
      max_all_t      = 52L
    )
  )

  # ---- Compute R expected values ----
  ffw_expected <- lapply(ffw_cases, function(cc) {
    r_find_first_week(
      cc$arr, cc$values, cc$min_run_length,
      cc$curr_visits, cc$forecast_time, cc$max_all_t
    )
  })
  fffc_expected <- lapply(fffc_cases, function(cc) {
    r_find_first_forecast_week(
      cc$arr, cc$values, cc$min_run_length,
      cc$forecast_time, cc$max_all_t
    )
  })

  # ---- Build Stan data ----
  N_CASES    <- length(ffw_cases)
  N_FC_CASES <- length(fffc_cases)

  MAX_ARR  <- max(sapply(ffw_cases,  function(x) length(x$arr)))
  MAX_VALS <- max(sapply(ffw_cases,  function(x) length(x$values)))
  MAX_OBS  <- max(1L, max(sapply(ffw_cases,  function(x) length(x$curr_visits))))
  MAX_FORE <- max(sapply(ffw_cases,  function(x) length(x$forecast_time)))

  MAX_ARR_FC  <- max(sapply(fffc_cases, function(x) length(x$arr)))
  MAX_VALS_FC <- max(sapply(fffc_cases, function(x) length(x$values)))
  MAX_FORE_FC <- max(sapply(fffc_cases, function(x) length(x$forecast_time)))

  arr_mat  <- matrix(0L, N_CASES, MAX_ARR)
  vals_mat <- matrix(0L, N_CASES, MAX_VALS)
  obs_mat  <- matrix(0L, N_CASES, MAX_OBS)
  fore_mat <- matrix(0L, N_CASES, MAX_FORE)
  n_arr    <- integer(N_CASES)
  n_vals   <- integer(N_CASES)
  n_obs    <- integer(N_CASES)
  n_fore   <- integer(N_CASES)
  min_run  <- integer(N_CASES)
  max_t    <- integer(N_CASES)

  for (i in seq_along(ffw_cases)) {
    cc         <- ffw_cases[[i]]
    n_arr[i]   <- length(cc$arr)
    n_vals[i]  <- length(cc$values)
    n_obs[i]   <- length(cc$curr_visits)
    n_fore[i]  <- length(cc$forecast_time)
    min_run[i] <- cc$min_run_length
    max_t[i]   <- cc$max_all_t
    arr_mat[i,  seq_len(n_arr[i])]  <- cc$arr
    vals_mat[i, seq_len(n_vals[i])] <- cc$values
    if (n_obs[i] > 0) obs_mat[i, seq_len(n_obs[i])] <- cc$curr_visits
    fore_mat[i, seq_len(n_fore[i])] <- cc$forecast_time
  }

  arr_fc_mat  <- matrix(0L, N_FC_CASES, MAX_ARR_FC)
  vals_fc_mat <- matrix(0L, N_FC_CASES, MAX_VALS_FC)
  fore_fc_mat <- matrix(0L, N_FC_CASES, MAX_FORE_FC)
  n_arr_fc    <- integer(N_FC_CASES)
  n_vals_fc   <- integer(N_FC_CASES)
  n_fore_fc   <- integer(N_FC_CASES)
  min_run_fc  <- integer(N_FC_CASES)
  max_t_fc    <- integer(N_FC_CASES)

  for (i in seq_along(fffc_cases)) {
    cc              <- fffc_cases[[i]]
    n_arr_fc[i]     <- length(cc$arr)
    n_vals_fc[i]    <- length(cc$values)
    n_fore_fc[i]    <- length(cc$forecast_time)
    min_run_fc[i]   <- cc$min_run_length
    max_t_fc[i]     <- cc$max_all_t
    arr_fc_mat[i,   seq_len(n_arr_fc[i])]  <- cc$arr
    vals_fc_mat[i,  seq_len(n_vals_fc[i])] <- cc$values
    fore_fc_mat[i,  seq_len(n_fore_fc[i])] <- cc$forecast_time
  }

  stan_data <- list(
    N_CASES    = N_CASES,
    MAX_ARR    = MAX_ARR,
    MAX_VALS   = MAX_VALS,
    MAX_OBS    = MAX_OBS,
    MAX_FORE   = MAX_FORE,
    n_arr      = n_arr,
    n_vals     = n_vals,
    n_obs      = n_obs,
    n_fore     = n_fore,
    min_run_length = min_run,
    max_all_t  = max_t,
    arr        = arr_mat,
    vals       = vals_mat,
    curr_visits     = obs_mat,
    forecast_time   = fore_mat,
    N_FC_CASES      = N_FC_CASES,
    MAX_ARR_FC      = MAX_ARR_FC,
    MAX_VALS_FC     = MAX_VALS_FC,
    MAX_FORE_FC     = MAX_FORE_FC,
    n_arr_fc        = n_arr_fc,
    n_vals_fc       = n_vals_fc,
    n_fore_fc       = n_fore_fc,
    min_run_fc      = min_run_fc,
    max_all_t_fc    = max_t_fc,
    arr_fc          = arr_fc_mat,
    vals_fc         = vals_fc_mat,
    forecast_time_fc = fore_fc_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_pfs_find_first_week_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  # ---- Assert find_first_week ----
  for (i in seq_along(ffw_cases)) {
    cc  <- ffw_cases[[i]]
    exp <- ffw_expected[[i]]
    expect_equal(
      get_val("out_week", i),
      exp$week,
      label = paste0("find_first_week case ", i, " (", cc$name, ") week")
    )
    expect_equal(
      get_val("out_rc", i),
      exp$right_censored,
      label = paste0("find_first_week case ", i, " (", cc$name, ") right_censored")
    )
  }

  # ---- Assert find_first_forecast_week ----
  for (i in seq_along(fffc_cases)) {
    cc  <- fffc_cases[[i]]
    exp <- fffc_expected[[i]]
    expect_equal(
      get_val("out_week_fc", i),
      exp$week,
      label = paste0("find_first_forecast_week case ", i, " (", cc$name, ") week")
    )
    expect_equal(
      get_val("out_rc_fc", i),
      exp$right_censored,
      label = paste0("find_first_forecast_week case ", i, " (", cc$name, ") right_censored")
    )
  }
})
