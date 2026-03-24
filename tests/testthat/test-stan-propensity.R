library(testthat)
library(here)

# 5 patients: 3 trial (group 1 at level 1), 2 RWD (group 2).
# Two hierarchy levels: (trial, patient). Patient level is identity.
make_propensity_data <- function(
    n_patients     = 5L,
    n_trial        = 3L,
    n_covar        = 2L,
    enable         = 1L,
    split_level    = 1L,
    target_group   = 1L,
    beta_intercept = 0.0,
    beta           = NULL) {
  n_rwd <- n_patients - n_trial
  # Level 1: trial group (1 = trial, 2 = RWD)
  # Level 2: patient identity (1..n_patients)
  trial_group  <- c(rep(1L, n_trial), rep(2L, n_rwd))
  patient_id   <- seq_len(n_patients)
  level_groups <- cbind(trial_group, patient_id)

  # Simple covariate matrix: each patient gets a unique row
  X <- matrix(seq_len(n_patients * n_covar) / 10, nrow = n_patients, ncol = n_covar)

  list(
    n_patients                  = n_patients,
    n_levels                    = 2L,
    patient_level_groups        = level_groups,
    n_covar                     = n_covar,
    covar_design_matrix         = X,
    enable_propensity_weighting = enable,
    propensity_split_level      = split_level,
    propensity_target_group     = target_group,
    propensity_intercept_sd     = 1.0,
    propensity_coef_sd          = 0.5,
    beta_propensity_intercept_1 = beta_intercept,
    beta_propensity             = if (is.null(beta)) rep(0.0, n_covar) else beta
  )
}

test_stan_propensity <- function(data) {
  test_stan_function(
    here("tests", "testthat", "stan", "test_propensity_weights.stan"),
    data
  )
}

test_that("propensity transformed_data: target range and indicator are correct", {
  data <- make_propensity_data()
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  expect_equal(get_stan_val(d, "out_target_start"), 1)
  expect_equal(get_stan_val(d, "out_target_end"),   3)
  expect_equal(get_stan_val(d, "out_n_target"),     3)

  # Trial patients (1-3) are target; RWD patients (4-5) are not
  for (i in 1:3) expect_equal(get_stan_val(d, "out_is_target", i), 1,
    label = paste0("out_is_target[", i, "]"))
  for (i in 4:5) expect_equal(get_stan_val(d, "out_is_target", i), 0,
    label = paste0("out_is_target[", i, "]"))
})

test_that("propensity weights: target patients always get weight 1.0", {
  data <- make_propensity_data(beta_intercept = 2.0, beta = c(0.5, -0.5))
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  # Target patients (1-3): weight must be exactly 1.0 regardless of betas
  for (i in 1:3) expect_equal(get_stan_val(d, "out_likelihood_weight", i), 1.0,
    tolerance = 1e-10, label = paste0("weight[", i, "]"))
})

test_that("propensity weights: non-target patients get inv_logit score", {
  beta_intercept <- 0.5
  beta           <- c(0.3, -0.2)
  data <- make_propensity_data(beta_intercept = beta_intercept, beta = beta)
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  X <- data$covar_design_matrix

  # RWD patients (4, 5): weight = inv_logit(intercept + X[i,] %*% beta)
  for (i in 4:5) {
    expected <- plogis(beta_intercept + sum(X[i, ] * beta))
    expect_equal(get_stan_val(d, "out_likelihood_weight", i), expected,
      tolerance = 1e-6, label = paste0("weight[", i, "]"))
  }
})

test_that("propensity weights: disabled mode gives all-ones weights", {
  data <- make_propensity_data(enable = 0L, beta_intercept = 99.0, beta = c(5, 5))
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  for (i in 1:5) expect_equal(get_stan_val(d, "out_likelihood_weight", i), 1.0,
    tolerance = 1e-10, label = paste0("weight[", i, "] when disabled"))
})

test_that("propensity transformed_data: correct when RWD patients come first", {
  # RWD patients first, trial patients last — contiguity still holds
  data <- make_propensity_data()
  # Flip patient ordering: RWD (group 2) first, trial (group 1) last
  data$patient_level_groups <- rbind(
    cbind(rep(2L, 2), 1:2),   # patients 1-2: RWD
    cbind(rep(1L, 3), 3:5)    # patients 3-5: trial
  )
  data$propensity_target_group <- 1L
  fit <- test_stan_propensity(data)
  d   <- posterior::as_draws_df(fit$draws())

  expect_equal(get_stan_val(d, "out_target_start"), 3)
  expect_equal(get_stan_val(d, "out_target_end"),   5)
  expect_equal(get_stan_val(d, "out_n_target"),     3)
  for (i in 1:2) expect_equal(get_stan_val(d, "out_is_target", i), 0,
    label = paste0("out_is_target[", i, "]"))
  for (i in 3:5) expect_equal(get_stan_val(d, "out_is_target", i), 1,
    label = paste0("out_is_target[", i, "]"))
})
