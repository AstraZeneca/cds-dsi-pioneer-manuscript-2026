// test_gr_decay_hierarchy.stan
// Backward-compat anchor for the gr_decay (Gompertz growth-rate decay) level-RE
// hierarchy. Proves the NEW frac-style hierarchy (population -> trial-arm ->
// patient) reduces EXACTLY to the OLD pooled population scalar kappa when
//   (A) all level intercept modes are NONE, and
//   (B) RE levels are enabled but their raw draws are exactly 0.
//
// Strategy 1: include the REAL module fragments
//   modules/gr_decay/transformed_data.stan       (flat-index + bucket-position assembly)
//   modules/gr_decay/parameters.stan             (pop intercept + level RE params)
//   modules/gr_decay/transformed_parameters.stan (kappa_i = exp(linear predictor))
// plus the minimal hierarchy scaffolding they depend on. kappa is computed in
// transformed parameters and copied out in generated quantities — deterministic
// given the (data-supplied) flags + (fixed_param) parameter draws.
//
// To force the RE raw draws to exactly 0 in Config B, the R harness passes
// init = 0 for every gr_decay parameter under fixed_param sampling, so the
// transformed-parameters linear predictor receives sd * 0 = 0 contributions.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
}

data {
  // ----- base hierarchy (subset of _hierarchy_data.stan) -----
  int<lower=1> n_trials;
  int<lower=0> n_patients;
  array[n_patients] int<lower=1, upper=n_trials> patient_trial;
  int<lower=1> n_levels;
  array[n_levels] int<lower=1> n_groups_per_level;
  array[n_patients, n_levels] int<lower=1> patient_level_groups;
  int<lower=0> n_covar;
  matrix[n_patients, n_covar] covar_design_matrix;

  // ----- forecast routing -----
  int<lower=0> n_forecast_patients;
  array[n_forecast_patients] int forecast_patient_idx;

  // ----- gr_decay flags (modules/gr_decay/flags.stan) -----
  int<lower=0,upper=1> enable_gr_decay;
  int<lower=0,upper=1> enable_pop_cov_gr_decay;
  array[n_levels] int<lower=0,upper=4> enable_level_intercept_gr_decay;
  array[n_levels] int<lower=0,upper=1> enable_level_cov_gr_decay;

  // ----- gr_decay hyperparams (modules/gr_decay/hyperparams.stan) -----
  // (only the FE intercept SDs are read in transformed parameters; the rest are
  // priors-only but declared so the module fragments resolve every name.)
  real gr_decay_log_loc_pop_mean;
  real<lower=0> gr_decay_log_loc_pop_sd;
  vector[n_covar] gr_decay_coef_qr_pop_mean;
  vector<lower=0>[n_covar] gr_decay_coef_qr_pop_sd;
  array[n_levels] real<lower=0> gr_decay_sd_level_intercept_sd;
  array[n_levels] real<lower=0> gr_decay_fe_sd_level_intercept;
  array[n_levels] row_vector<lower=0>[n_covar] gr_decay_sd_level_slope_sd;
}

transformed data {
  // LEVEL_MODE_* constants + hierarchy validation + n_total_groups/level_pos.
  #include "_hierarchy_transformed_data.stan"
  // n_forecast_groups_per_level + n_forecast_level_pos.
  #include "_forecast_routing_transformed_data.stan"

  // QR design matrix (n_covar = 0 in this test -> empty, but the module reads it).
  #include "_qr_decomposition.stan"

  // Flat indices + bucket positions for the level-indexed log(kappa) hierarchy.
  #include "modules/gr_decay/transformed_data.stan"
}

parameters {
  #include "modules/gr_decay/parameters.stan"
}

transformed parameters {
  #include "modules/gr_decay/transformed_parameters.stan"
}

generated quantities {
  // The per-patient Gompertz decay rate produced by the real module path.
  vector[n_forecast_patients] kappa_out = gr_decay_kappa;
}
