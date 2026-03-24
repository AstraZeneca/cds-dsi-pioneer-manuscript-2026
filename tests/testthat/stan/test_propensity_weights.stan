// Test harness for stan/modules/propensity/ transformed_data and
// transformed_parameters logic.
// betas are passed as data to fix their values without sampling.
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}
data {
  // Hierarchy (minimal: n_levels levels, identity mapping at last level)
  int<lower=1> n_patients;
  int<lower=1> n_levels;
  array[n_patients, n_levels] int patient_level_groups;

  int<lower=0> n_covar;
  matrix[n_patients, n_covar] covar_design_matrix;

  // Propensity module flags / hyperparams
  int<lower=0, upper=1> enable_propensity_weighting;
  int<lower=0> propensity_split_level;
  int<lower=0> propensity_target_group;
  real<lower=0> propensity_intercept_sd;
  real<lower=0> propensity_coef_sd;

  // Propensity betas — passed as data so values are fixed for testing
  real beta_propensity_intercept_1;  // will become beta_propensity_intercept[1]
  vector[n_covar] beta_propensity;
}
transformed data {
  // Replicate parameters.stan conditional sizing using data fields
  array[enable_propensity_weighting ? 1 : 0] real beta_propensity_intercept;
  if (enable_propensity_weighting) {
    beta_propensity_intercept[1] = beta_propensity_intercept_1;
  }

  // Module under test: sets propensity_target_start/end/n_target, propensity_is_target
  #include "modules/propensity/transformed_data.stan"
}
generated quantities {
  // Outputs from transformed_data.stan
  int out_target_start           = propensity_target_start;
  int out_target_end             = propensity_target_end;
  int out_n_target               = propensity_n_target;
  array[n_patients] int out_is_target = propensity_is_target;

  // transformed_parameters.stan logic inlined (block includes cannot go in GQ)
  vector[n_patients] out_likelihood_weight = ones_vector(n_patients);
  if (enable_propensity_weighting && propensity_split_level > 0) {
    out_likelihood_weight = inv_logit(
      beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
    );
    out_likelihood_weight[propensity_target_start:propensity_target_end] =
      ones_vector(propensity_n_target);
  }
}
