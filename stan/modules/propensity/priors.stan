// modules/propensity/priors.stan
// Joint propensity submodel: P(target | covariates).
// The bernoulli likelihood identifies beta_propensity from the data source split.
// The outcome likelihoods provide feedback via likelihood_weight.
if (enable_propensity_weighting && propensity_split_level > 0) {
  // Priors scaled to avoid U-shaped prior predictive: coef_sd = 2/sqrt(n_covar)
  beta_propensity_intercept ~ normal(0, propensity_intercept_sd);
  beta_propensity ~ normal(0, propensity_coef_sd);

  // Propensity likelihood
  propensity_is_target ~ bernoulli_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
}
