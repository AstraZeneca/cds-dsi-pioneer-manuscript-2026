// modules/propensity/priors.stan
// Priors scaled to avoid U-shaped prior predictive: coef_sd = 2/sqrt(n_covar)
if (enable_propensity_weighting && propensity_split_level > 0) {
  beta_propensity_intercept ~ normal(0, propensity_intercept_sd);
  beta_propensity ~ normal(0, propensity_coef_sd);
}
