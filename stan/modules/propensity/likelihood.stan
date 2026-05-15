// modules/propensity/likelihood.stan
// Propensity submodel likelihood: identifies beta_propensity from data source split.
// Gated by fit_propensity_data — skipped during prior predictive runs.
if (fit_propensity_data && enable_propensity_weighting && propensity_split_level > 0) {
  propensity_is_target ~ bernoulli_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
}
