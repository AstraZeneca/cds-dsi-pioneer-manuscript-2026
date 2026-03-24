// modules/propensity/hyperparams.stan
// Prior SDs for propensity logistic regression.
// With p standardized covariates, target logit-SD ≈ 2 (avoids U-shaped scores).
// Recommended defaults (set from R): intercept_sd=1.0, coef_sd=2/sqrt(n_covar).
real<lower=0> propensity_intercept_sd;
real<lower=0> propensity_coef_sd;
