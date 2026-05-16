// modules/propensity/parameters.stan
// Propensity logistic regression coefficients.
// Sized 0 when weighting is disabled to avoid sampling unused parameters.
// Explicit intercept needed because covar_design_matrix is centered (not QR).
array[enable_propensity_weighting ? 1 : 0] real beta_propensity_intercept;
vector[enable_propensity_weighting ? n_covar : 0] beta_propensity;
