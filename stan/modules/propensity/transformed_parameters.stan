// modules/propensity/transformed_parameters.stan
// Compute per-patient likelihood weights.
// Vectorized: matrix-vector multiply for all patients; range assignment for targets.
//
// NOTE: uses covar_design_matrix (original centered/scaled matrix) NOT Q_covar_design_matrix.
// All other modules use the QR-decomposed Q matrix for numerical stability, but propensity
// logistic regression is convex and doesn't need QR. Using the original matrix keeps
// beta_propensity coefficients interpretable (log-odds per SD change per covariate).
vector<lower=0, upper=1>[n_patients] likelihood_weight = ones_vector(n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
  likelihood_weight = inv_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
  // Override target patients to weight 1.0 (range assignment, no loop)
  likelihood_weight[propensity_target_start:propensity_target_end] =
    ones_vector(propensity_n_target);
}
