// modules/propensity/transformed_parameters.stan
// Compute per-patient likelihood weights as capped density ratios.
//
// The density ratio P(X|trial)/P(X|RWD) measures pure covariate similarity,
// free of the sample-size base rate that makes raw P(trial|X) uniformly low
// when N_rwd >> N_trial.
//
// Math: log(P(X|trial)/P(X|RWD)) = logit(P(trial|X)) - logit(P(trial))
//     = (beta_0 + X*beta) - log(N_target/N_nontarget)
//
// Capping at 1.0 ensures no RWD patient contributes more than a trial patient.
//
// NOTE: uses covar_design_matrix (original centered/scaled matrix) NOT Q_covar_design_matrix.
// All other modules use the QR-decomposed Q matrix for numerical stability, but propensity
// logistic regression is convex and doesn't need QR. Using the original matrix keeps
// beta_propensity coefficients interpretable (log-odds per SD change per covariate).
vector<lower=0, upper=1>[n_patients] likelihood_weight = ones_vector(n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
  // Log density ratio for all patients
  vector[n_patients] log_density_ratio =
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
    - propensity_log_marginal_odds;

  // Cap at 1.0: min(1, exp(log_dr)) = exp(min(0, log_dr))
  for (i in 1:n_patients) {
    likelihood_weight[i] = exp(fmin(0.0, log_density_ratio[i]));
  }

  // Override target patients to weight 1.0 (range assignment, no loop)
  likelihood_weight[propensity_target_start:propensity_target_end] =
    ones_vector(propensity_n_target);
}
