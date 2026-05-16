// ============================================================================
// STANDALONE PROPENSITY WEIGHTING MODEL
// ============================================================================
// Pure logistic regression for diagnosing propensity-weighted RWD borrowing.
// Runs in seconds (no PSA dynamics, no multistate hazard, no state-space).
//
// Reuses the same modular propensity includes as the full joint model —
// no duplicated logic. The only new code is the generated quantities block
// for diagnostic outputs.
//
// Outputs:
//   likelihood_weight[n_patients]  — capped density-ratio weights (0, 1]
//   log_density_ratio[n_patients]  — uncapped log P(X|trial)/P(X|RWD)
//   propensity_score[n_patients]   — P(trial|X) on probability scale
//   effective_n                    — sum of non-target weights
//   weighted_covar_mean[n_covar]   — weighted mean of each covariate (non-target)
//   target_covar_mean[n_covar]     — mean of each covariate (target patients)

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

data {
  #include "_hierarchy_data.stan"

  #include "modules/propensity/flags.stan"
  #include "modules/propensity/hyperparams.stan"
}

transformed data {
  #include "_hierarchy_transformed_data.stan"
  #include "modules/propensity/transformed_data.stan"
}

parameters {
  #include "modules/propensity/parameters.stan"
}

transformed parameters {
  #include "modules/propensity/transformed_parameters.stan"
}

model {
  #include "modules/propensity/priors.stan"
  #include "modules/propensity/likelihood.stan"
}

generated quantities {
  // Uncapped log density ratio: logit(P(trial|X)) - logit(P(trial))
  // Positive = more trial-like than base rate, negative = less trial-like.
  // The capped version (likelihood_weight) clips positive values to 0 (weight=1).
  vector[n_patients] log_density_ratio_gq = zeros_vector(n_patients);

  // P(trial|X) on probability scale — the raw propensity score
  vector[n_patients] propensity_score = rep_vector(0.5, n_patients);

  // Effective RWD sample size: sum of weights for non-target patients
  real effective_n = 0.0;

  // Covariate balance: weighted mean (non-target) vs unweighted mean (target)
  vector[n_covar] weighted_covar_mean = zeros_vector(n_covar);
  vector[n_covar] target_covar_mean = zeros_vector(n_covar);

  if (enable_propensity_weighting && propensity_split_level > 0) {
    // Log density ratio (uncapped) and propensity score
    log_density_ratio_gq =
      beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
      - propensity_log_marginal_odds;
    propensity_score = inv_logit(
      beta_propensity_intercept[1] + covar_design_matrix * beta_propensity);

    // Effective N and weighted covariate means for non-target patients
    {
      real weight_sum = 0.0;
      int n_target = 0;

      for (i in 1:n_patients) {
        if (propensity_is_target[i] == 0) {
          effective_n += likelihood_weight[i];
          weight_sum += likelihood_weight[i];
          weighted_covar_mean += likelihood_weight[i] * covar_design_matrix[i]';
        } else {
          n_target += 1;
          target_covar_mean += covar_design_matrix[i]';
        }
      }

      if (weight_sum > 0)
        weighted_covar_mean /= weight_sum;
      if (n_target > 0)
        target_covar_mean /= n_target;
    }
  }
}
