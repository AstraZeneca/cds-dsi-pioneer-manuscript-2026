vector[n_train_trials] trial_log_total_rate_effect = zeros_vector(n_train_trials);
vector[n_train_patients] patient_log_total_rate_effect = zeros_vector(n_train_patients);
vector[n_train_patients] patient_log_total_rate = rep_vector(pop_log_total_rate, n_train_patients);
// Fraction logit (population intercept only initially)
vector[n_train_patients] patient_decrease_frac_logit = rep_vector(pop_decrease_frac_logit, n_train_patients);
  
// Population coefficients are sampled in QR space (theta). Recover original-scale (beta) coefficients.
// R_covar_design_matrix is upper-triangular R from X = Q * R factorization of standardized design matrix.
// Solve R * beta = theta  => beta = R^{-1} * theta using triangular backsolve.
vector[n_covar] pop_decrease_frac_logit_coef = mdivide_right_tri_low(pop_decrease_frac_logit_coef_qr', R_covar_design_matrix')';
vector[n_covar] pop_decrease_prop_logis_coef = mdivide_right_tri_low(pop_decrease_prop_logis_coef_qr', R_covar_design_matrix')';

// Linear predictor uses theta directly (Q * theta) for numerical stability.
vector[n_train_patients] patient_decrease_frac_logit_linpred = Q_covar_design_matrix * pop_decrease_frac_logit_coef_qr;
  
// Trial- and patient-level hierarchical covariate effects (original beta scale then projected to QR space)
// matrix[n_train_trials, n_covar] trial_decrease_frac_logit_coef = rep_matrix(0, n_train_trials, n_covar);
matrix[add_patient_level_frac ? n_train_patients : 0, n_covar] patient_decrease_frac_logit_coef; 

if (!pop_covar_coef_only) {
  if (add_patient_level_frac) {
    patient_decrease_frac_logit_coef = raw_patient_decrease_frac_logit_coef .* rep_matrix(patient_decrease_frac_logit_coef_sd, n_train_patients);
    matrix[n_train_patients, n_covar] patient_decrease_frac_logit_coef_qr = patient_decrease_frac_logit_coef * R_covar_design_matrix';
    patient_decrease_frac_logit_linpred += rows_dot_product(Q_covar_design_matrix, patient_decrease_frac_logit_coef_qr);
  }
}

// Hierarchical random effects
if (!pop_rates_param_only) {
  if (add_trial_level_total_rate) {              // old flag add_trial_level_net_rate
    trial_log_total_rate_effect = trial_log_total_rate_sd * raw_trial_log_total_rate;
  }
  patient_log_total_rate_effect = patient_log_total_rate_sd * raw_patient_log_total_rate;
}

// Apply linear predictor & random effects
// old: patient_log_total_rate += patient_log_net_rate_linpred + trial_log_net_rate_effect[...] + patient_log_net_rate_effect;
// Removed total rate random effects from fraction logit to prevent double counting scale effects.
// Fraction logit should only reflect relative allocation (mix) independent of total rate.
patient_decrease_frac_logit += patient_decrease_frac_logit_linpred; // (add dedicated frac REs here later if needed)
if (add_patient_level_frac) {
  // Add patient-level intercept random effect on fraction logit (separate from covariate slopes)
  patient_decrease_frac_logit += patient_decrease_frac_logit_sd * raw_patient_decrease_frac_logit;
}

// Recompute population fraction pieces for each patient
vector[n_train_patients] patient_log_decrease_frac = log_inv_logit(patient_decrease_frac_logit);
vector[n_train_patients] patient_log_growth_frac   = log1m_inv_logit(patient_decrease_frac_logit);

// Update total rate only with random effects (no covariate shift)
patient_log_total_rate += trial_log_total_rate_effect[train_patient_trial] + patient_log_total_rate_effect;

// Derive patient-specific growth/decrease rates from total and fractions
vector[n_train_patients] patient_log_decrease_rate = patient_log_total_rate + patient_log_decrease_frac;
vector[n_train_patients] patient_log_growth_rate   = patient_log_total_rate + patient_log_growth_frac;
  
vector[n_train_patients] patient_log_growth_lag_effect = zeros_vector(n_train_patients);
vector[n_train_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_train_patients);
  
if (!pop_growth_lag_param_only) {
  patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
  patient_log_growth_lag += patient_log_growth_lag_effect;
}
  
vector[n_train_patients] patient_decrease_prop_logis_linpred = Q_covar_design_matrix * pop_decrease_prop_logis_coef_qr;
vector[n_trials] trial_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_trials);
vector[n_trials] trial_decrease_prop_logis_effect = zeros_vector(n_trials);
vector[n_train_patients] patient_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_train_patients) + patient_decrease_prop_logis_linpred; 
vector[n_train_patients] patient_decrease_prop_logis_effect = zeros_vector(n_train_patients);
  
if (!pop_initial_states_param_only ) {
  if (add_trial_level_prop) {
    trial_decrease_prop_logis_effect = trial_decrease_prop_logis_sd * raw_trial_decrease_prop_logis;
    trial_decrease_prop_logis += trial_decrease_prop_logis_effect;
  }
 
  patient_decrease_prop_logis_effect = patient_decrease_prop_logis_sd * raw_patient_decrease_prop_logis; 
  patient_decrease_prop_logis += trial_decrease_prop_logis_effect[patient_trial[train_patients_pos:train_patients_end]] 
                              + patient_decrease_prop_logis_effect;
} 
  
vector[n_train_patients] patient_log_decrease_prop = log_inv_logit(patient_decrease_prop_logis);
vector[n_train_patients] patient_log_growth_prop = log1m_inv_logit(patient_decrease_prop_logis);
  
vector[n_train_patients] patient_tumor_gp_rho = independ_long_process_noise ? zeros_vector(n_train_patients) : rep_vector(exp(log_pop_tumor_gp_rho), n_train_patients);  
  
matrix[n_total_train_visits, 2] states; 
  
profile("states") {
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_train_patients] log_patient_tumor_gp_rho_effect; 
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
  
  if (!independ_long_process_noise && !pop_rho_param_only) {
    log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
    log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
    patient_tumor_gp_rho = 
      exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial[train_patients_pos:train_patients_end]] + log_patient_tumor_gp_rho_effect);
  }

  states = calc_states(
    train_patient_visit_pos,
    train_patient_visits,
    patient_tumor_gp_rho,
    delta,
    pop_process_sd,
    independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
    // raw_patient_process_noise,
    independ_long_process_noise, independ_cross_process_noise,
    append_col(patient_log_decrease_prop, patient_log_growth_prop),
    exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
    rep_vector(0.0001, n_train_patients), // exp(patient_log_growth_lag), 
    0.0001, // exp(pop_log_growth_transition_rate),
    n_shards, // && !debug,
    0 // debug 
  );
}
