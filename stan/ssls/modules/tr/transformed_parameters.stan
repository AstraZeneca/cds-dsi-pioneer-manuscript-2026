// tr/transformed_parameters.stan — active linear predictor assembly for total rate module
vector[n_trials] tr_effect_trial_intercept = enable_trial_intercept_tr ? tr_sd_trial_intercept * tr_raw_trial_intercept : rep_vector(0, n_trials);
vector[n_patients] tr_effect_patient_intercept = enable_patient_intercept_tr ? tr_sd_patient_intercept * tr_raw_patient_intercept : rep_vector(0, n_patients);
vector[n_patients] tr_linpred_pop = enable_pop_cov_tr ? (Q_covar_design_matrix * tr_coef_qr_pop) : rep_vector(0, n_patients);

// Trial slope deviations: construct QR-space linear predictor contributions per patient
vector[n_patients] tr_linpred_trial = rep_vector(0, n_patients);

if (enable_trial_cov_tr) {
  matrix[n_trials, n_covar] tr_trial_slope_qr = (tr_raw_trial_slope .* rep_matrix(tr_sd_trial_slope', n_trials));
  tr_linpred_trial = rows_dot_product(Q_covar_design_matrix, tr_trial_slope_qr[patient_trial]);
}

// Patient slope deviations
vector[n_patients] tr_linpred_patient = rep_vector(0, n_patients);

if (enable_patient_cov_tr) {
  matrix[n_patients, n_covar] tr_patient_slope_qr = (tr_raw_patient_slope .* rep_matrix(tr_sd_patient_slope', n_patients));
  tr_linpred_patient = rows_dot_product(Q_covar_design_matrix, tr_patient_slope_qr);
}

vector[n_patients] tr_loc_patient = tr_loc_pop
  + tr_linpred_pop + tr_linpred_trial + tr_linpred_patient
  + tr_effect_trial_intercept[patient_trial]
  + tr_effect_patient_intercept; 

vector[enable_patient_process_noise_tr ? n_patients : 0] tr_log_sd_patient_process_noise; 
vector[enable_patient_process_noise_tr ? n_patients : 0] tr_phi_patient_process_noise; 
matrix[enable_patient_process_noise_tr ? n_patients : 0, max_t_width] tr_patient_process_noise; 

if (enable_patient_process_noise_tr) {
  // Non-centered parameterization: work in log-space, exponentiate once at the end
  // If patient hierarchy is disabled, all patients get population value
  if (enable_patient_process_noise_sd_tr) {
    tr_log_sd_patient_process_noise = tr_log_sd_pop_process_noise[1] + tr_sd_patient_log_sd_process_noise[1] * tr_raw_patient_log_sd_process_noise;
  } else {
    tr_log_sd_patient_process_noise = rep_vector(tr_log_sd_pop_process_noise[1], n_patients);
  }

  // Map to [0,1] via inv_logit link function
  // Population phi is inv_logit(logit_phi_pop), patient deviations on logit scale
  if (enable_patient_process_noise_phi_tr) {
    tr_phi_patient_process_noise = inv_logit(
      tr_logit_phi_pop_process_noise[1] + tr_sd_patient_phi_process_noise[1] * tr_raw_patient_phi_process_noise
    );
  } else {
    tr_phi_patient_process_noise = rep_vector(inv_logit(tr_logit_phi_pop_process_noise[1]), n_patients);
  }

  vector[n_patients] tr_sd_patient_process_noise = exp(tr_log_sd_patient_process_noise);

  // AR(1) process: deviation[t] = phi * deviation[t-1] + sigma * innovation[t]
  // This creates time-varying deviations with mean 0 that revert to baseline
  tr_patient_process_noise[, 1] = tr_raw_patient_process_noise[, 1] .* tr_sd_patient_process_noise;

  for (t in 2:max_t_width) {
    tr_patient_process_noise[, t] = tr_phi_patient_process_noise .* tr_patient_process_noise[, t - 1]
      + tr_raw_patient_process_noise[, t] .* tr_sd_patient_process_noise;
  }
}

// Population-level time-varying process noise (shared AR(1) across all patients)
row_vector[enable_pop_process_noise_tr ? max_t_width : 0] tr_pop_process_noise;

if (enable_pop_process_noise_tr) {
  real sigma_pop = exp(tr_log_sd_pop_process_noise_pop[1]);
  // log(inv_logit(x)) = x - log1p_exp(x), more stable than log(inv_logit(x))
  real log_phi_pop = tr_logit_phi_pop_process_noise_pop[1] - log1p_exp(tr_logit_phi_pop_process_noise_pop[1]);

  // Compute φ^t for t = 1..max_t_width (vectorized via exp/log)
  row_vector[max_t_width] phi_powers = exp(log_phi_pop * linspaced_row_vector(max_t_width, 1, max_t_width));

  // Scale raw innovations
  row_vector[max_t_width] z_pop = sigma_pop * tr_raw_pop_process_noise;

  // Vectorized AR(1): y[t] = Σ_{k=1}^{t} φ^(t-k) * z[k] = φ^t * cumsum(z/φ^k)
  tr_pop_process_noise = cumulative_sum(z_pop ./ phi_powers) .* phi_powers;
}