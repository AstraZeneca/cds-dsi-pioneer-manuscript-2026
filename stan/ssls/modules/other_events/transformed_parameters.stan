// ============================================================================
// Other Events Model Transformed Parameters
// ============================================================================

// --- Baseline Hazard: Population-level GP ---
array[n_causes] row_vector[max_all_t] log_pop_lambda;

for (k in 1:n_causes) {
  log_pop_lambda[k] = calc_gp_pred(
    all_tumor_measure_t, 
    log_lambda_gp_pop_intercept[k], 
    log_lambda_gp_pop_alpha[k], 
    log_lambda_gp_pop_rho[k], 
    delta, 
    log_lambda_gp_pop_eta[k]
  );
}

// --- Baseline Hazard: Trial-level GP (additive to population) ---
array[n_causes] matrix[oe_enable_trial_baseline_hazard ? n_trials : 0, max_all_t] log_trial_lambda_residual; 
array[n_causes] matrix[n_trials, max_all_t] log_trial_lambda; 

for (k in 1:n_causes) {
  log_trial_lambda[k] = rep_matrix(log_pop_lambda[k], n_trials); 
}

array[n_causes] vector[oe_enable_trial_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_intercept;

if (oe_enable_trial_baseline_hazard) {
  for (k in 1:n_causes) {
    log_lambda_gp_trial_intercept[k] = raw_log_lambda_gp_trial_intercept[k] * log_lambda_gp_trial_intercept_sd[k];
    for (s in 1:n_trials) {
      log_trial_lambda_residual[k, s] = calc_gp_pred(
        all_tumor_measure_t, 
        log_lambda_gp_trial_intercept[k, s], 
        log_lambda_gp_trial_alpha[k], 
        log_lambda_gp_trial_rho[k], 
        delta, 
        log_lambda_gp_trial_eta[k, s]
      );
      log_trial_lambda[k, s] += log_trial_lambda_residual[k, s];
    }
  }
}

// --- Proportional Hazard: Covariate Effects ---
array[n_causes] vector[n_patients] oe_time_invariant_log_hazard_ratio = rep_array(rep_vector(0, n_patients), n_causes);

// Only compute if any covariate flag is enabled
if (oe_enable_pop_cov || oe_enable_pop_tumor_cov || 
    oe_enable_trial_cov || oe_enable_trial_tumor_cov) {
  for (k in 1:n_causes) {
    vector[n_patients] tumor_linpred = rep_vector(0, n_patients);
    vector[n_patients] covar_linpred = rep_vector(0, n_patients);
    
    // Population-level tumor effects (QR space) - SCAFFOLDED, not used yet
    if (oe_enable_pop_tumor_cov) {
      tumor_linpred = Q_tumor_sum_covar * oe_tumor_coef_qr_pop[k];
    }
    
    // Population-level non-tumor covariate effects (QR space)
    if (oe_enable_pop_cov) {
      covar_linpred = Q_covar_design_matrix * oe_covar_coef_qr_pop[k];
    }
    
    // Trial-level random tumor slopes (additive) - SCAFFOLDED, not used yet
    if (oe_enable_trial_tumor_cov) {
      matrix[n_trials, n_tumor_covar] trial_tumor_slope_qr = oe_raw_trial_tumor_slope[k] .* rep_matrix(oe_sd_trial_tumor_slope[k], n_trials);
      tumor_linpred += rows_dot_product(Q_tumor_sum_covar, 
                                        trial_tumor_slope_qr[patient_trial]);
    }
    
    // Trial-level random non-tumor slopes (additive)
    if (oe_enable_trial_cov) {
      matrix[n_trials, n_covar] trial_slope_qr = oe_raw_trial_slope[k] .* rep_matrix(oe_sd_trial_slope[k], n_trials);
      covar_linpred += rows_dot_product(Q_covar_design_matrix, 
                                        trial_slope_qr[patient_trial]);
    }
    
    // Total log hazard ratio
    oe_time_invariant_log_hazard_ratio[k] = tumor_linpred + covar_linpred;
  }
}

// --- Combined: log_cond_prob_surv ---
array[n_causes] matrix<upper=0>[n_patients, max_all_t] log_cond_prob_surv;

for (s in 1:n_trials) {
  int patient_start, patient_end; 
  (patient_start, patient_end) = get_pos(trial_patient_pos, s);

  for (k in 1:n_causes) {
    // Start with baseline hazard (population + trial GP)
    log_cond_prob_surv[k, patient_start:patient_end] = rep_matrix(log_trial_lambda[k, s], get_pos_size(trial_patient_pos, s));
    
    // Add proportional hazard from covariates (always include, will be zeros if disabled)
    log_cond_prob_surv[k, patient_start:patient_end] += rep_matrix(oe_time_invariant_log_hazard_ratio[k, patient_start:patient_end], max_all_t);
    
    // Transform to log conditional survival probability
    log_cond_prob_surv[k, patient_start:patient_end] = - exp(log_cond_prob_surv[k, patient_start:patient_end]); 
  }
}
