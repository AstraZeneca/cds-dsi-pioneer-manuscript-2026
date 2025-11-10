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

// Time-varying tumor burden covariate: log(SLD) at each time point
// states_full_grid[1] = log(regression), states_full_grid[2] = log(growth)
// log(SLD) = log(regression + growth) = log_sum_exp(log_regression, log_growth)
array[n_causes] matrix[n_patients, max_all_t] oe_time_varying_log_hazard_ratio = rep_array(rep_matrix(0, n_patients, max_all_t), n_causes);

// Compute time-varying tumor burden from states if enabled
if (oe_enable_pop_tumor_cov) {
  for (k in 1:n_causes) {
    real tumor_coef = oe_tumor_coef_pop[k][1];
    
    for (i in 1:n_patients) {
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, i);
      
      // Calculate offset: where does absolute time 1 map to in this patient's states grid?
      // states_full_grid[*, i, 1] corresponds to patient's first visit
      // Absolute time 1 → column index = 1 - first_visit + 1 = 2 - first_visit
      int first_visit = t_patient_visits[visit_start];
      int states_start_col = max(1, 2 - first_visit);  // Start of absolute time range [1, max_all_t]
      int states_end_col = max_all_t - first_visit + 1;  // End of absolute time range
      
      // Extract log(SLD) for absolute times [1, max_all_t] from this patient's states grid
      // states_full_grid gives log(normalized_SLD) where normalized = ratio to baseline
      // Convert to absolute SLD: log(SLD_absolute) = log(baseline) + log(normalized)
      row_vector[max_all_t] log_sld_normalized = 
        log_sum_exp(states_full_grid[1][i, states_start_col:states_end_col], 
                    states_full_grid[2][i, states_start_col:states_end_col]);
      
      // Add baseline to get absolute SLD in cm
      row_vector[max_all_t] log_sld_absolute = log_baseline_sld[i] + log_sld_normalized;
      
      // Z-score normalize using distribution of ALL observed SLD values
      // This makes coefficient interpretable as log HR per 1-SD change in log(SLD in cm)
      row_vector[max_all_t] log_sld_z = (log_sld_absolute - mean_log_sld_all) / sd_log_sld_all;
      
      // Apply tumor coefficient to get time-varying log hazard ratio
      oe_time_varying_log_hazard_ratio[k, i] = tumor_coef * log_sld_z;
    }
  }
}

// Compute time-INVARIANT non-tumor covariate effects if enabled
if (oe_enable_pop_cov || oe_enable_trial_cov) {
  for (k in 1:n_causes) {
    vector[n_patients] covar_linpred = rep_vector(0, n_patients);
    
    // Population-level non-tumor covariate effects (QR space)
    if (oe_enable_pop_cov) {
      covar_linpred = Q_covar_design_matrix * oe_covar_coef_qr_pop[k];
    }
    
    // Trial-level random non-tumor slopes (additive)
    if (oe_enable_trial_cov) {
      matrix[n_trials, n_covar] trial_slope_qr = oe_raw_trial_slope[k] .* rep_matrix(oe_sd_trial_slope[k], n_trials);
      covar_linpred += rows_dot_product(Q_covar_design_matrix, 
                                        trial_slope_qr[patient_trial]);
    }
    
    // Total log hazard ratio (only non-tumor covariates, no tumor effects here)
    oe_time_invariant_log_hazard_ratio[k] = covar_linpred;
  }
}

// --- Combined: log_cond_prob_surv ---
array[n_causes] matrix<upper=0>[n_patients, max_all_t] log_cond_prob_surv = oe_time_varying_log_hazard_ratio;

for (s in 1:n_trials) {
  int patient_start, patient_end; 
  (patient_start, patient_end) = get_pos(trial_patient_pos, s);

  for (k in 1:n_causes) {
    // Start with baseline hazard (population + trial GP)
    log_cond_prob_surv[k, patient_start:patient_end] += rep_matrix(log_trial_lambda[k, s], get_pos_size(trial_patient_pos, s));
    
    // Add time-invariant covariate effects (broadcast to all times)
    log_cond_prob_surv[k, patient_start:patient_end] += rep_matrix(oe_time_invariant_log_hazard_ratio[k, patient_start:patient_end], max_all_t);
    
    // Add time-varying tumor burden effect
      // log_cond_prob_surv[k, patient_start:patient_end] += oe_time_varying_log_hazard_ratio[k, patient_start:patient_end];
    
    // Transform to log conditional survival probability
    log_cond_prob_surv[k, patient_start:patient_end] = - exp(log_cond_prob_surv[k, patient_start:patient_end]); 
  }
}
