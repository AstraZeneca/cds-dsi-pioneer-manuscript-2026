// ============================================================================
// Other Events Model Transformed Parameters
// ============================================================================

// --- Baseline Hazard: Population-level GP ---
array[n_causes] row_vector[max_all_t] log_pop_lambda;

if (n_causes > 0) {
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
}

// --- Baseline Hazard: Trial-level GP (additive to population) ---
array[n_causes] matrix[oe_enable_trial_baseline_hazard ? n_trials : 0, max_all_t] log_trial_lambda_residual; 
array[n_causes] matrix[n_trials, max_all_t] log_trial_lambda; 

if (n_causes > 0) {
  for (k in 1:n_causes) {
    log_trial_lambda[k] = rep_matrix(log_pop_lambda[k], n_trials); 
  }
}

array[n_causes] vector[oe_enable_trial_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_intercept;

if (n_causes > 0 && oe_enable_trial_baseline_hazard) {
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
// Patient rates (patient_log_decrease_rate, patient_log_growth_rate) are used directly
array[n_causes] matrix[n_patients, max_all_t] oe_time_varying_log_hazard_ratio = rep_array(rep_matrix(0, n_patients, max_all_t), n_causes);

// Compute time-varying tumor burden from states if enabled
if (n_causes > 0 && oe_enable_pop_tumor_cov) {
  // Normalization constants (median and IQR) for SLD are computed in transformed_data from observed data
  // This provides fixed, iteration-stable normalization for consistent prior interpretation
  
  // Compute tumor covariate effects using normalized values
  for (k in 1:n_causes) {
    // Tumor covariate coefficients:
    // [1] = log(SLD) effect
    // [2] = log(decrease rate) effect  
    // [3] = log(growth rate) effect
    // [4] = SLD velocity effect (d/dt log(SLD))
    
    for (i in 1:n_patients) {
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, i);
      
      // Patient's first visit time in absolute weeks
      int first_visit = t_patient_visits[visit_start];
      
      // Mapping: absolute time t → patient-relative column = t - first_visit + 1
      // Data invariant: first_visit <= 0, so states_start_col = 2 - first_visit >= 2
      int states_start_col = 2 - first_visit;  // Column for absolute time 1
      int states_end_col = states_start_col + max_all_t - 1;  // Column for absolute time max_all_t
      
      // Extract log(SLD) for absolute times [1, max_all_t] from this patient's states grid
      // states_full_grid gives log(normalized_SLD) where normalized = ratio to baseline
      row_vector[max_all_t] log_sld_normalized = log_sum_exp(
        states_full_grid[1][i, states_start_col:states_end_col],
        states_full_grid[2][i, states_start_col:states_end_col]
      );

      // Add baseline to get absolute SLD in cm
      row_vector[max_all_t] log_sld_absolute = log_baseline_sld[i] + log_sld_normalized;
      
      // Median-center and IQR-scale using constants from observed data (computed in transformed_data)
      row_vector[max_all_t] log_sld_standardized = (log_sld_absolute - median_log_sld_obs) / iqr_log_sld_obs;

      // Initialize with log(SLD) effect
      oe_time_varying_log_hazard_ratio[k, i] = oe_tumor_coef_pop[k][1] * log_sld_standardized;
      
      // Add patient-level rate effects if coefficients are provided
      if (n_tumor_covar >= 3) {
        if (enable_patient_process_noise_tr) {
          // Time-varying rates: extract contiguous slice using same indexing as states
          row_vector[max_all_t] patient_log_decrease_rate_abs = patient_log_decrease_rate[i, states_start_col:states_end_col];
          row_vector[max_all_t] patient_log_growth_rate_abs = patient_log_growth_rate[i, states_start_col:states_end_col];
          oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][2] * patient_log_decrease_rate_abs;
          oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][3] * patient_log_growth_rate_abs;
        } else {
          // Constant rates: use scalar broadcasting (same as main branch behavior)
          // scalar * row_vector broadcasts efficiently without extra allocation
          oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][2] * patient_log_decrease_rate[i, 1];
          oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][3] * patient_log_growth_rate[i, 1];
        }
      }

      // Add SLD velocity (time-varying first derivative)
      if (n_tumor_covar >= 4) {
        // Compute velocity as weekly change in log(SLD)
        // velocity[t] = log_sld[t] - log_sld[t-1]
        row_vector[max_all_t] sld_velocity = rep_row_vector(0, max_all_t);
        
        // Vectorized computation: velocity[2:T] = log_sld[2:T] - log_sld[1:T-1]
        // First time point remains 0 (no prior measurement)
        sld_velocity[2:max_all_t] = log_sld_absolute[2:] - log_sld_absolute[:(max_all_t - 1)];
        
        // Z-score normalize velocity
        // Use statistics from t >= 2 (exclude first point which is zero by construction)
        real mean_velocity = mean(sld_velocity[2:]);
        real sd_velocity = sd(sld_velocity[2:]);
        row_vector[max_all_t] sld_velocity_z = (sld_velocity - mean_velocity) / sd_velocity;
        
        // Add velocity effect to hazard ratio
        oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][4] * sld_velocity_z;
      }
    }
  }
}

// Compute time-INVARIANT non-tumor covariate effects if enabled
if (n_causes > 0 && (oe_enable_pop_cov || oe_enable_trial_cov)) {
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

if (n_causes > 0) {
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
}
