array[n_causes] matrix[n_base_separate_trials, max_confresp_week] log_crcr_lambda; // log baseline hazard

profile("crcr population baseline hazards") {
  for (s in 1:n_base_separate_trials) {
    for (k in 1:n_causes) {
      log_crcr_lambda[k, s] = 
        calc_gp_pred(confresp_range, log_crcr_lambda_gp_intercept[k, s], log_crcr_lambda_gp_alpha[k, s], log_crcr_lambda_gp_rho[k, s], delta, log_crcr_lambda_gp_eta[k, s]);
    }
  }
}

// Trial level variation in log baseline hazard
array[n_causes] matrix[add_trial_level_baseline_hazard ? n_trials : 0, max_confresp_week] log_crcr_trial_lambda_residual; 
array[n_causes] matrix[n_trials, max_confresp_week] log_crcr_trial_lambda; 

for (k in 1:n_causes) {
  if (separate_baseline_hazard) {
    log_crcr_trial_lambda[k] = log_crcr_lambda[k];
  } else {
    log_crcr_trial_lambda[k] = rep_matrix(log_crcr_lambda[k, 1], n_trials); 
  }
}

array[n_causes] vector[add_trial_level_baseline_hazard ? n_trials : 0] log_crcr_lambda_gp_trial_intercept;
array[n_causes, no_prop_hazard ? 0 : n_trials] vector[n_tumor_covar + n_covar] crcr_covar_trial_coef;

// Trial level variation in prop hazard parameters
array[n_causes, add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_tumor_covar + n_covar] crcr_covar_trial_coef_residual;

if (add_trial_level_baseline_hazard) {
  for (k in 1:n_causes) {
    log_crcr_lambda_gp_trial_intercept[k] = raw_log_crcr_lambda_gp_trial_intercept[k] * log_crcr_lambda_gp_trial_intercept_sd[k];
    
    for (s in 1:n_trials) {
      log_crcr_trial_lambda_residual[k, s] = calc_gp_pred(
        confresp_range,
        log_crcr_lambda_gp_trial_intercept[k, s], log_crcr_lambda_gp_trial_alpha[k], log_crcr_lambda_gp_trial_rho[k], delta, log_crcr_lambda_gp_trial_eta[k, s]
      );
      
      log_crcr_trial_lambda[k, s] += log_crcr_trial_lambda_residual[k, s];
    }
  }
}

// cholesky_factor_cov[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar : 1] L_crcr_covar_trial_cov;
// 
// if (no_prop_hazard || !add_trial_level_prop_hazard) {
//   L_crcr_covar_trial_cov[1, 1] = 1; 
// }

if (!no_prop_hazard) {
  // if (add_trial_level_prop_hazard) {
  //   L_crcr_covar_trial_cov = diag_pre_multiply(crcr_covar_trial_sd, L_crcr_covar_trial_corr);
  // }
  
  for (k in 1:n_causes) {
    for (s in 1:n_trials) {
      if (add_trial_level_prop_hazard) {
        // crcr_covar_trial_coef_residual[k, s] = L_crcr_covar_trial_cov * raw_crcr_covar_trial_coef[k, s];
        crcr_covar_trial_coef_residual[k, s] = crcr_covar_trial_sd .* raw_crcr_covar_trial_coef[k, s];
        crcr_covar_trial_coef[k, s] = append_row(crcr_tumor_stim_pop_coef[k, 1], crcr_covar_effect[k, 1]) + crcr_covar_trial_coef_residual[k, s];
      } else if (separate_prop_hazard) {
        crcr_covar_trial_coef[k, s] = append_row(crcr_tumor_stim_pop_coef[k, s], crcr_covar_effect[k, s]);
      } else {
        crcr_covar_trial_coef[k, s] = append_row(crcr_tumor_stim_pop_coef[k, 1], crcr_covar_effect[k, 1]);
      }
    }
  }
} 

array[n_causes] vector[no_prop_hazard ? 0 : n_patients] patient_log_crcr_hazard_ratio; // Log proportional hazard
array[n_causes] matrix<upper = 0>[n_patients, max_confresp_week] log_crcr_cond_prob_surv; // Log conditional probability of remaining unclassified

profile("log_crcr_cond_prob_surv") { // Calculate patient-interval conditional probability of survival (without a classified response) 
  int patient_pos = 1;
  
  for (s in 1:n_trials) {
    int patient_end = patient_pos + n_trial_patients[s] - 1;
  
    for (k in 1:n_causes) {
      log_crcr_cond_prob_surv[k, patient_pos:patient_end] = rep_matrix(log_crcr_trial_lambda[k, s], n_trial_patients[s]);
        
      if (!no_prop_hazard) {
        patient_log_crcr_hazard_ratio[k, patient_pos:patient_end] =
          tumor_sum_covar[patient_pos:patient_end] * crcr_covar_trial_coef[k, s, :n_tumor_covar] + 
          covar_design_matrix[patient_pos:patient_end] * crcr_covar_trial_coef[k, s, (n_tumor_covar + 1):];
        
        log_crcr_cond_prob_surv[k, patient_pos:patient_end] += rep_matrix(patient_log_crcr_hazard_ratio[k, patient_pos:patient_end], max_confresp_week);
      }
      
      log_crcr_cond_prob_surv[k, patient_pos:patient_end] = - exp(log_crcr_cond_prob_surv[k, patient_pos:patient_end]); 
    }
    
    patient_pos = patient_end + 1;
  }
}

array[n_causes] matrix<upper = 1e-6>[n_patients, max_confresp_week] log_cif = calc_log_cif(log_crcr_cond_prob_surv);
