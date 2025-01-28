#include "../baseline_hazard/baseline_hazard_transformed_parameters.stan"
#include "../crcr/crcr_transformed_parameters.stan"

// Log conditional probability of survival at each interval. Separate columns for confirmed responders and non-responders, so we can calculate the mixture log likelihood. 
array[no_prop_hazard || pfs_only ? 1 : n_causes] matrix<upper = 0>[n_patients, max_all_t] log_cond_prob_surv; 
array[n_causes] vector[no_prop_hazard ? 0 : n_patients] time_invariant_log_hazard_ratio; // Log proportional hazard 

array[add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_tumor_covar + n_covar + 1] covar_trial_coef_residual; // Multilevel variations
array[no_prop_hazard ? 0 : n_trials] vector[n_tumor_covar + n_covar + 1] covar_trial_coef;

cholesky_factor_cov[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar + 1 : 1] L_covar_trial_cov;

if (add_trial_level_prop_hazard && !no_prop_hazard) {
  L_covar_trial_cov = diag_pre_multiply(covar_trial_sd, L_covar_trial_corr);
 
  for (s in 1:n_trials) { 
    covar_trial_coef_residual[s] = L_covar_trial_cov * raw_covar_trial_coef[s];
  }
} else {
  L_covar_trial_cov[1, 1] = 1; 
}

{ // Calculate patient-interval conditional probability of disease progression.
  int patient_pos = 1;
  int n_used_causes = no_prop_hazard ? 1 : n_causes;
  
  profile("log surv loop") {
    for (s in 1:n_trials) {
      int patient_end = patient_pos + n_trial_patients[s] - 1; 
      
      if (!no_prop_hazard) { 
        // Tumor size effect
        covar_trial_coef[s, :n_tumor_covar] = no_tumor_effects ? rep_vector(0, n_tumor_covar) : tumor_stim_pop_coef[separate_prop_hazard ? s : 1];
        // Othe patient level covariates
        covar_trial_coef[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)] = covar_effect[separate_prop_hazard ? s : 1];
        // Confirmed response status, to impute if not observed
        covar_trial_coef[s, n_tumor_covar + n_covar + 1] = pfs_only ? 0 : conf_resp_effect[separate_prop_hazard ? s : 1]; 
        
        if (add_trial_level_prop_hazard) {
          covar_trial_coef[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)] += covar_trial_coef_residual[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)];
          
          if (!no_tumor_effects) {
            covar_trial_coef[s, :n_tumor_covar] += covar_trial_coef_residual[s, :n_tumor_covar];
          }
          
          if (!pfs_only) {
            covar_trial_coef[s, n_tumor_covar + n_covar + 1] += covar_trial_coef_residual[s, n_tumor_covar + n_covar + 1];
          }
        }
       
        for (k in 1:n_causes) { 
          time_invariant_log_hazard_ratio[k, patient_pos:patient_end] = tumor_sum_covar[patient_pos:patient_end] * covar_trial_coef[s, :n_tumor_covar] + 
            covar_design_matrix[patient_pos:patient_end] * covar_trial_coef[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)];
        }
        
        if (!pfs_only) { 
          time_invariant_log_hazard_ratio[n_causes, patient_pos:patient_end] += covar_trial_coef[s, n_tumor_covar + n_covar + 1];
        }
      } 
      
      for (k in 1:n_used_causes) {
        log_cond_prob_surv[k, patient_pos:patient_end] = rep_matrix(log_trial_lambda[s], n_trial_patients[s]);
          
        if (!no_prop_hazard && !pfs_only) {
          log_cond_prob_surv[k, patient_pos:patient_end] += rep_matrix(time_invariant_log_hazard_ratio[k, patient_pos:patient_end], max_all_t);
        }
        
        log_cond_prob_surv[k, patient_pos:patient_end] = - exp(log_cond_prob_surv[k, patient_pos:patient_end]); 
      }
      
      patient_pos = patient_end + 1; 
    }
  }
}

