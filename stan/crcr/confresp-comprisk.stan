/*
 * Standalone script for the confirmed response competing risk survival (CRCR) model. This is mostly helpful when debugging the CRCR model in isolation. Most of the code
 * is included from other files.
 */

functions {
  #include "../util.stan"
  #include "../pfs_functions.stan"
  #include "crcr_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> crcr_ignore_interval_censoring; // Treat observed intervals as true intervals 
  int<lower = 0, upper = 1> gen_log_lik; // Generate log likelihood
  int<lower = 0, upper = 1> prior_sense; // Calculate prior sensitivity info for {priorsense}
  int<lower = 0, upper = 1> no_prop_hazard; 
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level_baseline_hazard;
  int<lower = 0, upper = 1> add_trial_level_prop_hazard;
  int<lower = 0, upper = 1> separate_baseline_hazard;
  int<lower = 0, upper = 1> separate_prop_hazard;

  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  #include "crcr_data.stan"
 
  // Calculating log likelihood for a single trial. Useful if you want to compare the preformance of a model using a single trial with one that is multilevel. 
  int<lower = 0, upper = n_trials> log_lik_trial; 
 
  // Hyperparam
  #include "crcr_hyperparam.stan"
}

transformed data {
  int n_causes = 2; // We only have confirmed response and non-response
  
  #include "../base_transformed_data.stan" 
  #include "crcr_transformed_data.stan"
  
  int grain_size = 83; // For reduce_sum()
  int leave_out_trial = 0;
  int n_training_patients = n_patients;
}

parameters {
  #include "crcr_parameters.stan"
}

transformed parameters {
  #include "crcr_transformed_parameters.stan"
}

model {
  profile("priors") {
    #include "crcr_priors.stan" 
  }
  
  if (fit_data) {
    profile("loglik") {
      // last_unclassified_response_week ~ comp_risk_pch(confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week);
      
      if (leave_out_trial > 0) {
        fatal_error("This code is out of date and needs to be fixed. The structure of log_crcr_cond_prob_surv has changed.");
        // target += reduce_sum(
        //   partial_sum_crcr_lupmf, last_unclassified_response_week[training_patients], grain_size,
        //   confirmed_response_cause[training_patients], 
        //   confirmed_response_censored[training_patients], 
        //   crcr_ignore_interval_censoring ? zeros_int_array(n_training_patients) : confirmed_response_interval_censored[training_patients], 
        //   log_crcr_cond_prob_surv[, training_crcr_intervals]
        // );
      } else {
        target += reduce_sum(
          partial_sum_crcr_lupmf, last_unclassified_response_week, grain_size,
          confirmed_response_cause, 
          confirmed_response_censored, 
          crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
          log_crcr_cond_prob_surv
        );
      }
    }
  }
}

generated quantities {
  #include "crcr_gen_quants.stan"
  
  vector[gen_log_lik || prior_sense ? n_training_patients : 0] log_lik = rep_vector(0, gen_log_lik || prior_sense ? n_training_patients : 0);
  real lprior = 0;

  #include "crcr_log_lik_prior_sense.stan"
  
  vector[(gen_log_lik || prior_sense) && log_lik_trial > 0 && leave_out_trial == 0 ? n_trial_patients[log_lik_trial] : 0] trial_log_lik;
  
  if (rows(trial_log_lik) > 0) {
    int pos = trial_patient_pos[log_lik_trial];
    int end = trial_patient_pos[log_lik_trial + 1] - 1;
    
    trial_log_lik = log_lik[pos:end];
  }
}
