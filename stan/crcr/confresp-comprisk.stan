functions {
  #include "../util.stan"
  #include "../pfs_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;
  int<lower = 0, upper = 1> add_tumor_location_level;

  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  
  array[n_patients] int<lower = 0, upper = 1> confirmed_response;
  array[n_patients] int<lower = 1> confirmed_response_week;
  array[n_patients] int<lower = 0, upper = 1> confirmed_response_censored;
  
  array[n_patients] int<lower = 1> experiment_start_week;
  int<lower = 1> prediction_week; // At what week are starting our analysis
  
  int<lower = 0> n_covar; 
  matrix[n_patients, n_covar] covar_design_matrix;
 
  // Hyperparam
  #include "crcr_hyperparam.stan"
}

transformed data {
  #include "../tumor/tumor_transformed_data.stan" 
  #include "crcr_transformed_data.stan"
}

parameters {
  #include "crcr_parameters.stan"
}

transformed parameters {
  #include "crcr_transformed_parameters.stan"
}

model {
  // Priors
 
  #include "crcr_priors.stan" 
  
  // Likelihood
  
  if (fit_data) {
    last_unclassified_response_week ~ comp_risk_pch(confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week);
  }
}

generated quantities {
  array[n_patients] matrix[max_confresp_week, n_causes] cif; // cumulative incidence function
  array[n_patients] simplex[n_causes] prob_cause; 
  
  for (i in 1:n_patients) { 
    int patient_prob_pos = 1 + (i - 1) * max_confresp_week; 
    
    for (t in 1:max_confresp_week) {
      cif[i, t] = 
        exp(sum(log_crcr_cond_prob_surv[patient_prob_pos:(patient_prob_pos + t - 2)]) + log1m_exp(log_crcr_cond_prob_surv[patient_prob_pos + t - 1])); 
    }
  
    for (k in 1:n_causes) {
      cif[i, , k] = cumulative_sum(cif[i, , k]);
    }
    
    prob_cause[i] = cif[i, max_confresp_week]';
    prob_cause[i] /= sum(prob_cause[i]); 
  }
}
