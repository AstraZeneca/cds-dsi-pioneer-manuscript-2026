functions {
  #include "../extern_util.stan"
  #include "../util.stan"
  #include "../extern_pfs_functions.stan"
  #include "../pfs_functions.stan"
  #include "crcr_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;

  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  #include "crcr_data.stan"

  // int<lower = 1> prediction_week; // At what week are starting our analysis
  
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
  #include "crcr_gen_quants.stan"
}
