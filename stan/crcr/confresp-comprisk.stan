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
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed intervals as true intervals 
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;

  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  #include "crcr_data.stan"
 
  // Hyperparam
  #include "crcr_hyperparam.stan"
}

transformed data {
  #include "../base_transformed_data.stan" 
  #include "crcr_transformed_data.stan"
  
  int grain_size = 83;
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
      target += reduce_sum(
        partial_sum_crcr_lupmf, last_unclassified_response_week, grain_size,
        confirmed_response_cause, 
        early_confirmed_response_censored, ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
        log_crcr_cond_prob_surv, max_confresp_week
      );
    }
  }
}

generated quantities {
  #include "crcr_gen_quants.stan"
}
