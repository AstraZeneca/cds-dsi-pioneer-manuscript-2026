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
  int<lower = 0, upper = 1> crcr_ignore_interval_censoring; // Treat observed intervals as true intervals 
  int<lower = 0, upper = 1> gen_log_lik;
  int<lower = 0, upper = 1> prior_sense;
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level_baseline_hazard;
  int<lower = 0, upper = 1> add_trial_level_prop_hazard;
  int<lower = 0, upper = 1 - add_trial_level_baseline_hazard> separate_baseline_hazard;
  int<lower = 0, upper = 1 - add_trial_level_prop_hazard> separate_prop_hazard;

  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  #include "crcr_data.stan"
  #include "../bootstrap/leave_out_trial_bootstrap_data.stan"
 
  // Hyperparam
  #include "crcr_hyperparam.stan"
}

transformed data {
  #include "../base_transformed_data.stan" 
  #include "crcr_transformed_data.stan"
  #include "../bootstrap/leave_out_trial_bootstrap_transformed_data.stan"
  
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
      
      if (leave_out_trial > 0) {
        target += reduce_sum(
          partial_sum_crcr_lupmf, last_unclassified_response_week[training_patients], grain_size,
          confirmed_response_cause[training_patients], 
          early_confirmed_response_censored[training_patients], 
          crcr_ignore_interval_censoring ? zeros_int_array(n_training_patients) : confirmed_response_interval_censored[training_patients], 
          log_crcr_cond_prob_surv[training_crcr_intervals], max_confresp_week
        );
      } else {
        target += reduce_sum(
          partial_sum_crcr_lupmf, last_unclassified_response_week, grain_size,
          confirmed_response_cause, 
          early_confirmed_response_censored, crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
          log_crcr_cond_prob_surv, max_confresp_week
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
}
