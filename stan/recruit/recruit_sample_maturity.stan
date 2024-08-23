functions {
  #include "../extern_util.stan"
  #include "../util.stan"
  #include "../extern_pfs_functions.stan"
  #include "../pfs_functions.stan"
  #include "../crcr/crcr_functions.stan"
}

data {
  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  #include "../crcr/crcr_data.stan"
  
  int<lower = 0> n_covar; 
  matrix[n_patients, n_covar] covar_design_matrix;
  
  // Recruitment simulation data
  
  int<lower = 1> n_lambda;
  vector<lower = 0>[n_lambda] lambda;
  
  int<lower = 1> n_pred_week;
  array[n_pred_week] int<lower = 1> pred_week;
  
  real<lower = 0> phi;
 
  // Hyperparam
  #include "../crcr/crcr_hyperparam.stan"
}

transformed data {
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level = 0;
  int<lower = 0, upper = 1> add_tumor_location_level = 0;
  
  #include "../tumor/tumor_transformed_data.stan" 
  #include "../crcr/crcr_transformed_data.stan"
}

parameters {
  #include "../crcr/crcr_parameters.stan"
}

transformed parameters {
  #include "../crcr/crcr_transformed_parameters.stan"
}

model {
  // Priors
 
  #include "../crcr/crcr_priors.stan" 
  
  // Likelihood
  
  last_unclassified_response_week ~ comp_risk_pch(confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week);
}

generated quantities {
  array[n_lambda, n_pred_week] int<lower = 0> n_sample = rep_array(0, n_lambda, n_pred_week);
  array[n_lambda, n_pred_week] real<lower = 0, upper = 1> maturity_rate;
  
  for (l in 1:n_lambda) {
    for (p in 1:n_pred_week) {
      array[n_patients] int experiment_start = neg_binomial_2_rng(rep_vector(lambda[l], n_patients), rep_vector(phi, n_patients));
     
      int n_matured_patients = 0;
      
      for (i in 1:n_patients) {
        int n_intervals = max_confresp_week; 
        
        int censored; 
        int maturity_week;
        int cause;
        
        (maturity_week, censored, cause) = competing_risks_survival_time_rng(-exp(log_crcr_trial_lambda[patient_trial[i], 1:n_intervals])); 
        
        if (experiment_start[i] <= pred_week[p]) {
          n_sample[l, p] += 1;
        } 
        
        if (experiment_start[i] + maturity_week <= pred_week[p] && !censored) {
          n_matured_patients += 1;
        }
      }
     
      if (n_sample[l, p] > 0) { 
        maturity_rate[l, p] = n_matured_patients * 1.0 / n_sample[l, p];
      } else {
        maturity_rate[l, p] = 0; 
      }
    }
  }
}
