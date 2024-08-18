functions {
  #include "util.stan"
  #include "pfs_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> gen_pfs; // Generate simulated data
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  int<lower = 0, upper = 1> time_varying_conf_resp;
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;
  int<lower = 0, upper = 1> add_tumor_location_level;

  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;
 
  int<lower = 0> n_covar; 
  matrix[n_patients, n_covar] covar_design_matrix;
  
  array[n_patients] int<lower = 1> experiment_start_week; // Week 1 is the first week of the experiment
  int<lower = 1> prediction_week; // At what week are starting our prediction 

  // I use these to generate 2d hazard ratios over tumor sizes observed over time. 
  int<lower = 0> n_grid_tumors;
  array[n_grid_tumors] int<lower = 1, upper = sum(n_patient_tumors)> grid_tumors;
  
  // Hyperparam
  #include "baseline_hazard_hyperparam.stan"
  
  vector<lower = 0>[2] tumor_stim_pop_coef_sd;
  real<lower = 0> conf_resp_effect_sd;
  vector<lower = 0>[n_covar] covar_effect_sd;
}

transformed data {
  #include "tumor/tumor_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
}

parameters {
  #include "baseline_hazard_parameters.stan"
}


transformed parameters {
  #include "baseline_hazard_transformed_parameters.stan"
  
  vector<upper = 0>[n_time_periods] log_cond_prob_surv;
  
  profile("log_cond_prob") { // Calculate patient-interval conditional probability of disease progression.
    int pfs_interval_pos = 1;
    int conresp_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int n_intervals = gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i];
      int pfs_interval_end = pfs_interval_pos + n_intervals - 1;
      
      vector[n_intervals] hazard_ratio_pred = log_trial_lambda[patient_trial[i], 1:n_intervals];  

      log_cond_prob_surv[pfs_interval_pos:pfs_interval_end] = - exp(hazard_ratio_pred);

      pfs_interval_pos = pfs_interval_end + 1;
    }
  }
}

model {
  // Priors
  
  #include "baseline_hazard_priors.stan"
  
  // Likelihood
  
  profile("likelihood") {
    if (fit_data) {
      pfs ~ pch2(right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv, gen_pfs ? max_all_t : 0, rep_array(1, n_patients));
    }
  }
}

generated quantities {
  array[gen_pfs ? n_patients : 0] int<lower = 0> sim_pfs; 
  array[gen_pfs ? n_patients : 0] int<lower = 0, upper = 1> sim_censored; 
  
  real<lower = 0> sim_median_pfs = 0; 
  vector<lower = 0, upper = 1>[gen_pfs ? max_all_t + 1 : 0] km_est; 
  
  profile("gen_pfs") {
    if (gen_pfs) {
      int pfs_interval_pos = 1;
      
      for (i in 1:n_patients) {
        int n_intervals = max_all_t;
        int pfs_interval_end = pfs_interval_pos + n_intervals - 1;
       
        (sim_pfs[i], sim_censored[i]) = survival_time_rng(log_cond_prob_surv[pfs_interval_pos:pfs_interval_end]); 
        
        pfs_interval_pos = pfs_interval_end + 1;
      }
      
      sim_median_pfs = survival_median(sim_pfs, max_all_t).1; 
      km_est = estimate_kaplan_meier(sim_pfs, sim_censored, max_all_t).1; 
    }
  }
}

