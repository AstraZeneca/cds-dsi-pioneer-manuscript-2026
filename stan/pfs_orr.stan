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
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;
  
  #include "base_data.stan"
  #include "tumor_data.stan"
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0, upper = 1> orr;
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  
  // Hyperparam
  
  #include "baseline_hazard_hyperparam.stan"

  real<lower = 0> orr_pop_coef_sd;
}

transformed data {
  #include "tumor_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
  
  array[n_patients] int patient_max_2nd_tumor_t = rep_array(1, n_patients); // Using this for the loglik calculation in "generated quantities".
}

parameters {
  #include "baseline_hazard_parameters.stan"
  
  real orr_coef;
}

transformed parameters {
  #include "baseline_hazard_transformed_parameters.stan"
  
  vector[n_time_periods] disease_progress_pred; // DP predictor for all patients at all observed and unobserved intervals.
  
  { // Calculate patient-interval conditional probability of disease progression.
    int pfs_interval_pos = 1;

    for (i in 1:n_patients) {
      int n_intervals = gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i];
      int pfs_interval_end = pfs_interval_pos + n_intervals - 1;

      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = log_trial_lambda[patient_trial[i], 1:n_intervals] + orr_coef * orr[i];

      pfs_interval_pos = pfs_interval_end + 1;
    }
  }
  
  vector<lower = 0, upper = 1>[n_time_periods] disease_progress_prob = inv_cloglog(disease_progress_pred); // DP conditional probability for all patients (same as above)
}

model {
  // Priors
  
  #include "baseline_hazard_priors.stan"
  
  orr_coef ~ normal(0, orr_pop_coef_sd);
  
  // Likelihood
  
  if (fit_data) {
    pfs ~ pch(right_uncensored, interval_censored, ignore_interval_censoring, disease_progress_prob, gen_pfs ? max_all_t : 0, rep_array(1, n_patients)); 
  }
}

generated quantities {
  #include "pfs_generated_quant.stan"
   
  real<lower = 0, upper = 1> base_pfs_6mon_prob;
  real<lower = 0, upper = 1> base_pfs_9mon_prob;
  real<lower = 0, upper = 1> orr_pfs_6mon_prob;
  real<lower = 0, upper = 1> orr_pfs_9mon_prob;
  
  {
    vector[9 * 4] base_pf_cond_lp = -exp(log_lambda[:(9 * 4)]); 
    vector[9 * 4] orr_pf_cond_lp = -exp(log_lambda[:(9 * 4)] + orr_coef); 
   
    vector[9 * 4] base_pf_marginal_prob; 
    vector[9 * 4] orr_pf_marginal_prob; 
    
    real current_base_pf_lp = 0; // Marginal
    real current_orr_pf_lp = 0; // Marginal
    
    for (t in 1:(9 * 4)) {
      base_pf_marginal_prob[t] = exp(current_base_pf_lp + log1m_exp(base_pf_cond_lp[t]));
      orr_pf_marginal_prob[t] = exp(current_orr_pf_lp + log1m_exp(orr_pf_cond_lp[t]));
      
      current_base_pf_lp += base_pf_cond_lp[t];
      current_orr_pf_lp += orr_pf_cond_lp[t];
    } 
    
    base_pfs_6mon_prob = 1 - sum(base_pf_marginal_prob[:(6 * 4)]);
    base_pfs_9mon_prob = 1 - sum(base_pf_marginal_prob);
    orr_pfs_6mon_prob = 1 - sum(orr_pf_marginal_prob[:(6 * 4)]);
    orr_pfs_9mon_prob = 1 - sum(orr_pf_marginal_prob);
  }
  
  real<lower = -1, upper = 1> pfs_6mon_prob_diff = orr_pfs_6mon_prob - base_pfs_6mon_prob;
  real<lower = -1, upper = 1> pfs_9mon_prob_diff = orr_pfs_9mon_prob - base_pfs_9mon_prob;
}
