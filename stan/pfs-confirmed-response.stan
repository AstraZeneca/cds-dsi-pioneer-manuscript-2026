functions {
  #include "extern_util.stan"
  #include "util.stan"
  #include "pfs_functions.stan"
  #include "extern_pfs_functions.stan"
  #include "crcr/crcr_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  // int<lower = 0, upper = 1> time_varying_conf_resp;
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;

  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;
 
  #include "crcr/crcr_data.stan" 
  
  int<lower = 0> n_bootstrap_param;
  array[n_bootstrap_param] int<lower = 1> prediction_week; // At what week are starting our prediction 
  // array[n_bootstrap_param] real<lower = 0> recruit_lambda; // neg binom rate
  // real<lower = 0> recruit_phi; // neg binom dispersion 

  // Hyperparam
  #include "baseline_hazard_hyperparam.stan"
  #include "crcr/crcr_hyperparam.stan"
  
  vector<lower = 0>[2] tumor_stim_pop_coef_sd;
  real<lower = 0> conf_resp_effect_sd;
  vector<lower = 0>[n_covar] covar_effect_sd;
 
  real<lower = 0> covar_trial_sd_sd;
  real<lower = 0> covar_trial_corr_eta; 
}

transformed data {
  int gen_pfs = 1;
  
  #include "base_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
  #include "crcr/crcr_transformed_data.stan"
  
  int n_missing_confirmed_response = sum(confirmed_response_censored); 
  int n_obs_confirmed_response = n_patients - n_missing_confirmed_response; 
  array[n_obs_confirmed_response] int<lower = 1, upper = n_patients> obs_confirmed_response;
  array[n_missing_confirmed_response] int<lower = 1, upper = n_patients> missing_confirmed_response;
  
  (obs_confirmed_response, missing_confirmed_response) = get_mask_idx(confirmed_response_censored);
  
  array[n_patients] int<lower = 1> sorted_experiment_start_week = sort_asc(experiment_start_week);
  
  int n_marg_prob_t = 5;
  array[n_marg_prob_t] int marg_prob_t = { 1, 5, 10, 50, 100 };
  int n_marg_prob_patients = 5;
  
  int crcr_grain_size = 83;
}

parameters {
  #include "baseline_hazard_parameters.stan"
  #include "crcr/crcr_parameters.stan"
  #include "recruit/recruit_parameters.stan"
  
  vector[n_tumor_covar] tumor_stim_pop_coef;
  vector[n_covar] covar_effect;  
  real conf_resp_effect;
  
  vector<lower = 0>[add_trial_level ? n_tumor_covar + n_covar + 1 : 0] covar_trial_sd;
  // cholesky_factor_corr[add_trial_level ? n_tumor_covar + n_covar + 1 : 0] L_covar_trial_corr;
  matrix[n_tumor_covar + n_covar + 1, add_trial_level ? n_trials : 0] raw_covar_trial_coef;
}


transformed parameters {
  #include "baseline_hazard_transformed_parameters.stan"
  #include "crcr/crcr_transformed_parameters.stan"
  
  matrix<lower = 0, upper = 1>[n_patients, n_causes] conf_resp_prob; 
  
  matrix<upper = 0>[n_time_periods, n_causes] log_cond_prob_surv;
  matrix[n_patients, n_causes] time_invariant_log_hazard_ratio = rep_matrix(tumor_sum_covar * tumor_stim_pop_coef + covar_design_matrix * covar_effect, n_causes);
  
  time_invariant_log_hazard_ratio[, n_causes] += conf_resp_effect; 
  
  matrix[n_tumor_covar + n_covar + 1, add_trial_level ? n_trials : 0] covar_trial_coef;
  
  if (add_trial_level) {
    covar_trial_coef = diag_pre_multiply(covar_trial_sd, raw_covar_trial_coef); 
    // covar_trial_coef = diag_pre_multiply(covar_trial_sd, L_covar_trial_corr) * raw_covar_trial_coef; 
  } 
  
  profile("log_cond_prob_surv") { // Calculate patient-interval conditional probability of disease progression.
    array[n_patients] matrix[max_confresp_week, n_causes] cif;
    
    (cif, conf_resp_prob) = calc_cif(n_patients, log_crcr_cond_prob_surv, max_confresp_week);
    
    int pfs_interval_pos = 1;
    int patient_pos = 1;
  
    for (s in 1:n_trials) {
      int patient_end = patient_pos + n_trial_patients[s] - 1; 
     
      if (add_trial_level) { 
        time_invariant_log_hazard_ratio[patient_pos:patient_end] += rep_matrix(
          tumor_sum_covar[patient_pos:patient_end] * covar_trial_coef[:n_tumor_covar, s] + 
          covar_design_matrix[patient_pos:patient_end] * covar_trial_coef[(n_tumor_covar + 1):(n_tumor_covar + n_covar), s], 
          n_causes 
        );
          
        time_invariant_log_hazard_ratio[patient_pos:patient_end, n_causes] += covar_trial_coef[n_tumor_covar + n_covar + 1, s]; 
      }
  
      for (i in patient_pos:patient_end) {
        int n_intervals = max_all_t;
        int pfs_interval_end = pfs_interval_pos + n_intervals - 1;
        
        log_cond_prob_surv[pfs_interval_pos:pfs_interval_end] = 
          - exp(rep_matrix(log_trial_lambda[s, 1:n_intervals], n_causes) + rep_matrix(time_invariant_log_hazard_ratio[i], n_intervals));
  
        pfs_interval_pos = pfs_interval_end + 1;
      }
      
      patient_pos = patient_end + 1; 
    }
  }
}

model {
  // Priors
  
  profile("crcr priors") {
    #include "crcr/crcr_priors.stan"
  }
  
  #include "baseline_hazard_priors.stan"
  #include "recruit/recruit_priors.stan"
  
  profile("pfs priors") {
    tumor_stim_pop_coef ~ normal(0, tumor_stim_pop_coef_sd);
    conf_resp_effect ~ normal(0, conf_resp_effect_sd); 
    covar_effect ~ normal(0, covar_effect_sd);
    
    if (add_trial_level) {
      covar_trial_sd ~ normal(0, covar_trial_sd_sd);
      // L_covar_trial_corr ~ lkj_corr_cholesky(covar_trial_corr_eta);
      to_vector(raw_covar_trial_coef) ~ std_normal(); 
    }
  }
  
  // Likelihood
  
  if (fit_data) {
    profile("crcr loglik") {
      target += reduce_sum(
        partial_sum_crcr_lupmf, last_unclassified_response_week, crcr_grain_size,
        confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week
      );
    }
    
    // last_unclassified_response_week ~ comp_risk_pch(confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week);
    
    experiment_start_week ~ neg_binomial_2(recruit_lambda[patient_trial], recruit_phi[patient_trial]);  
  
    profile("pfs loglik") {  
      matrix[n_patients, n_causes] response_lp = append_col( 
        calc_pch_loglik2(
          pfs, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv[, 1], max_all_t, rep_array(1, n_patients)
        ),
        calc_pch_loglik2(
          pfs, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv[, 2], max_all_t, rep_array(1, n_patients)
        )
      );
      
      for (i in 1:n_patients) {
        if (confirmed_response_censored[i]) { // Unclassified
          target += log_mix(conf_resp_prob[i, 1], response_lp[i, 1], response_lp[i, 2]);
        } else {
          target += response_lp[i, confirmed_response_cause[i]];
        }
      }
    }
  }
}

generated quantities {
  #include "crcr/crcr_gen_quants.stan"
  
  array[n_marg_prob_patients, n_causes] vector<lower = 0, upper = 1>[n_marg_prob_t] marginal_exit_prob;
  
  {
    int pfs_interval_pos = 1;
    
    for (i in 1:n_marg_prob_patients) {
      int pfs_interval_end = pfs_interval_pos + max(marg_prob_t) - 1;
      
      for (k in 1:n_causes) {
        marginal_exit_prob[i, k] = calculate_marginal_exit_prob(log_cond_prob_surv[pfs_interval_pos:pfs_interval_end, k], max(marg_prob_t))[marg_prob_t];
      }
      
      pfs_interval_pos += max_all_t;
    }
  }
  
  array[n_patients] int<lower = 0, upper = 1> sim_confirmed_response = confirmed_response;
  
  sim_confirmed_response[missing_confirmed_response] = bernoulli_rng(conf_resp_prob[missing_confirmed_response, 2]);
  
  array[n_patients] int<lower = 0> sim_pfs; 
  array[n_patients] int<lower = 0, upper = 1> sim_censored; 
  
  real<lower = 0> sim_median_pfs;
  
  vector<lower = 0, upper = 1>[max_all_t + 1] km_est; 
  
  vector<lower = 0>[n_bootstrap_param] bootstrap_median_pfs;
  array[n_bootstrap_param] int<lower = 0, upper = n_patients> n_bootstrap_sample;
  vector<lower = 0, upper = 1>[n_bootstrap_param] bootstrap_maturity_rate;
 
  { 
    int pfs_interval_pos = 1;
    array[n_causes] matrix[max_all_t, n_patients] mat_log_cond_prob_surv;
    
    profile("gen_pfs") {
      for (i in 1:n_patients) {
        int n_intervals = max_all_t;
        int pfs_interval_end = pfs_interval_pos + n_intervals - 1;
        
        for (k in 1:n_causes) {
          mat_log_cond_prob_surv[k, , i] = log_cond_prob_surv[pfs_interval_pos:pfs_interval_end, k];
        }
        
        (sim_pfs[i], sim_censored[i]) = survival_time_rng(mat_log_cond_prob_surv[sim_confirmed_response[i] + 1, , i]); 
        
        pfs_interval_pos = pfs_interval_end + 1;
      }
      
      sim_median_pfs = survival_median(sim_pfs, max_all_t).1; 
      
      km_est = estimate_kaplan_meier(sim_pfs, sim_censored, max_all_t).1; 
    }
    
    profile("bootstrap") {
      for (p in 1:n_bootstrap_param) {
        bootstrap_maturity_rate[p] = 0; 
        
        array[n_patients] int experiment_start = sort_asc(neg_binomial_2_rng(rep_vector(recruit_lambda[1], n_patients), rep_vector(recruit_phi[1], n_patients)));
        array[n_patients] int current_bootstrap_pfs;
        
        int i = 0;
        
        while ((i + 1 <= n_patients) && (experiment_start[i + 1] <= prediction_week[p])) {
          i += 1;
          
          int current_patient = discrete_range_rng(1, n_patients);
          int bootstrap_confirmed_response; 
          
          if (experiment_start[i] + confirmed_response_week[current_patient] - 1 <= prediction_week[p]) {
            bootstrap_confirmed_response = confirmed_response[current_patient];
          } else {
            // Not observed yet, so let's estimate it.
            bootstrap_confirmed_response = bernoulli_rng(conf_resp_prob[current_patient, 2]);
          }
          
          bootstrap_maturity_rate[p] += bootstrap_confirmed_response;
          current_bootstrap_pfs[i] = survival_time_rng(mat_log_cond_prob_surv[bootstrap_confirmed_response + 1, , current_patient]).1; 
        }
        
        n_bootstrap_sample[p] = i;
        
        if (n_bootstrap_sample[p] > 0) { 
          bootstrap_median_pfs[p] = survival_median(current_bootstrap_pfs[:n_bootstrap_sample[p]], max_all_t).1;  
          bootstrap_maturity_rate[p] /= n_bootstrap_sample[p];
        } else {
          bootstrap_median_pfs[p] = 0; 
          bootstrap_maturity_rate[p] = 0;
        }
      } 
    }
  }
}

