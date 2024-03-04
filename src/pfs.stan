functions {
  #include "util.stan"
  #include "pfs_functions.stan"
}

data {
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_pfs;
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  int<lower = 0, upper = 1> use_tumor_model;
  int<lower = 0, upper = 1> add_trial_level;
  
  // 0: None; 1: First two, no interaction; 2: First two, interaction; 3: First two percentage difference; 4: tumor intercept only
  // 5: quadratic
  int<lower = 0, upper = 5> tumor_hazard_type; 
  
  #include "base_data.stan"
  #include "tumor_data.stan"
    
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  
  // Hyperparam
  
  real log_lambda_gp_intercept_mean;
  real<lower = 0> log_lambda_gp_intercept_sd;
  real<lower = 0> log_lambda_gp_alpha_sd;
  real<lower = 0> log_lambda_gp_trial_alpha_sd;
  real<lower = 0> log_lambda_gp_rho_alpha;
  real<lower = 0> log_lambda_gp_rho_beta;
  real<lower = 0> log_lambda_gp_trial_intercept_sd_sd;
  real<lower = 0> tumor_stim_intercept_sd;
  vector<lower = 0>[5] tumor_stim_coef_sd;
  real<lower = 0> tumor_stim_trial_coef_sd_sd;
}

transformed data {
  #include "tumor_transformed_data.stan" 
  
  if (min(n_patient_screening_t) <= 0) {
    reject("We need at least one screening measurement per patient, for now: ", min(n_patient_screening_t));
  }
  
  
  array[n_trials] int<lower = 0, upper = n_patients> n_trial_patients = rep_array(0, n_trials);
  array[n_patients] int<lower = 1, upper = n_patients> ordered_trial_patients;
 
  {
    array[n_trials] int trial_patient_pos;
    trial_patient_pos[1] = 1;
 
    for (i in 1:n_patients) {
      n_trial_patients[patient_trial[i]] += 1;
    }
    
    
    for (s in 2:n_trials) {
      trial_patient_pos[s] = sum(n_trial_patients[:(s - 1)]) + 1;
    }
    
    for (i in 1:n_patients) {
      ordered_trial_patients[trial_patient_pos[patient_trial[i]]] = i;
      trial_patient_pos[patient_trial[i]] += 1;
    }
  } 
  
  int<lower = 0, upper = max(pfs) * n_patients> n_total_pfs = sum(pfs);
  
  // These are used for the time interval distance between base hazard
  array[max_all_t] real pfs_range;
  array[max_all_t] int pfs_range_int;
  vector[max_all_t] pfs_range_vec;
  
  for (i in 1:max_all_t) {
    pfs_range[i] = i / 12.0;
    pfs_range_int[i] = i;
    pfs_range_vec[i] = i - 1; 
  }
  
  int n_covar_col = 0;
  
  if (tumor_hazard_type == 1) {
    n_covar_col = 2;
  } else if (tumor_hazard_type == 2) {
    n_covar_col = 3;
  } else if (tumor_hazard_type == 3) {
    n_covar_col = 1;
  } else if (tumor_hazard_type == 5) {
    n_covar_col = 5;
  }
  
  matrix[sum(n_patient_tumors), n_covar_col] tumor_covar;
  
  if (tumor_hazard_type > 0) {
    if (tumor_hazard_type == 1 || tumor_hazard_type == 2 || tumor_hazard_type == 5) {  
      tuple(real, real, vector[sum(n_measures)]) standardize_results = standardize_nonzero_tumor_sizes(tumor_size);
      vector[sum(n_measures)] standardized_tumor_size = standardize_results.3;
      tumor_covar[, 1:2] = prepare_early_tumors_design_matrix(standardized_tumor_size, n_patient_tumors, n_measures, n_screening_t);
      
      if (tumor_hazard_type == 2 || tumor_hazard_type == 5) {
        tumor_covar[, 3] = tumor_covar[, 1] .* tumor_covar[, 2];
      }
      
      if (tumor_hazard_type == 5) {
        tumor_covar[, 4] = tumor_covar[, 1]^2;
        tumor_covar[, 5] = tumor_covar[, 2]^2;
      }
    } else if (tumor_hazard_type == 3) {
      matrix[sum(n_patient_tumors), 2] covar = prepare_early_tumors_design_matrix(tumor_size, n_patient_tumors, n_measures, n_screening_t);
      
      tumor_covar[, 1] = (covar[, 2] - covar[, 1]) ./ covar[, 1];
    } 
  }
  
  // Censoring information calculated from PFS and t_measure; no need to pass it in. 
  array[n_patients] int<lower = 0> interval_censored = rep_array(0, n_patients);
  array[n_patients] int<lower = 0, upper = 1> right_uncensored = rep_array(0, n_patients);

  {
    tuple(array[n_patients] int, array[n_patients] int) censoring_res = identify_censoring(pfs, death_week, n_patient_tumors, n_measures, t_measure);
    interval_censored = censoring_res.1;
    
    for (i in 1:n_patients) {
        right_uncensored[i] = 1 - right_censored[i];
        
        if (interval_censored[i] > 0 && right_censored[i]) {
          reject("Cannot be both interval and right censored.");
        }
    }
      
    print("Number of interval censored observations: ", sum(interval_censored));
    print("Number of right censored observations: ", sum(right_censored));
  }
 
  // If generating PFS we need to calculate probs for all possible time intervals, otherwise only up to observed PFS. 
  int<lower = 0> n_time_periods = gen_pfs ? n_patients * max_all_t : n_total_pfs + sum(right_uncensored) + sum(interval_censored);
}

parameters {
  #include "tumor_parameters.stan"
  
  // Base hazard GP parameters
  real<lower = 0> log_lambda_gp_alpha;
  real<lower = 0> log_lambda_gp_rho;
  vector[max_all_t] log_lambda_gp_eta;
  real log_lambda_gp_intercept;
  
  vector[add_trial_level ? n_trials : 0] raw_log_lambda_gp_trial_intercept;
  real<lower = 0> log_lambda_gp_trial_intercept_sd;
  
  real<lower = 0> log_lambda_gp_trial_alpha;
  real<lower = 0> log_lambda_gp_trial_rho;
  array[add_trial_level ? n_trials : 0] vector[max_all_t] log_lambda_gp_trial_eta;
 
  // Tumor level influence on hazard 
  real<lower = 0> tumor_stim_intercept; // DO NOT REMOVE; this is a per tumor intercept and not per patient intercept which is included in lambda.
  vector<lower = 0>[n_covar_col] tumor_stim_coef;
  
  array[add_trial_level ? n_trials : 0] vector[n_covar_col + 1] raw_tumor_stim_trial_coef_mult;
  vector<lower = 0>[add_trial_level ? n_covar_col + 1 : 0] tumor_stim_trial_coef_mult_sd;
}

transformed parameters {
  #include "tumor_transformed_parameters.stan"
  
  // Base hazard
  vector[max_all_t] log_lambda = calc_gp_pred(pfs_range, log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, delta, log_lambda_gp_eta);
  vector[add_trial_level ? n_trials : 0] log_lambda_gp_trial_intercept;
  array[n_trials] vector[max_all_t] log_trial_lambda; 
  
  vector<lower = 0, upper = 1>[n_time_periods] disease_progress_prob;
  
  vector[n_trials] tumor_stim_trial_intercept;
  array[n_trials] vector[n_covar_col] tumor_stim_trial_coef;
  
  if (add_trial_level) {
    log_lambda_gp_trial_intercept = raw_log_lambda_gp_trial_intercept * log_lambda_gp_trial_intercept_sd;
  }
  
  for (s in 1:n_trials) {
    tumor_stim_trial_intercept[s] = tumor_stim_intercept * exp(add_trial_level ? raw_tumor_stim_trial_coef_mult[s, 1] * tumor_stim_trial_coef_mult_sd[1] : 0);
   
    if (n_covar_col > 0) { 
      // TODO Add correlation between intercept and the coefs
      tumor_stim_trial_coef[s] = 
        tumor_stim_coef .* (add_trial_level ? exp(raw_tumor_stim_trial_coef_mult[s, 2:] .* tumor_stim_trial_coef_mult_sd[2:]) : rep_vector(1, n_covar_col));
    }
    
    log_trial_lambda[s] = log_lambda + ( 
      add_trial_level ? 
      calc_gp_pred(pfs_range, log_lambda_gp_trial_intercept[s], log_lambda_gp_trial_alpha, log_lambda_gp_trial_rho, delta, log_lambda_gp_trial_eta[s]): 
      rep_vector(0, max_all_t));
  }
  
  vector[n_patients] total_time_invar_tumor_stim;
  vector[n_time_periods] disease_progress_pred;
  
  { // Calculate patient-interval conditional probability of disease progression.
    
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
  
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int n_intervals = gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i];
      int pfs_interval_end = pfs_interval_pos + n_intervals - 1; 
      
      vector[n_patient_tumors[i]] tumor_stim = linear_tumor_stimulus(
        tumor_stim_trial_intercept[patient_trial[i]], 
        tumor_stim_trial_coef[patient_trial[i]], 
        tumor_covar[tumor_pos:tumor_end]
      );
      
      total_time_invar_tumor_stim[i] = tumor_hazard_type > 0 ? sum(tumor_stim) - mean(n_patient_tumors) * tumor_stim_intercept : 0;
     
      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = log_trial_lambda[patient_trial[i], 1:n_intervals] + total_time_invar_tumor_stim[i];
      
      tumor_pos = tumor_end + 1;
      pfs_interval_pos = pfs_interval_end + 1;
    }
      
    disease_progress_prob = inv_cloglog(disease_progress_pred); 
  }
}

model {
  #include "tumor_model.stan"
  
  // Priors
  
  log_lambda_gp_alpha ~ normal(0, log_lambda_gp_alpha_sd);
  log_lambda_gp_rho ~ inv_gamma(log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);
  log_lambda_gp_eta ~ std_normal();
  log_lambda_gp_intercept ~ normal(log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd);
  
  log_lambda_gp_trial_intercept_sd ~ normal(0, log_lambda_gp_trial_intercept_sd_sd);
  raw_log_lambda_gp_trial_intercept ~ std_normal(); 
  
  tumor_stim_intercept ~ normal(0, tumor_stim_intercept_sd);
  
  // TODO separate hyperparam for these parameters 
  log_lambda_gp_trial_alpha ~ normal(0, log_lambda_gp_trial_alpha_sd);
  log_lambda_gp_trial_rho ~ inv_gamma(log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);
 
  if (add_trial_level) { 
    tumor_stim_trial_coef_mult_sd ~ normal(0, tumor_stim_trial_coef_sd_sd);
    
    for (s in 1:n_trials) {
      raw_tumor_stim_trial_coef_mult[s] ~ std_normal();
      log_lambda_gp_trial_eta[s] ~ std_normal();
    }
  }
  
  if ((tumor_hazard_type > 0 && tumor_hazard_type < 4) || tumor_hazard_type > 4) {
    tumor_stim_coef[1] ~ normal(0, tumor_stim_coef_sd[1]);

    if (tumor_hazard_type == 1 || tumor_hazard_type == 2) {
      tumor_stim_coef[2] ~ normal(0, tumor_stim_coef_sd[2]);

      if (tumor_hazard_type == 2 || tumor_hazard_type == 5) {
        tumor_stim_coef[3] ~ normal(0, tumor_stim_coef_sd[3]);
      }
      
      if (tumor_hazard_type == 5) {
        tumor_stim_coef[4] ~ normal(0, tumor_stim_coef_sd[4]);
        tumor_stim_coef[5] ~ normal(0, tumor_stim_coef_sd[5]);
      }
    }
  }
  
  if (fit_data) {
    pfs ~ pch(right_uncensored, interval_censored, ignore_interval_censoring, disease_progress_prob, gen_pfs ? max_all_t : 0);
  }
}

generated quantities {
  vector[fit_data ? n_patients : 0] log_lik; 
  
  if (fit_data) {
    // Do not ignore censoring when calculating this
    log_lik = calc_pch_loglik(pfs, right_uncensored, interval_censored, 0, disease_progress_prob, gen_pfs ? max_all_t : 0);
  }
  
  vector<lower = 0, upper = 1>[max_all_t] base_pf_cond_prob = 1 - inv_cloglog(log_lambda); // Progression free conditional prob if not using covar
    // Progress free conditional probability if only 1 tumor per patient fixed at size = 1 
  vector<lower = 0, upper = 1>[max_all_t] one_tumor_pf_cond_prob; 
  
  array[add_trial_level ? n_trials : 0] vector<lower = 0, upper = 1>[max_all_t] trial_base_pf_cond_prob; 
  array[add_trial_level ? n_trials : 0] vector<lower = 0, upper = 1>[max_all_t] trial_one_tumor_pf_cond_prob; 
  
  if (add_trial_level) {
    for (s in 1:n_trials) {
      trial_base_pf_cond_prob[s] = 1 - inv_cloglog(log_trial_lambda[s]);
      
      if (tumor_hazard_type > 0) {
        trial_one_tumor_pf_cond_prob[s] = 
          1 - calculate_progress_linear_prob(
            n_patient_tumors, log_trial_lambda[s], tumor_stim_trial_intercept[s], tumor_stim_trial_coef[s], [rep_row_vector(1, n_covar_col)]
          ); 
      } else {
        trial_one_tumor_pf_cond_prob[s] = trial_base_pf_cond_prob[s];
      } 
    }
  }
 
  if (tumor_hazard_type > 0) {
    one_tumor_pf_cond_prob = 
      1 - calculate_progress_linear_prob(n_patient_tumors, log_lambda, tumor_stim_intercept, tumor_stim_coef, [rep_row_vector(1, n_covar_col)]); 
  } else {
    one_tumor_pf_cond_prob = base_pf_cond_prob;
  } 
  
  vector<lower = 0, upper = 1>[max_all_t] base_survival;
  vector<lower = 0, upper = 1>[max_all_t] one_tumor_survival;
  real<lower = 0, upper = max_all_t> base_cond_expected_pfs;
  real<lower = 0, upper = max_all_t> one_tumor_cond_expected_pfs;
  real<lower = 0, upper = max_all_t> base_cond_median_pfs;
  real<lower = 0, upper = max_all_t> one_tumor_cond_median_pfs;
  
  row_vector<lower = 0, upper = 1>[max_all_t] base_dp_prob;
  row_vector<lower = 0, upper = 1>[max_all_t] one_tumor_dp_prob;
  
  row_vector<lower = 0, upper = 1>[max_all_t] base_dp_prob_not_censored;
  row_vector<lower = 0, upper = 1>[max_all_t] one_tumor_dp_prob_not_censored;
  
  array[gen_pfs ? n_patients : 0] int<lower = 0> rep_pfs;
  array[gen_pfs ? n_patients : 0] int<lower = 0, upper = 1> rep_right_censored;
  array[gen_pfs ? n_patients : 0] int<lower = 0> rep_interval_censored;
 
  // Kaplan-Meier survival probability, aggregated over generated patients' data.  
  vector<lower = 0, upper = 1>[gen_pfs ? max(t_measure) + 1 : 0] km_est; 
  array[add_trial_level && gen_pfs ? n_trials : 0] vector<lower = 0, upper = 1>[max(t_measure) + 1] trial_km_est; 
  
  if (gen_pfs) {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
    int t_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + max_all_t - 1;
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int t_end = t_pos + n_measures[tumor_pos] - 1;
     
      tuple(int, int, int, int) pfs_res = pfs_rng(
        disease_progress_prob[pfs_interval_pos:pfs_interval_end], 
        gen_interval_censored ? t_measure[(t_pos + n_screening_t[tumor_pos]):t_end] : pfs_range_int
      );
      
      rep_interval_censored[i] = pfs_res.1;
      rep_right_censored[i] = pfs_res.2;
      rep_pfs[i] = pfs_res.3; 
      
      pfs_interval_pos = pfs_interval_end + 1;
      t_pos += sum(n_measures[tumor_pos:tumor_end]);
      tumor_pos = tumor_end + 1;
    }
    
    km_est = estimate_kaplan_meier(rep_pfs, rep_right_censored, max_all_t).1; 
    
    
    if (add_trial_level) {
      int trial_patient_pos = 1;
      
      for (s in 1:n_trials) {
        int trial_patient_end = trial_patient_pos + n_trial_patients[s] - 1;
        
        trial_km_est[s] = estimate_kaplan_meier(rep_pfs[trial_patient_pos:trial_patient_end], rep_right_censored[trial_patient_pos:trial_patient_end], max_all_t).1; 
        
        trial_patient_pos = trial_patient_end + 1;
      }
    }
  } 
  
  {
    tuple(vector[max_all_t], vector[max_all_t]) base_marginal_prob_res = calculate_marginal_dp_prob(base_pf_cond_prob, max_all_t);  
    tuple(vector[max_all_t], vector[max_all_t]) one_tumor_marginal_prob_res = calculate_marginal_dp_prob(one_tumor_pf_cond_prob, max_all_t);  
    
    base_dp_prob = base_marginal_prob_res.1';
    one_tumor_dp_prob = one_tumor_marginal_prob_res.1';
    
    base_survival = base_marginal_prob_res.2;
    one_tumor_survival = one_tumor_marginal_prob_res.2;
    
    base_dp_prob_not_censored = base_dp_prob / (1 - base_survival[max_all_t]);
    one_tumor_dp_prob_not_censored = one_tumor_dp_prob / (1 - one_tumor_survival[max_all_t]);
    
    base_cond_expected_pfs = base_dp_prob_not_censored * pfs_range_vec[:max_all_t];
    one_tumor_cond_expected_pfs = one_tumor_dp_prob_not_censored * pfs_range_vec[:max_all_t];
    
    base_cond_median_pfs = pfs_quantiles_from_prob(base_dp_prob_not_censored', { 0.5 })[1];
    one_tumor_cond_median_pfs = pfs_quantiles_from_prob(one_tumor_dp_prob_not_censored', { 0.5 })[1];
  }
}
