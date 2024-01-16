functions {
  #include "util.stan"
  #include "pfs_functions.stan"
}

data {
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_pfs;
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> no_tumor_stim; // Don't use tumor characteristics to predict survival
  int<lower = 0, upper = 1> early_tumors_only; // Only time-invariant covariates used: tumor size from t = 1,2.
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  
  #include "base_data.stan"
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive
  array[n_patients] int<lower = 0> death_week; 
  
  // Hyperparam
  
  real log_lambda_gp_intercept_mean;
  real<lower = 0> log_lambda_gp_intercept_sd;
  real<lower = 0> log_lambda_gp_alpha_sd;
  real<lower = 0> log_lambda_gp_rho_alpha;
  real<lower = 0> log_lambda_gp_rho_beta;
  real<lower = 0> tumor_stim_intercept_sd;
  vector<lower = 0>[2] tumor_stim_coef_sd;
}

transformed data {
  #include "tumor_transformed_data.stan" 
  
  if (min(n_patient_screening_t) <= 0) {
    reject("We need at least one screening measurement per patient, for now: ", min(n_patient_screening_t));
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
  
  matrix[early_tumors_only && !no_tumor_stim ? sum(n_patient_tumors) : 0, 2] standardized_tumor_covar;
  
  if (!no_tumor_stim) { 
    if (early_tumors_only) {
      tuple(real, real, vector[sum(n_measures)]) standardize_results = standardize_nonzero_tumor_sizes(tumor_size);
      vector[sum(n_measures)] standardized_tumor_size = standardize_results.3;
      
      standardized_tumor_covar = prepare_early_tumors_design_matrix(standardized_tumor_size, n_patient_tumors, n_measures, n_screening_t);
    } else {
      reject("Not supported yet.");
    }
  }
 
  // Censoring information calculated from PFS and t_measure; no need to pass it in. 
  array[n_patients] int<lower = 0> interval_censored = rep_array(0, n_patients);
  array[n_patients] int<lower = 0, upper = 1> right_censored = rep_array(0, n_patients);
  array[n_patients] int<lower = 0, upper = 1> right_uncensored = rep_array(0, n_patients);

  if (fit_data) {
    tuple(array[n_patients] int, array[n_patients] int) censoring_res = identify_censoring(pfs, death_week, n_patient_tumors, n_measures, t_measure);
    interval_censored = censoring_res.1;
    right_censored = censoring_res.2;
    
    for (i in 1:n_patients) {
        right_uncensored[i] = 1 - right_censored[i];
        
        if (interval_censored[i] > 0 && right_censored[i]) {
          reject("Cannot be both interval and right censored.");
        }
    }
      
    print("Number of interval censored observations: ", sum(interval_censored));
  }
 
  // If generating PFS we need to calculate probs for all possible time intervals, otherwise only up to observed PFS. 
  int<lower = 0> n_time_periods = gen_pfs ? n_patients * max_all_t : n_total_pfs + sum(right_uncensored) + sum(interval_censored);
}

parameters {
  // Base hazard GP parameters
  real<lower = 0> log_lambda_gp_alpha;
  real<lower = 0> log_lambda_gp_rho;
  vector[max_all_t] log_lambda_gp_eta;
  real log_lambda_gp_intercept;
 
  // Tumor level influence on hazard 
  real<lower = 0> tumor_stim_intercept; // DO NOT REMOVE; this is a per tumor intercept and not per patient intercept which is included in lambda.
  vector<lower = 0>[no_tumor_stim ? 0 : 2] tumor_stim_coef;
}

transformed parameters {
  // Base hazard
  vector[max_all_t] log_lambda = calc_gp_pred(pfs_range, log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, delta, log_lambda_gp_eta);
  vector<lower = 0, upper = 1>[n_time_periods] disease_progress_prob;
  
  vector[n_time_periods] disease_progress_pred;
  
  { // Calculate patient-interval conditional probability of disease progression.
    vector[n_patients] total_time_invar_tumor_stim;
    
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
  
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int pfs_interval_end = pfs_interval_pos + (gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i]) - 1; 
      
      total_time_invar_tumor_stim[i] = 
        no_tumor_stim ? 0 : sum(linear_tumor_stimulus(tumor_stim_intercept, tumor_stim_coef, standardized_tumor_covar[tumor_pos:tumor_end]));
      
      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = 
        log_lambda[1:(gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i])] + total_time_invar_tumor_stim[i];
      
      tumor_pos = tumor_end + 1;
      pfs_interval_pos = pfs_interval_end + 1;
    }
      
    disease_progress_prob = inv_cloglog(disease_progress_pred); 
  }
}

model {
  // Priors
  
  log_lambda_gp_alpha ~ normal(0, log_lambda_gp_alpha_sd);
  log_lambda_gp_rho ~ inv_gamma(log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);
  log_lambda_gp_eta ~ std_normal();
  log_lambda_gp_intercept ~ normal(log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd);
  
  tumor_stim_intercept ~ normal(0, tumor_stim_intercept_sd);
  
  if (!no_tumor_stim) {
    tumor_stim_coef[1] ~ normal(0, tumor_stim_coef_sd[1]);
    tumor_stim_coef[2] ~ normal(0, tumor_stim_coef_sd[2]);
  }
  
  if (fit_data) {
    int pfs_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + pfs[i] + right_uncensored[i] + interval_censored[i] - 1; 
      int curr_interval_censored = ignore_interval_censoring ? 0 : interval_censored[i];
      vector[curr_interval_censored + right_uncensored[i]] interval_lp = rep_vector(0, curr_interval_censored + right_uncensored[i]);
      int observed_pfs_interval_end = pfs_interval_end - right_uncensored[i] - curr_interval_censored;
    
      // print("i = ", i); 
      // print("pfs[i] = ", pfs[i]); 
      // print("pfs_interval_pos = ", pfs_interval_pos);
      // print("pfs_interval_end = ", pfs_interval_end);
      // print("curr_interval_censored = ", curr_interval_censored);
      // print("observed_pfs_interval_end = ", observed_pfs_interval_end);
      // print("right_uncensored[i] = ", right_uncensored[i]);
      // print("disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end] = ", disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end]);
    
      // These are the time intervals we are sure that the patient has progression free 
      target += bernoulli_lupmf(0 | disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end]);
      
      for (t in 1:(curr_interval_censored + right_uncensored[i])) {
        if (t > 1) { // We need to add more possible intervals that the patient remained progression free.
          interval_lp[t] = bernoulli_lupmf(0 | disease_progress_prob[(observed_pfs_interval_end + 1):(observed_pfs_interval_end + t - 1)]);
        }
        
        // If not right censored add pdf of disease progression. 
        interval_lp[t] += right_uncensored[i] * bernoulli_lupmf(1 | disease_progress_prob[observed_pfs_interval_end + t]);
      }
     
      if (curr_interval_censored > 0) {
        // There are more than one candidate true PFS: sum of the probabilities and then log.
        target += log_sum_exp(interval_lp); 
      } else if (right_uncensored[i]) {
        target += interval_lp[1]; // PFS not observed because of right censoring.
      }
      
      // If generating PFS jump ahead to the beginning of the next patient's probs.
      pfs_interval_pos = (gen_pfs ? pfs_interval_pos + max_all_t - 1 : pfs_interval_end) + 1;
    }
  }
}

generated quantities {
  vector<lower = 0, upper = 1>[max_all_t] base_pf_cond_prob = 1 - inv_cloglog(log_lambda); // Progression free conditional prob if not using covar
    // Progress free conditional probability if only 1 tumor per patient fixed at size = 1 
  vector<lower = 0, upper = 1>[max_all_t] one_tumor_pf_cond_prob = no_tumor_stim ? 
    base_pf_cond_prob : 1 - calculate_progress_linear_prob(log_lambda, no_tumor_stim ? 0 : tumor_stim_intercept, tumor_stim_coef, [[1, 1]]);
  
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
  
  if (gen_pfs) {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
    int t_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + max_all_t - 1;
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int t_end = t_pos + n_measures[tumor_pos] - 1;
     
      // print("n_screening_t[i] = ", n_screening_t[i]); 
      // print("t_measure[(t_pos + n_screening_t[tumor_pos]):t_end] = ", t_measure[(t_pos + n_screening_t[tumor_pos]):t_end]);
      
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
    
    km_est = estimate_kaplan_meier(rep_pfs, rep_right_censored, max_all_t); 
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
