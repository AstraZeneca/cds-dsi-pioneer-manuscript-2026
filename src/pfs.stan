functions {
  #include "util.stan"
 
  // Simple regression model for the influence of tumors on surival. 
  vector linear_tumor_stimulus(real intercept, vector coef, matrix covar) {
    return intercept + covar * coef; 
  } 
 
  // Combine influence of all tumors on survival and calculate probability of survival using a cloglog link function. 
  vector calculate_progress_linear_prob(vector log_lambda, real tumor_intercept, vector tumor_coef, matrix tumor_covar) {
      real total_time_invar_tumor_stim = sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar));
      
      return inv_cloglog(log_lambda + total_time_invar_tumor_stim);
  }
 
  // Hazard function given a base hazard and time-invariant covariates.  
  vector calculate_linear_hazard(vector log_lambda, real tumor_intercept, vector tumor_coef, matrix tumor_covar) {
      real total_time_invar_tumor_stim = sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar));
      
      return exp(log_lambda + total_time_invar_tumor_stim);
  }
 
  // Calculate the last observed measure for each patient. 
  array[] int get_max_t(array[] int t_measure, array[] int n_measures) {
    int n_patients = size(n_measures);
    int t_pos = 1;
    array[n_patients] int max_t;
      
    for (i in 1:n_patients) {
      int t_end = t_pos + n_measures[i] - 2; // The baseline measure is not included in t_measures.
      
      max_t[i] = max(t_measure[t_pos:t_end]);
    }
    
    return max_t;
  }
 
  // Given PFS and tumor measures data, determine interval and right censoring for each patient. 
  tuple(array[] int, array[] int) identify_censoring(array[] int pfs, array[] int n_measures, array[] int t_measure) { 
    int n_patients = size(n_measures);
    array[n_patients] int interval_censored = rep_array(0, n_patients);
    array[n_patients] int right_censored = rep_array(1, n_patients);
    int t_pos = 1;
      
      for (i in 1:n_patients) {
        int t_end = t_pos + n_measures[i] - 2; // The baseline measure is not included in t_measures.
        
        for (t_index in 1:(n_measures[i] - 1)) {
          int curr_t = t_measure[t_pos + t_index - 1]; 
          
          if (curr_t > pfs[i]) {
            interval_censored[i] = curr_t - pfs[i] - 1; // Progression actually happened between pfs[i] and the next measured interval. 
            right_censored[i] = 0; // Found an observation after pfs[i]
            
            break;
          }
        }
        
        t_pos = t_end + 1;
      }
      
      return (interval_censored, right_censored);
  }
 
  // Random PFS generator given conditional progression probability and obseration intervals. 
  tuple(int, int, int, int) pfs_rng(vector prob, array[] int t) {
    int n_prob = rows(prob);
    int n_t = size(t);
    int actual_pfs = 0;
    int observed_pfs = 0;
    int right_censored = 0;
    int interval_censored = 0;
    int pfs_measure_index = 0;
    
    if (n_prob < t[n_t]) {
      reject("Insufficient number of probabilities provided.");
    }
   
    // Iterate over intervals until disease progression occurs 
    while (actual_pfs < n_prob && !bernoulli_rng(prob[actual_pfs + 1])) {
      actual_pfs += 1;
    }
   
    // Find the index in the measurement t array that corresponds to the interval observed after last progression free interval.
    while (pfs_measure_index < n_t && t[pfs_measure_index + 1] <= actual_pfs) {
      pfs_measure_index += 1;
    }
    
    right_censored = pfs_measure_index >= n_t; // We passed beyond the measurement t array so we must be right censored.
    observed_pfs = pfs_measure_index > 0 ? t[min(pfs_measure_index, n_t)] : 0; // Figure out which of the observed intervals would be the observed PFS
    interval_censored = !right_censored ? t[pfs_measure_index + 1] - observed_pfs - 1 : 0; // Figure how many intervals forward could be the true PFS 
    
    return (interval_censored, right_censored, observed_pfs, actual_pfs);
  }
 
  // Survival aggregated over all patients, S(t) = Pr[T > t], t \in {0,..., N} 
  vector estimate_kaplan_meier(array[] int pfs, array[] int right_censored, int max_t) {
    vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
    int n_pfs = size(pfs); // How many patients
    array[n_pfs] int sorted_pfs_idx = sort_indices_asc(pfs);
    int pfs_pos = 1;
    int n = n_pfs; // How many patients still haven't seen disease progression. 
    
    // For each time interval in 0..max_t see how many patiented exited and calculate proportion surviving.
    
    for (t in 0:max_t) {
      int ex = 0; // How many saw disease progression (exited) in current interval.
      real prev_s = t > 0 ? s[t] : 1.0;
      
      while ((n > 0) && (pfs_pos <= n_pfs) && (right_censored[sorted_pfs_idx[pfs_pos]] || (pfs[sorted_pfs_idx[pfs_pos]] <= t))) {
        ex += !right_censored[sorted_pfs_idx[pfs_pos]];
        pfs_pos += 1;
      }
     
      s[t + 1] = n > 0 ? prev_s * (n - ex) / n : 0.0;
      n -= ex;
    }
    
    return s; 
  }  

  // Calculate marginal probability of disease progression at every time interval, given conditional probabilities.
  tuple(vector, vector) calculate_marginal_dp_prob(vector cond_pf_prob, int max_all_t) {
    vector[max_all_t] marginal_dp_prob;
    vector[max_all_t] dp_cdf;
  
    for (m in 1:max_all_t) {
      if (m > 1) {
        marginal_dp_prob[m] = (1 - cond_pf_prob[m]) * prod(cond_pf_prob[1:(m - 1)]);
        dp_cdf[m] = marginal_dp_prob[m] + dp_cdf[m - 1]; 
      } else {
        marginal_dp_prob[m] = 1 - cond_pf_prob[m];
        dp_cdf[m] = marginal_dp_prob[m];
      }
    }
    
    return(marginal_dp_prob, 1 - dp_cdf);
  }  
 
  // Create (n_patients * n_tumors) x 2 matrix of each tumor's covariates from t = 1, 2. 
  matrix prepare_early_tumors_design_matrix(vector tumor_size, array[] int n_patient_tumors, array[] int n_measures) {
    matrix[sum(n_patient_tumors), 2] tumor_covar;
    int n_patients = size(n_patient_tumors);
    int tumor_pos = 1;
    int covar_pos = 1;
    
    for (i in 1:n_patients) {
      int n_current_tumors = n_patient_tumors[i];
      int n_current_measures = n_measures[i];
      
      if (n_current_measures < 2) {
        reject("Two measures minimum needed per tumor.");
      }
      
      for (j in 1:n_current_tumors) {
        int tumor_end = tumor_pos + n_current_measures - 1; 
        
        tumor_covar[covar_pos, ] = tumor_size[tumor_pos:(tumor_pos + 1)]';
        covar_pos += 1;
        
        tumor_pos = tumor_end + 1;
      }
    }
    
    return tumor_covar;
  }
}

data {
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_pfs;
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> early_tumors_only; // Only time-invariant covariates used: tumor size from t = 1,2.
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  
  #include "base_data.stan"
  
  // [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
  vector<lower = 0>[to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures))] tumor_size;
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive
  
  // Hyperparam
  
  real log_lambda_gp_intercept_mean;
  real<lower = 0> log_lambda_gp_intercept_sd;
  real<lower = 0> tumor_stim_intercept_sd; 
  vector<lower = 0>[2] tumor_stim_coef_sd;
}

transformed data {
  real delta = 1e-9;
  int<lower = 0> n_total_measures = to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures));
  int<lower = 0> max_pfs = max(pfs);
  array[n_patients] int<lower = 1> max_t = get_max_t(t_measure, n_measures);
  int<lower = 1> max_all_t = max(max_t);
  int<lower = 0, upper = max_pfs * n_patients> n_total_pfs = sum(pfs);
  int<lower = 0> n_time_periods;
  
  // These are used for the time interval distance between base hazard
  array[max_all_t] real pfs_range;
  array[max_all_t] int pfs_range_int;
  vector[max_all_t] pfs_range_vec;
  
  int<lower = 0> n_all_tumors = sum(n_patient_tumors);
  matrix[early_tumors_only ? n_all_tumors : 0, 2] standardized_tumor_covar;
 
  // Censoring information calculated from PFS and t_measure; no need to pass it in. 
  array[n_patients] int<lower = 0> interval_censored = rep_array(0, n_patients);
  array[n_patients] int<lower = 0, upper = 1> right_censored = rep_array(0, n_patients);
  array[n_patients] int<lower = 0, upper = 1> right_uncensored = rep_array(0, n_patients);
  
  for (i in 1:max_all_t) {
    pfs_range[i] = i;
    pfs_range_int[i] = i;
  }
  
  pfs_range_vec = to_vector(pfs_range);
 
  if (early_tumors_only) {
    tuple(real, real, vector[n_total_measures]) standardize_results = standardize_nonzero_tumor_sizes(tumor_size);
    vector[n_total_measures] standardized_tumor_size = standardize_results.3;
    
    standardized_tumor_covar = prepare_early_tumors_design_matrix(standardized_tumor_size, n_patient_tumors, n_measures);
  } else {
    reject("Not supported yet.");
  }
  
  if (fit_data) {
    tuple(array[n_patients] int, array[n_patients] int) censoring_res = identify_censoring(pfs, n_measures, t_measure);
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
  n_time_periods = gen_pfs ? n_patients * max_all_t : n_total_pfs + sum(right_uncensored) + sum(interval_censored);
}

parameters {
  // Base hazard GP parameters
  real<lower = 0> log_lambda_gp_alpha;
  real<lower = 0> log_lambda_gp_rho;
  vector[max_all_t] log_lambda_gp_eta;
  real log_lambda_gp_intercept;
 
  // Tumor level influence on hazard 
  real<lower = 0> tumor_stim_intercept;
  vector<lower = 0>[2] tumor_stim_coef;
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
      
      total_time_invar_tumor_stim[i] = sum(linear_tumor_stimulus(
        tumor_stim_intercept, tumor_stim_coef, standardized_tumor_covar[tumor_pos:tumor_end]
      ));
      
      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = 
        log_lambda[1:(gen_pfs ? max_t[i] : pfs[i] + right_uncensored[i] + interval_censored[i])] + total_time_invar_tumor_stim[i];
      
      tumor_pos = tumor_end + 1;
      pfs_interval_pos = pfs_interval_end + 1;
    }
      
    disease_progress_prob = inv_cloglog(disease_progress_pred); 
  }
}

model {
  // Priors
  
  log_lambda_gp_alpha ~ normal(0, 0.25);
  log_lambda_gp_rho ~ inv_gamma(5, 5);
  log_lambda_gp_eta ~ std_normal();
  log_lambda_gp_intercept ~ normal(log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd);
  
  tumor_stim_intercept ~ normal(0, tumor_stim_intercept_sd); 
  tumor_stim_coef[1] ~ normal(0, tumor_stim_coef_sd[1]);
  tumor_stim_coef[2] ~ normal(0, tumor_stim_coef_sd[2]);
  
  if (fit_data) {
    int pfs_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + pfs[i] + right_uncensored[i] + interval_censored[i] - 1; 
      int curr_interval_censored = ignore_interval_censoring ? 0 : interval_censored[i];
      vector[curr_interval_censored + right_uncensored[i]] interval_lp = rep_vector(0, curr_interval_censored + right_uncensored[i]);
      int observed_pfs_interval_end = pfs_interval_end - right_uncensored[i] - curr_interval_censored;
    
      // These are the time intervals we are sure that the patient has progression free 
      target += bernoulli_lupmf(0 | disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end]);
      
      for (j in 1:(curr_interval_censored + right_uncensored[i])) {
        if (j > 1) { // We need to add more possible intervals that the patient remained progression free.
          interval_lp[j] = bernoulli_lupmf(0 | disease_progress_prob[(observed_pfs_interval_end + 1):(observed_pfs_interval_end + j - 1)]);
        }
        
        // If not right censored add pdf of disease progression. 
        interval_lp[j] += right_uncensored[i] * bernoulli_lupmf(1 | disease_progress_prob[observed_pfs_interval_end + j]);
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
  vector<lower = 0, upper = 1>[max_all_t] one_tumor_pf_cond_prob = 
    // Progress free conditional probability if only 1 tumor per patient fixed at size = 1 
    1 - calculate_progress_linear_prob(log_lambda, tumor_stim_intercept, tumor_stim_coef, [[1, 1]]); 
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
      int t_end = t_pos + n_measures[i] - 2; // Baseline measure not included in t_measures
     
      tuple(int, int, int, int) pfs_res = pfs_rng(
        disease_progress_prob[pfs_interval_pos:pfs_interval_end], 
        gen_interval_censored ? t_measure[t_pos:t_end] : pfs_range_int
      );
      
      rep_interval_censored[i] = pfs_res.1;
      rep_right_censored[i] = pfs_res.2;
      rep_pfs[i] = pfs_res.3; 
      
      pfs_interval_pos = pfs_interval_end + 1;
      t_pos = t_end + 1;
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
