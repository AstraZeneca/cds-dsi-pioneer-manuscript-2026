functions {
  #include "util.stan"
  
  vector linear_tumor_stimulus(real intercept, vector coef, matrix covar) {
    return intercept + covar * coef; 
  } 
  
  vector calculate_progress_linear_prob(vector log_lambda, real tumor_intercept, vector tumor_coef, matrix tumor_covar) {
      real total_time_invar_tumor_stim = sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar));
      
      return inv_cloglog(log_lambda + total_time_invar_tumor_stim);
  }
  
  vector calculate_linear_hazard(vector log_lambda, real tumor_intercept, vector tumor_coef, matrix tumor_covar) {
      real total_time_invar_tumor_stim = sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar));
      
      return exp(log_lambda + total_time_invar_tumor_stim);
  }
  
  int pfs_rng(vector prob) {
    int n_prob = rows(prob);
    int pfs = 0;
    
    while (pfs < n_prob && !bernoulli_rng(prob[pfs + 1])) {
      pfs += 1;
    }
    
    return pfs;
  }
  
  vector estimate_kaplan_meier(array[] int pfs, int max_pfs) {
    vector[max_pfs + 1] s = rep_vector(1.0, max_pfs + 1);
    
    for (r in 1:max_pfs) {
      int n = 0; 
      int ex = 0;
      
      if (s[r] > 0) {
        for (i in 1:size(pfs)) {
          n += (pfs[i] + 1 >= r);
          ex += (pfs[i] + 1 == r);
        }
      }
      
      s[r + 1] = n > 0 ? s[r] * (n - ex) / n : 0;
    }
    
    return s[2:]; 
  }  
  
  matrix prepare_early_tumors_design_matrix(vector tumor_size, array[] int n_patient_tumors, array[] int n_measures, real tumor_sd) {
    matrix[sum(n_patient_tumors), 2] scaled_tumor_covar;
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
        
        scaled_tumor_covar[covar_pos, ] = tumor_size[tumor_pos:(tumor_pos + 1)]' / tumor_sd;
        covar_pos += 1;
        
        tumor_pos = tumor_end + 1;
      }
    }
    
    return scaled_tumor_covar;
  }
}

data {
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_pfs;
  int<lower = 0, upper = 1> early_tumors_only;
  
  #include "base_data.stan"
 
  // BUGBUG temporary; need to use the same shape as the tumor model. 
  // matrix<lower = 0>[sum(n_patient_tumors), 2] tumor_size;
  
  // [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
  vector<lower = 0>[to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures))] tumor_size; // cm
  
  array[n_patients] int<lower = 0> pfs;
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  array[n_patients] int<lower = 0> interval_censored;
  
  // Hyperparam
  
  real log_lambda_gp_intercept_mean;
  real<lower = 0> log_lambda_gp_intercept_sd;
  real<lower = 0> tumor_stim_intercept_sd; 
  vector<lower = 0>[2] tumor_stim_coef_sd;
}

transformed data {
  real delta = 1e-9;
  int<lower = 0> max_pfs = max(pfs);
  int<lower = 0, upper = max_pfs * n_patients> n_total_pfs = sum(pfs);
  int<lower = 0, upper = max_pfs * n_patients> n_time_periods;
  array[max_pfs + 1] real pfs_range;
  vector[max_pfs + 1] pfs_range_vec;
  int<lower = 0> n_all_tumors = sum(n_patient_tumors);
  real<lower = 0> tumor_covar_sd;
  matrix[early_tumors_only ? n_all_tumors : 0, 2] scaled_tumor_covar;
  array[n_patients] int<lower = 0, upper = 1> right_uncensored;
  
  for (i in 1:(max_pfs + 1)) {
    pfs_range[i] = i;
  }
  
  pfs_range_vec = to_vector(pfs_range);
 
  if (early_tumors_only) { 
    tumor_covar_sd = sd(tumor_size);
    
    scaled_tumor_covar = prepare_early_tumors_design_matrix(tumor_size, n_patient_tumors, n_measures, tumor_covar_sd); 
  } else {
    reject("Not supported yet.");
  }
  
  for (i in 1:n_patients) {
    right_uncensored[i] = 1 - right_censored[i];
    
    if (pfs[i] + interval_censored[i] > max_pfs) {
      reject("Interval censoring cannot exceed max_pfs.");
    }
    
    if (interval_censored[i] > 0 && right_censored[i]) {
      reject("Cannot be both interval and right censored.");
    }
  }
  
  n_time_periods = gen_pfs ? max_pfs * n_patients : n_total_pfs + sum(right_uncensored) + sum(interval_censored);
}

parameters {
  real<lower = 0> log_lambda_gp_alpha;
  real<lower = 0> log_lambda_gp_rho;
  vector[max_pfs + 1] log_lambda_gp_eta;
  real log_lambda_gp_intercept;
  
  real<lower = 0> tumor_stim_intercept;
  vector<lower = 0>[2] tumor_stim_coef;
}

transformed parameters {
  vector[max_pfs + 1] log_lambda = calc_gp_pred(
    pfs_range, log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, delta, log_lambda_gp_eta);
  vector[n_time_periods] disease_progress_pred;
  vector<lower = 0, upper = 1>[n_time_periods] disease_progress_prob;
  vector[n_patients] total_time_invar_tumor_stim;
  
  {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
  
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int pfs_interval_end = pfs_interval_pos + (gen_pfs ? max_pfs - 1 : pfs[i] - right_censored[i] + interval_censored[i]); 
      
      total_time_invar_tumor_stim[i] = sum(linear_tumor_stimulus(tumor_stim_intercept, tumor_stim_coef, scaled_tumor_covar[tumor_pos:tumor_end]));
      
      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = 
        log_lambda[1:(gen_pfs ? max_pfs : pfs[i] + right_uncensored[i] + interval_censored[i])] + total_time_invar_tumor_stim[i];
      
      tumor_pos = tumor_end + 1;
      pfs_interval_pos = pfs_interval_end + 1;
    }
      
    disease_progress_prob = inv_cloglog(disease_progress_pred); 
  }
}

model {
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
      int pfs_interval_end = pfs_interval_pos + pfs[i] - right_censored[i]; 
      vector[interval_censored[i] + 1] interval_lp = rep_vector(0, interval_censored[i] + 1);
      
      for (j in 0:(interval_censored[i])) {
        int candidate_pfs = pfs[i] + j;
        
        if (candidate_pfs > 0) {
          interval_lp[j + 1] = bernoulli_lupmf(0 | disease_progress_prob[pfs_interval_pos:(pfs_interval_end - right_uncensored[i])]);
        }
        
        interval_lp[j + 1] += right_uncensored[i] * bernoulli_lupmf(1 | disease_progress_prob[pfs_interval_end]);
        
        pfs_interval_end += 1;
      }
     
      if (interval_censored[i] > 0) {
        target += log_sum_exp(interval_lp);
      } else {
        target += interval_lp[1];
      }
      
      pfs_interval_pos = gen_pfs ? pfs_interval_pos + max_pfs : pfs_interval_end;
    }
  }
}

generated quantities {
  vector<lower = 0, upper = 1>[max_pfs + 1] base_pf_cond_prob = 1 - inv_cloglog(log_lambda);
  vector<lower = 0, upper = 1>[max_pfs + 1] one_tumor_pf_cond_prob = 1 - calculate_progress_linear_prob(log_lambda, tumor_stim_intercept, tumor_stim_coef, [[1, 1]]); 
  vector<lower = 0, upper = 1>[max_pfs] base_survival;
  vector<lower = 0, upper = 1>[max_pfs] one_tumor_survival;
  real<lower = 0, upper = max_pfs> base_cond_expected_pfs;
  real<lower = 0, upper = max_pfs> one_tumor_cond_expected_pfs;
  
  array[gen_pfs ? n_patients : 0] int<lower = 0> rep_pfs;
  array[gen_pfs ? n_patients : 0] int<lower = 0, upper = 1> rep_right_censored;
  
  vector<lower = 0, upper = 1>[gen_pfs ? max_pfs : 0] km_est; 
  
  
  if (gen_pfs) {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + max_pfs - 1;
      int t = 1;
      int disease_progress = 0;
      
      rep_pfs[i] = pfs_rng(disease_progress_prob[pfs_interval_pos:pfs_interval_end]);
      rep_right_censored[i] = rep_pfs[i] >= max_pfs;
      
      pfs_interval_pos = pfs_interval_end + 1;
    }
    
    km_est = estimate_kaplan_meier(rep_pfs, max_pfs); 
  } 
  
  {
    row_vector[max_pfs] base_dp_prob;
    row_vector[max_pfs] one_tumor_dp_prob;
  
    for (m in 1:max_pfs) {
      if (m > 1) {
        base_dp_prob[m] = (1 - base_pf_cond_prob[m]) * prod(base_pf_cond_prob[1:(m - 1)]);
        base_survival[m] = base_dp_prob[m] + base_survival[m - 1]; 
        one_tumor_dp_prob[m] = (1 - one_tumor_pf_cond_prob[m]) * prod(one_tumor_pf_cond_prob[1:(m - 1)]);
        one_tumor_survival[m] = one_tumor_dp_prob[m] + one_tumor_survival[m - 1]; 
      } else {
        base_dp_prob[m] = 1 - base_pf_cond_prob[m];
        base_survival[m] = base_dp_prob[m];
        one_tumor_dp_prob[m] = 1 - one_tumor_pf_cond_prob[m];
        one_tumor_survival[m] = one_tumor_dp_prob[m];
      }
    }
    
    base_survival = 1 - base_survival;
    one_tumor_survival = 1 - one_tumor_survival;
    
    base_cond_expected_pfs = (base_dp_prob / (1 - base_survival[max_pfs])) * pfs_range_vec[:max_pfs];
    one_tumor_cond_expected_pfs = (one_tumor_dp_prob / (1 - one_tumor_survival[max_pfs])) * pfs_range_vec[:max_pfs];
  }
}
