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
    int disease_progress = 0;
    int pfs = 0;
    
    while (pfs < n_prob && !bernoulli_rng(prob[pfs + 1])) {
      pfs += 1;
    }
    
    return pfs;
  }
}

data {
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_pfs;
  
  #include "base_data.stan"
  
  matrix<lower = 0>[sum(n_patient_tumors), n_measures] tumor_size;
  
  int<lower = 0> max_pfs;
  array[n_patients] int<lower = 0> pfs;
  array[n_patients] int<lower = 0, upper = 1> censored;
}

transformed data {
  real delta = 1e-9;
  int<lower = 0, upper = max_pfs * n_patients> n_total_pfs = sum(pfs);
  array[max_pfs] real pfs_range;
  int<lower = 0> n_all_tumors = sum(n_patient_tumors);
  matrix[n_all_tumors, n_measures] tumor_covar = tumor_size;
  array[n_patients] int<lower = 0, upper = 1> uncensored;
  
  for (i in 1:max_pfs) {
    pfs_range[i] = i;
  }
  
  if (n_measures != 2) {
    reject("Only supporting two measures for now.");
  }
  
  tumor_covar[, 2] -= tumor_covar[, 1];
  
  for (i in 1:n_patients) {
    uncensored[i] = 1 - censored[i];
  }
}

parameters {
  real<lower = 0> log_lambda_gp_alpha;
  real<lower = 0> log_lambda_gp_rho;
  vector[max_pfs] log_lambda_gp_eta;
  real log_lambda_gp_intercept;
  
  real<lower = 0> tumor_stim_intercept;
  vector<lower = 0>[2] tumor_stim_coef;
}

transformed parameters {
  vector[max_pfs] log_lambda = calc_gp_pred(pfs_range, log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, delta, log_lambda_gp_eta);
  vector[n_total_pfs + sum(uncensored)] disease_progress_pred;
  vector<lower = 0, upper = 1>[n_total_pfs + sum(uncensored)] disease_progress_prob;
  vector[n_patients] total_time_invar_tumor_stim;
  
  {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
  
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int pfs_interval_end = pfs_interval_pos + pfs[i] - censored[i]; 
      
      total_time_invar_tumor_stim[i] = sum(linear_tumor_stimulus(tumor_stim_intercept, tumor_stim_coef, tumor_covar[tumor_pos:tumor_end]));
      
      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = log_lambda[1:(pfs[i] + uncensored[i])] + total_time_invar_tumor_stim[i];
      
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
  log_lambda_gp_intercept ~ normal(-2, 0.5);
  
  tumor_stim_intercept ~ normal(0, 0.25);
  tumor_stim_coef[1] ~ normal(0, 0.125);
  tumor_stim_coef[2] ~ normal(0, 0.125);
  
  if (fit_data) {
    int pfs_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + pfs[i] - censored[i]; 
   
      if (pfs[i] > 0) {
        target += bernoulli_lupmf(0 | disease_progress_prob[pfs_interval_pos:(pfs_interval_end - uncensored[i])]);
      }

      target += uncensored[i] * bernoulli_lupmf(1 | disease_progress_prob[pfs_interval_end]);
      
      pfs_interval_pos = pfs_interval_end + 1;
    }
  }
}

generated quantities {
  vector<lower = 0, upper = 1>[max_pfs] base_pf_prob = 1 - inv_cloglog(log_lambda);
  
  array[gen_pfs ? n_patients : 0] int<lower = 0> rep_pfs;
  array[gen_pfs ? n_patients : 0] int<lower = 0, upper = 1> rep_censored;
  
  if (gen_pfs) {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + max_pfs - 1;
      int t = 1;
      int disease_progress = 0;
      
      rep_pfs[i] = pfs_rng(disease_progress_prob[pfs_interval_pos:pfs_interval_end]);
      rep_censored[i] = rep_pfs[i] >= max_pfs;
      
      pfs_interval_pos = pfs_interval_end + 1;
    }
  }
}
