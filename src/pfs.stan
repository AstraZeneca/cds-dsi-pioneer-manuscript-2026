functions {
  #include "util.stan"
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
  
  real tumor_stim_intercept;
  vector<lower = 0>[2] tumor_stim_coef;
}

transformed parameters {
  vector[max_pfs] log_lambda = calc_gp_pred(pfs_range, log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, delta, log_lambda_gp_eta);
  vector[n_total_pfs + sum(uncensored)] disease_progress_pred;
  vector<lower = 0, upper = 1>[n_total_pfs + sum(uncensored)] disease_progress_prob;
  
  {
    int tumor_pos = 1;
    int pfs_interval_pos = 1;
  
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      real total_time_invar_tumor_stim = sum(tumor_stim_intercept + tumor_covar[tumor_pos:tumor_end] * tumor_stim_coef);
     
      int pfs_interval_end = pfs_interval_pos + pfs[i] - censored[i]; 
      
      disease_progress_pred[pfs_interval_pos:pfs_interval_end] = log_lambda[1:(pfs[i] + uncensored[i])] + total_time_invar_tumor_stim;
      
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
  log_lambda_gp_intercept ~ normal(0, 0.5);
  
  tumor_stim_intercept ~ normal(-0.5, 0.125);
  tumor_stim_coef ~ normal(0, 0.05);
  
  if (fit_data) {
    int pfs_interval_pos = 1;
    
    for (i in 1:n_patients) {
      int pfs_interval_end = pfs_interval_pos + pfs[i] - censored[i]; 
      
      target += bernoulli_lupmf(0 | disease_progress_prob[pfs_interval_pos:(pfs_interval_end - 1)]) + (censored[i] ? 0 : bernoulli_lupmf(1 | disease_progress_prob[pfs_interval_end]));
      
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
      int disease_progress = 0;
      rep_censored[i] = 0;
      rep_pfs[i] = 0;
      
      for (t in 1:max_pfs) {
        disease_progress = bernoulli_rng(disease_progress_prob[pfs_interval_pos]);
        
        pfs_interval_pos += 1;
        
        if (disease_progress) {
          rep_pfs[i] = t - 1;
          rep_censored[i] = t >= max_pfs;
          
          break;
        }
      }
    }
  }
}
