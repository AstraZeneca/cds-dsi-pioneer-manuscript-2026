functions {
  #include "util.stan"
}

data {
  int<lower = 0, upper = 1> fit_data; // Sample of prior prediction only
  int<lower = 0, upper = 1> predict_missing_sizes; 
  int<lower = 0, upper = 1> gen_tumor_sizes; // Generate data given prior/posterior of model parameters
  int<lower = 0, upper = 1> multilevel_patient;
  int<lower = 0, upper = 1> multilevel_tumor;
  
  #include "base_data.stan"
  
  // Hyperparam
  real<lower = 0> pop_tumor_gp_rho_alpha;
  real<lower = 0> pop_tumor_gp_rho_beta;
}

transformed data {
  #include "tumor_transformed_data.stan"
}

parameters {
  real<lower = 0> tumor_mean; // Actual tumor mean, not the lognormal mean
  real<lower = 0> tumor_sd; // Actual tumor SD, not the lognormal one
 
  real<lower = 0> pop_tumor_gp_rho;
 
  // Multilevel intercepts
  
  real<lower = 0> patient_tumor_gp_intercept_sd; 
  vector[multilevel_patient ? n_patients : 0] patient_tumor_gp_intercept_effect;
  
  vector<lower = 0>[multilevel_tumor ? n_patients : 0] tumor_gp_intercept_sd; 
  vector[multilevel_tumor ? sum(n_patient_tumors) : 0] tumor_gp_intercept_effect;
  
}

transformed parameters {
  // Given priors on actual tumor mean and SD calculate the corresponding lognormal ones
  real<lower = 0> pop_tumor_gp_alpha = sqrt(log((tumor_sd / tumor_mean)^2 + 1));
  real pop_tumor_gp_intercept = log(tumor_mean) - pop_tumor_gp_alpha^2 / 2.0;
  
  vector[sum(n_patient_tumors)] tumor_gp_intercept = rep_vector(pop_tumor_gp_intercept, sum(n_patient_tumors));
  
  if (multilevel_patient) {
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        tumor_gp_intercept[tumor_pos] += patient_tumor_gp_intercept_effect[i];
        
        tumor_pos += 1;
      }
    }
  }
  
  if (multilevel_tumor) {
    tumor_gp_intercept += tumor_gp_intercept_effect;
  }
}

model {
  tumor_mean ~ normal(2.8, 0.1);
  // tumor_sd ~ normal(0, 1.25);
  tumor_sd ~ normal(0, 2);
  
  patient_tumor_gp_intercept_sd ~ normal(0, 0.25);
 
  if (multilevel_patient) { 
    patient_tumor_gp_intercept_effect ~ normal(0, patient_tumor_gp_intercept_sd);
  }
  
  if (multilevel_tumor) {
    tumor_gp_intercept_sd ~ normal(0, 0.1);
    
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        tumor_gp_intercept_effect ~ normal(0, tumor_gp_intercept_sd[i]);
        
        tumor_pos += 1;
      }
    }
  }
  
  pop_tumor_gp_rho ~ inv_gamma(pop_tumor_gp_rho_alpha, pop_tumor_gp_rho_beta);

  if (fit_data) {
    int tumor_pos = 1;
    int t_measure_pos = 1;
    int t_missing_measure_pos = 1;

    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
        
        matrix[n_measures[tumor_pos], n_measures[tumor_pos]] L_current_tumor_vcov = 
          calc_gp_cholesky_vcov(all_measure_t[patient_t_measure_idx[t_measure_pos:t_measure_end]], pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
        
        log(tumor_size[t_measure_pos:t_measure_end]) ~ multi_normal_cholesky(rep_vector(tumor_gp_intercept[tumor_pos], n_measures[tumor_pos]), L_current_tumor_vcov);

        tumor_pos += 1;
        t_measure_pos = t_measure_end + 1;
      }
    }
  }  
}

generated quantities {
  vector<lower = 0>[predict_missing_sizes ? n_all_tumor_measures : 0] all_tumor_size; // This includes both observed and missing measurements
  vector<lower = 0>[gen_tumor_sizes ? sum(n_measures) : 0] rep_tumor_size; // Drawing *new* data
 
  if (predict_missing_sizes || gen_tumor_sizes) { 
    matrix[max(patient_max_t_width), max(patient_max_t_width)] pop_tumor_vcov = calc_gp_vcov(all_measure_t, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
    
    int tumor_pos = 1;
    int full_measure_pos = 1;
    int obs_measure_pos = 1;
    int missing_measure_pos = 1;
    
    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        int full_measure_end = full_measure_pos + patient_max_t_width[i] - 1;
        int obs_measure_end = obs_measure_pos + n_measures[tumor_pos] - 1;
        int missing_measure_end = missing_measure_pos + n_missing_measures[tumor_pos] - 1;
        
        array[n_measures[tumor_pos]] int current_t_measure_idx = patient_t_measure_idx[obs_measure_pos:obs_measure_end];
        
        if (predict_missing_sizes) {
          array[n_missing_measures[tumor_pos]] int current_t_missing_measure_idx = patient_t_missing_measure_idx[missing_measure_pos:missing_measure_end];
          
          vector[patient_max_t_width[i]] current_tumor_size;
          
          current_tumor_size[current_t_measure_idx] = tumor_size[obs_measure_pos:obs_measure_end];
          
          if (n_missing_measures[tumor_pos] > 0) {
            // Here we are not simply drawing tumor sizes from the GP distribution conditional on observed data. We need to make sure that the new
            // imputations are jointly determined with the observed data, i.e., are correlated to the observed data.
            
            current_tumor_size[current_t_missing_measure_idx] = 
              exp(gp_pred_rng( // This function I took from the Stan docs that makes it easy to make these draws.
                all_measure_t[current_t_missing_measure_idx],
                log(current_tumor_size[current_t_measure_idx]), 
                all_measure_t[current_t_measure_idx], 
                pop_tumor_vcov[current_t_measure_idx, current_t_measure_idx],
                pop_tumor_gp_alpha,
                pop_tumor_gp_rho,
                delta
              )); 
          }
          
          all_tumor_size[full_measure_pos:full_measure_end] = current_tumor_size;
        }
        
        if (gen_tumor_sizes) {
          rep_tumor_size[obs_measure_pos:obs_measure_end] = exp(multi_normal_rng(
            rep_vector(tumor_gp_intercept[tumor_pos], n_measures[tumor_pos]), 
            pop_tumor_vcov[current_t_measure_idx, current_t_measure_idx]
          ));
        }
        
        full_measure_pos = full_measure_end + 1;
        obs_measure_pos = obs_measure_end + 1;
        missing_measure_pos = missing_measure_end + 1;
        tumor_pos += 1;
      }
    }
  }
}
