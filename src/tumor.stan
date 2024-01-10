functions {
  #include "util.stan"
}

data {
  int<lower = 0, upper = 1> fit_data; // Sample of prior prediction only
  int<lower = 0, upper = 1> predict_missing_sizes; 
  int<lower = 0, upper = 1> gen_tumor_sizes; // Generate data given prior/posterior of model parameters
  
  #include "base_data.stan"
 
  // If generating data no input needed. 
  // [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
  vector<lower = 0>[sum(n_measures)] tumor_size; // cm 
  
  // Hyperparam
  real<lower = 0> pop_tumor_gp_rho_alpha;
  real<lower = 0> pop_tumor_gp_rho_beta;
}

transformed data {
  real delta = 1e-9;
  
  int<lower = 0> n_tumor_measures = sum(n_measures); 
  int<lower = 0> n_all_tumor_measures;
  
  int min_all_t = min(t_measure);
  int<lower = min_all_t> max_all_t = max(t_measure);
  array[sum(n_measures)] int<lower = 1> patient_t_measure_idx;
  
  array[n_patients] int<lower = 0> patient_max_t_width;
  array[sum(n_patient_tumors)] int<lower = 0> n_missing_measures = calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors); 
  array[sum(n_missing_measures)] int<lower = 1> patient_t_missing_measure_idx; 
  
  {
    int tumor_pos = 1;
    int t_measure_pos = 1;
    
    array[sum(n_measures)] int t_measure_idx;
    
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      int save_t_measure_pos = t_measure_pos;
      
      array[n_patient_tumors[i]] int min_t_idx;
      array[n_patient_tumors[i]] int max_t_idx;
      
      for (j in 1:n_patient_tumors[i]) {
        int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1;
        
        for (tp in t_measure_pos:t_measure_end) {
          t_measure_idx[tp] = t_measure[tp] - min_all_t + 1;
        }
        
        min_t_idx[j] = min(t_measure_idx[t_measure_pos:t_measure_end]);
        max_t_idx[j] = max(t_measure_idx[t_measure_pos:t_measure_end]);
        
        t_measure_pos = t_measure_end + 1;
      }
      
      patient_max_t_width[i] = max(max_t_idx) - min(min_t_idx) + 1;
      
      t_measure_pos = save_t_measure_pos;
      
      for (j in 1:n_patient_tumors[i]) {
        int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1; 
        
        for (tp in t_measure_pos:t_measure_end) {
          patient_t_measure_idx[tp] = t_measure_idx[tp] - min(min_t_idx) + 1;
        }
        
        t_measure_pos = t_measure_end + 1;
      }
      
      tumor_pos = tumor_end + 1;
    }
  }
  
  array[max(patient_max_t_width)] real all_measure_t; // This is for the GP "proximity" between size measurement time intervals.
  
  for (t in 1:max(patient_max_t_width)) {
    all_measure_t[t] = t / 12.0; 
  } 
  
  n_all_tumor_measures = to_int(to_row_vector(patient_max_t_width) * to_vector(n_patient_tumors));
  
  if (sum(n_missing_measures) > 0) { 
    patient_t_missing_measure_idx = calculate_t_missing_measure(n_measures, n_missing_measures, patient_t_measure_idx, n_patient_tumors);
  }
}

parameters {
  real<lower = 0> tumor_mean;
  real<lower = 0> tumor_sd;
  
  real<lower = 0> pop_tumor_gp_rho;
}

transformed parameters {
  real<lower = 0> pop_tumor_gp_alpha = sqrt(log((tumor_sd / tumor_mean)^2 + 1));
  real pop_tumor_gp_intercept = log(tumor_mean) - pop_tumor_gp_alpha^2 / 2.0;
}

model {
  tumor_mean ~ normal(2.8, 0.1);
  tumor_sd ~ normal(0, 1.25);
  
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
        
        log(tumor_size[t_measure_pos:t_measure_end]) ~ multi_normal_cholesky(rep_vector(pop_tumor_gp_intercept, n_measures[tumor_pos]), L_current_tumor_vcov);

        tumor_pos += 1;
        t_measure_pos = t_measure_end + 1;
      }
    }
  }  
}

generated quantities {
  vector<lower = 0>[predict_missing_sizes ? n_all_tumor_measures : 0] all_tumor_size;
  vector<lower = 0>[gen_tumor_sizes ? sum(n_measures) : 0] rep_tumor_size;
  cov_matrix[max(patient_max_t_width)] pop_tumor_vcov = calc_gp_vcov(all_measure_t, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
  
  if (predict_missing_sizes) {
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
        array[n_missing_measures[tumor_pos]] int current_t_missing_measure_idx = patient_t_missing_measure_idx[missing_measure_pos:missing_measure_end];
        
        vector[patient_max_t_width[i]] current_tumor_size;
        
        current_tumor_size[current_t_measure_idx] = tumor_size[obs_measure_pos:obs_measure_end];
        
        if (n_missing_measures[tumor_pos] > 0) {
          current_tumor_size[current_t_missing_measure_idx] = 
            exp(gp_pred_rng(
              all_measure_t[current_t_missing_measure_idx],
              current_tumor_size[current_t_measure_idx], 
              all_measure_t[current_t_measure_idx], 
              pop_tumor_vcov[current_t_measure_idx, current_t_measure_idx],
              pop_tumor_gp_alpha,
              pop_tumor_gp_rho,
              delta
            )); 
        }
        
        all_tumor_size[full_measure_pos:full_measure_end] = current_tumor_size;
        
        full_measure_pos = full_measure_end + 1;
        obs_measure_pos = obs_measure_end + 1;
        missing_measure_pos = missing_measure_end + 1;
        tumor_pos += 1;
      }
    }
  }
 
  if (gen_tumor_sizes) {
    int tumor_pos = 1;
    int t_measure_pos = 1;
    
    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
        
        rep_tumor_size[t_measure_pos:t_measure_end] = exp(multi_normal_rng(
          rep_vector(pop_tumor_gp_intercept, n_measures[tumor_pos]), 
          pop_tumor_vcov[patient_t_measure_idx[t_measure_pos:t_measure_end], patient_t_measure_idx[t_measure_pos:t_measure_end]]
        ));
        
        tumor_pos += 1; 
        t_measure_pos = t_measure_end + 1;
      }
    }
  }
}
