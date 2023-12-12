functions {
  #include "util.stan"
}

data {
  int<lower = 0, upper = 1> fit_data; // Sample of prior prediction only
  int<lower = 0, upper = 1> gen_tumor_sizes; // Generate data given prior/posterior of model parameters
  
  #include "base_data.stan"
 
  // If generating data no input needed. 
  // [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
  vector<lower = 0>[gen_tumor_sizes ? 0 : to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures))] tumor_size; // cm 
  
  // Hyperparam
  real<lower = 0> pop_tumor_gp_rho_alpha;
  real<lower = 0> pop_tumor_gp_rho_beta;
}

transformed data {
  real delta = 1e-9;
  array[sum(n_measures) - n_patients] int<lower = 2> t_measure_p1 = to_int(to_array_1d(to_vector(t_measure) + 1));  
  int<lower = 0> n_tumor_measures = to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures));
  int n_all_tumor_measures = sum(n_patient_tumors); // This represents all _potential_ measures.
  int<lower = 0> max_all_t = max(t_measure); 
  array[n_patients] int<lower = 0, upper = max_all_t> max_t;
  array[max_all_t + 1] real all_measure_idx; // This is for the GP "proximity" between size measurement time intervals.
  
  // Not using this yet. This would be used to identify which of all the potential measures are actually observed in the data.
  array[n_tumor_measures] int<lower = 1, upper = n_tumor_measures> obs_tumor_measures_idx; 
  
  for (t in 1:(max_all_t + 1)) {
    all_measure_idx[t] = t / 12.0; // 1 = year
  } 
 
  {
    int measure_pos = 1;
    int obs_measure_idx_pos = 1;
    int measure_offset = 0;
    
    for (i in 1:n_patients) {
      int measure_end = measure_pos + n_measures[i] - 2; // Base measurement not included in t_measure
      int n_current_tumors = n_patient_tumors[i];
      
      max_t[i] = max(t_measure[measure_pos:measure_end]);
      n_all_tumor_measures += max_t[i] * n_current_tumors;
      
      for (j in 1:n_current_tumors) {
        int obs_measure_idx_end = obs_measure_idx_pos + n_measures[i] - 1;
        
        obs_tumor_measures_idx[obs_measure_idx_pos] = 1 + measure_offset;
        
        for (t in 1:(n_measures[i] - 1)) {
          obs_tumor_measures_idx[obs_measure_idx_pos + t] = t_measure[measure_pos + t - 1] + 1 + measure_offset; 
        }
        
        obs_measure_idx_pos = obs_measure_idx_end + 1;
        measure_offset += max_t[i] + 1;
      }
      
      measure_pos = measure_end + 1;
    }
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
  
  // GP vcov matrix
  matrix[max_all_t + 1, max_all_t + 1] L_tumor_vcov = 
    calc_gp_cholesky_vcov(all_measure_idx, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
}

model {
  tumor_mean ~ normal(2.8, 0.1);
  // tumor_sd ~ normal(1, 0.1);
  tumor_sd ~ normal(0, 1.25);
  
  pop_tumor_gp_rho ~ inv_gamma(pop_tumor_gp_rho_alpha, pop_tumor_gp_intercept);

  if (fit_data) {
    int tumor_pos = 1;
    int t_measure_pos = 1;

    for (i in 1:n_patients) {
      int n_current_tumors = n_patient_tumors[i];
      int t_measure_end = t_measure_pos + n_measures[i] - 2;

      for (j in 1:n_current_tumors) {
        int tumor_end = tumor_pos + n_measures[i] - 1;
        
        // I need to add the baseline measure here
        array[n_measures[i]] int curr_full_t_measure = append_array({ 1 }, t_measure_p1[t_measure_pos:t_measure_end]);

        log(tumor_size[tumor_pos:tumor_end]) ~ multi_normal_cholesky(
          rep_vector(pop_tumor_gp_intercept, n_measures[i]), 
          L_tumor_vcov[curr_full_t_measure, curr_full_t_measure]
        );

        tumor_pos = tumor_end + 1;
      }

      t_measure_pos = t_measure_end + 1;
    }
  }  
}

generated quantities {
  vector<lower = 0>[gen_tumor_sizes ? n_all_tumor_measures : 0] rep_tumor_size;
 
  if (gen_tumor_sizes) {
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        int tumor_end = tumor_pos + max_t[i];
        
        rep_tumor_size[tumor_pos:tumor_end] = exp(multi_normal_cholesky_rng(
          rep_vector(pop_tumor_gp_intercept, max_t[i] + 1), L_tumor_vcov[:(max_t[i] + 1), :(max_t[i] + 1)]
        ));
        
        tumor_pos = tumor_end + 1;
      }
    }
  }
}
