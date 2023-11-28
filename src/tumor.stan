functions {
  #include "util.stan"
  
  // Scale tumor sizes by the standard deviation of all non-zero tumors (a size of zero means the tumor doesn't exist yet/anymore).
  tuple(real, vector) scale_nonzero_tumor_sizes(vector tumor_size) {
    int n_tumor_measures = rows(tumor_size);
    array[n_tumor_measures] int nonzero_tumor_idx;
    int measured_pos = 1;
    real tumor_sd;
    
    for (t in 1:n_tumor_measures) {
      if (tumor_size[t] > 0) {
        nonzero_tumor_idx[measured_pos] = t;
        measured_pos += 1;
      }
    }
    
    tumor_sd = sd(tumor_size[nonzero_tumor_idx[:(measured_pos - 1)]]);
    
    return (tumor_sd, tumor_size / tumor_sd); 
  }
}

data {
  int<lower = 0, upper = 1> fit_data; // Sample of prior prediction only
  int<lower = 0, upper = 1> gen_tumor_sizes; // Generate data given prior/posterior of model parameters
  
  #include "base_data.stan"
 
  // If generating data no input needed. 
  // [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
  vector<lower = 0>[gen_tumor_sizes ? 0 : to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures))] tumor_size; 
}

transformed data {
  array[sum(n_measures) - n_patients] int<lower = 2> t_measure_p1 = to_int(to_array_1d(to_vector(t_measure) + 1));  
  int<lower = 0> n_tumor_measures = to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures));
  int n_all_tumor_measures = sum(n_patient_tumors); // This represents all _potential_ measures.
  int<lower = 0> max_all_t = max(t_measure); 
  array[n_patients] int<lower = 0, upper = max_all_t> max_t;
  array[max_all_t + 1] real all_measure_idx; // This is for the GP "proximity" between size measurement time intervals.
  real delta = 1e-9;
  // real<lower = 0> tumor_sd = 0;
  vector<lower = 0>[rows(tumor_size) > 0 ? n_tumor_measures : 0] scaled_tumor_size;
  
  // Not using this yet. This would be used to identify which of all the potential measures are actually observed in the data.
  array[n_tumor_measures] int<lower = 1, upper = n_tumor_measures> obs_tumor_measures_idx; 
  
  if (rows(tumor_size) > 0) {
    tuple(real, vector[n_tumor_measures]) scale_results = scale_nonzero_tumor_sizes(tumor_size);
    
    // tumor_sd = scale_results.1;
    // scaled_tumor_size = scale_results.2;
    scaled_tumor_size = tumor_size; // BUGBUG not scaling yet
  }
 
  for (t in 1:(max_all_t + 1)) {
    all_measure_idx[t] = t;
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
  real pop_tumor_gp_intercept;
  real<lower = 0> pop_tumor_gp_alpha;
  real<lower = 0> pop_tumor_gp_rho;
}

transformed parameters {
  // GP vcov matrix
  matrix[max_all_t + 1, max_all_t + 1] L_tumor_vcov = calc_gp_cholesky_vcov(all_measure_idx, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
}

model {
  pop_tumor_gp_intercept ~ normal(0, 0.25);
  pop_tumor_gp_alpha ~ normal(0, 0.25);
  pop_tumor_gp_rho ~ inv_gamma(5, 5);

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

        log(scaled_tumor_size[tumor_pos:tumor_end]) ~ multi_normal_cholesky(
          rep_vector(pop_tumor_gp_intercept, max_t[i] + 1), 
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
