functions {
  #include "util.stan"
  
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
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_tumor_sizes;
  
  #include "base_data.stan"
  
  // [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
  vector<lower = 0>[gen_tumor_sizes ? 0 : to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures))] tumor_size; // cm
}

transformed data {
  int<lower = 0> n_tumor_measures = to_int(to_row_vector(n_patient_tumors) * to_vector(n_measures));
  int<lower = 0, upper = sum(n_patient_tumors) * max(t_measure)> n_all_tumor_measures = 0;
  array[n_patients] int<lower = 0, upper = max(t_measure)> max_t;
  array[max(t_measure)] real all_measure_idx;
  real delta = 1e-9;
  real<lower = 0> tumor_sd = 0;
  vector<lower = 0>[rows(tumor_size) > 0 ? n_tumor_measures : 0] scaled_tumor_size;
  array[n_tumor_measures] int<lower = 1, upper = n_tumor_measures> obs_tumor_measures_idx;
  
  if (rows(tumor_size) > 0) {
    tuple(real, vector[n_tumor_measures]) scale_results = scale_nonzero_tumor_sizes(tumor_size);
    
    tumor_sd = scale_results.1;
    scaled_tumor_size = scale_results.2;
  }
 
  for (t in 1:max(t_measure)) {
    all_measure_idx[t] = t;
  } 
 
  {
    int measure_pos = 1;
    int obs_measure_idx_pos = 1;
    int tumor_count = 0;
    
    for (i in 1:n_patients) {
      int measure_end = measure_pos + n_measures[i] - 1;
      int n_current_tumors = n_patient_tumors[i];
      
      max_t[i] = t_measure[measure_end];
      n_all_tumor_measures += max_t[i] * n_current_tumors;
      
      for (j in 1:n_current_tumors) {
        int obs_measure_idx_end = obs_measure_idx_pos + n_measures[i] - 1;
        
        for (t in 1:n_measures[i]) {
          obs_tumor_measures_idx[obs_measure_idx_pos + t - 1] = t_measure[measure_pos + t - 1] + tumor_count; 
        }
        
        obs_measure_idx_pos = obs_measure_idx_end + 1;
        tumor_count += 1;
      }
      
      measure_pos = measure_end + 1;
    }
  }
}

parameters {
  real pop_tumor_gp_intercept;
  real<lower = 0> pop_tumor_gp_alpha;
  real<lower = 0> pop_tumor_gp_rho;
  real<lower = 0> pop_tumor_sigma;

  vector[n_all_tumor_measures] eta;  
}

transformed parameters {
  vector[n_all_tumor_measures] tumor_pred; 
  
  {
    int tumor_pos = 1;
    int t_measure_pos = 1;
    
    for (i in 1:n_patients) {
      int n_current_tumors = n_patient_tumors[i];
      int n_current_measures = n_measures[i];
      int t_measure_end = t_measure_pos + n_current_measures - 1;
      
      for (j in 1:n_current_tumors) {
        int tumor_end = tumor_pos + max_t[i] - 1; 
        
        tumor_pred[tumor_pos:tumor_end] = calc_gp_pred(
          all_measure_idx[:max_t[i]], pop_tumor_gp_intercept, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta, eta[tumor_pos:tumor_end]); 
        
        tumor_pos = tumor_end + 1;
      }
      
      t_measure_pos = t_measure_end + 1;
    }
  }
}

model {
  pop_tumor_gp_intercept ~ normal(0, 0.25);
  pop_tumor_gp_alpha ~ normal(0, 0.25);
  pop_tumor_gp_rho ~ inv_gamma(5, 5);
  pop_tumor_sigma ~ normal(0, 0.5);
  
  eta ~ std_normal();

  if (fit_data) {
    // to_vector(scaled_tumor_size) ~ lognormal(to_vector(tumor_pred), pop_tumor_sigma);
    
    scaled_tumor_size ~ lognormal(tumor_pred[obs_tumor_measures_idx], pop_tumor_sigma);
   
    // int tumor_pos = 1;
    // int t_measure_pos = 1;
    // 
    // for (i in 1:n_patients) {
    //   int n_current_tumors = n_patient_tumors[i];
    //   int n_current_measures = n_measures[i];
    //   int t_measure_end = t_measure_pos + n_current_measures - 1;
    //   int current_max_t = t_measure[t_measure_end];
    //   
    //   for (j in 1:n_current_tumors) {
    //     int tumor_end = tumor_pos + current_max_t - 1; 
    //     
    //     scaled_tumor_size[tumor_pos:tumor_end][t_measure[t_measure_pos:t_measure_end]] ~ 
    //       lognormal(tumor_pred[tumor_pos:tumor_end][t_measure[t_measure_pos:t_measure_end]], pop_tumor_sigma); 
    //     
    //     tumor_pos = tumor_end + 1;
    //   }
    //   
    //   t_measure_pos = t_measure_end + 1;
    // }
  }  
}

generated quantities {
  array[gen_tumor_sizes ? n_all_tumor_measures : 0] real<lower = 0> rep_tumor_size;
  
  if (gen_tumor_sizes) {
    rep_tumor_size = lognormal_rng(tumor_pred, pop_tumor_sigma);
  }
}

