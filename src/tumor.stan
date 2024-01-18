functions {
  #include "util.stan"
}

data {
  #include "base_data.stan"
  #include "tumor_data.stan"
}

transformed data {
  int use_tumor_model = 1;
  
  #include "tumor_transformed_data.stan"
}

parameters {
  #include "tumor_parameters.stan"
}

transformed parameters {
  #include "tumor_transformed_parameters.stan"
}

model {
  #include "tumor_model.stan"
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
