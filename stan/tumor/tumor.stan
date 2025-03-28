functions {
  #include "../extern_util.stan"
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../state_space.stan"
}

data {
  #include "../base_data.stan"
  #include "tumor_data.stan"
}

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
}

parameters {
  #include "tumor_parameters.stan"
}

transformed parameters {
  #include "tumor_transformed_parameters.stan"
}

model {
  #include "tumor_priors.stan"
  #include "tumor_model.stan"
}

generated quantities {
  vector[predict_missing_sizes ? sum(n_patient_full_visits) : 0] all_latent_tumor_gp;
  vector<lower = 0>[predict_missing_sizes ? sum(n_patient_full_visits) : 0] all_sum_tumor_size; // This includes both observed and missing measurements

  if (predict_missing_sizes) {
    for (s in 1:n_trials) {
      int curr_trial_patient_pos, curr_trial_patient_end;
      (curr_trial_patient_pos, curr_trial_patient_end) = get_pos(trial_patient_pos, s);
      
      int actual_s = min(s, n_tumor_separate_trials);
      
      matrix[max_t_width, max_t_width] K_pop = calc_gp_vcov(all_measure_t, pop_tumor_gp_alpha[actual_s], pop_tumor_gp_rho[actual_s], delta);
  
      for (i in curr_trial_patient_pos:curr_trial_patient_end) {
        // print("i = ", i);
        
        int full_visit_pos, full_visit_end;
        (full_visit_pos, full_visit_end) = get_pos(patient_full_visits_pos, i);
        
        array[n_patient_visits[i] - n_patient_non_measured_tumors[i]] int curr_measured_tumors = get_int_sub_array(global_measured_tumors_idx, patient_measured_tumors_pos, i);
        array[n_patient_non_measured_tumors[i]] int curr_non_measured_tumors = get_int_sub_array(global_non_measured_tumors_idx, patient_non_measured_tumors_pos, i);
        
        array[n_patient_visits[i]] int curr_visit_idx = get_int_sub_array(t_patient_visit_idx, patient_visit_pos, i); 
        array[n_patient_visits[i]] int curr_visit_full_idx = get_int_sub_array(patient_visit_full_idx, patient_visit_pos, i); 
        all_latent_tumor_gp[curr_visit_full_idx] = latent_pop_tumor_gp[actual_s, pop_unique_visits_idx_dict[sort_asc(append_array(curr_measured_tumors, curr_non_measured_tumors))]];
        all_sum_tumor_size[curr_visit_full_idx] = get_sub_vector(sum_tumor_size, patient_visit_pos, i);
        
        if (n_patient_missing_visits[i] > 0) {
          array[n_patient_missing_visits[i]] int curr_missing_visit_idx = get_int_sub_array(t_patient_missing_visits_idx, patient_missing_visits_pos, i);
          array[n_patient_missing_visits[i]] int curr_missing_visit_full_idx = get_int_sub_array(patient_missing_visit_full_idx, patient_missing_visits_pos, i);
        
          all_latent_tumor_gp[curr_missing_visit_full_idx] = gp_pred_rng(
            all_measure_t[curr_missing_visit_idx],
            all_latent_tumor_gp[curr_visit_full_idx],
            all_measure_t[curr_visit_idx],
            K_pop[curr_visit_idx, curr_visit_idx],
            K_pop[curr_missing_visit_idx, curr_missing_visit_idx],
            pop_tumor_gp_alpha[actual_s],
            pop_tumor_gp_rho[actual_s],
            delta
          );
          
          all_sum_tumor_size[curr_missing_visit_full_idx] = fmax( 
            to_vector(normal_rng(tumor_mean[actual_s] + all_latent_tumor_gp[curr_missing_visit_full_idx], tumor_sd[actual_s])),
            zeros_vector(n_patient_missing_visits[i])
          );
        }
      }
    }
  }
  
  /*
  // vector<lower = 0>[gen_tumor_sizes ? sum(n_measures) : 0] rep_tumor_size; // Drawing *new* data
  // 
  if (predict_missing_sizes) { // || gen_tumor_sizes) {
  //   matrix[max(patient_max_t_width), max(patient_max_t_width)] pop_tumor_vcov = calc_gp_vcov(all_measure_t, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
  //   
  //   int tumor_pos = 1;
    int full_visit_pos = 1;
  //   int obs_measure_pos = 1;
  //   int missing_measure_pos = 1;
  //   
  
    for (s in 1:n_trials) {
      int curr_patient_pos = trial_patient_pos[s];
      int curr_patient_end = trial_patient_pos[s + 1] - 1;
      
      for (i in curr_patient_pos:curr_patient_end) {
        int full_visit_end = full_visit_pos + n_patient_full_visits[i] - 1;
        int curr_patient_visits_pos = patient_visit_pos[i];
        int curr_patient_visits_end = patient_visit_pos[i + 1] - 1;
        int curr_patient_missing_visits_pos = patient_missing_visit_pos[i];
        int curr_patient_missing_visits_end = patient_missing_visit_pos[i + 1] - 1;
        
        int n_visits = n_patient_visits[i], n_missing_visits = n_patient_missing_visits[i];
        array[n_visits] int visit_idx = patient_t_visit_idx[curr_patient_visits_pos:curr_patient_visits_end];  
        array[n_missing_visits] int missing_visit_idx = patient_t_missing_visits_idx[curr_patient_missing_visits_pos:curr_patient_missing_visits_end];  
        
        vector[n_patient_full_visits[i]] current_log_sum_tumor_size;
        current_log_sum_tumor_size[visit_idx] = log_sum_tumor_size[curr_patient_visits_pos:curr_patient_visits_end];
        
        matrix[n_visits, n_visits] patient_tumor_gp_cov = trial_tumor_gp_cov[s, visit_idx, visit_idx];
        matrix[n_missing_visits, n_missing_visits] patient_missing_tumor_gp_cov = trial_tumor_gp_cov[s, missing_visit_idx, missing_visit_idx];
        
        if (n_missing_visits > 0) {
          current_log_sum_tumor_size[missing_visit_idx] = 
              gp_pred_rng( // This function I took from the Stan docs that makes it easy to make these draws.
                all_measure_t[missing_visit_idx],
                current_log_sum_tumor_size[visit_idx],
                all_measure_t[visit_idx],
                patient_tumor_gp_cov,
                patient_missing_tumor_gp_cov,
                pop_tumor_gp_alpha[s], 
                pop_tumor_gp_rho[s] 
                // pop_tumor_sigma[s]^2
              );
        }
        
        all_sum_tumor_size[full_visit_pos:full_visit_end] = exp(current_log_sum_tumor_size);
      
        full_visit_pos = full_visit_end + 1;
      }
    }
    // for (i in 1:n_patients) {
  //     for (j in 1:n_patient_tumors[i]) {
  //       int full_measure_end = full_measure_pos + n_full_measures[tumor_pos] - 1;
  //       int obs_measure_end = obs_measure_pos + n_measures[tumor_pos] - 1;
  //       int missing_measure_end = missing_measure_pos + n_missing_measures[tumor_pos] - 1;
  //       
  //       array[n_measures[tumor_pos]] int current_t_measure_idx = patient_t_measure_idx[obs_measure_pos:obs_measure_end];
  //       
  //       if (predict_missing_sizes) {
  //         array[n_missing_measures[tumor_pos]] int current_t_missing_measure_idx = patient_t_missing_measure_idx[missing_measure_pos:missing_measure_end];
  //         
  //         // vector[patient_max_t_width[i]] current_tumor_size;
  //         vector[n_full_measures[tumor_pos]] current_tumor_size;
  //         
  //         // print("i = ", i, ", j = ", j, " current_t_measure_idx = ", current_t_measure_idx);
  //         
  //         current_tumor_size[current_t_measure_idx] = tumor_size[obs_measure_pos:obs_measure_end];
  //         
  //         if (n_missing_measures[tumor_pos] > 0) {
  //           // Here we are not simply drawing tumor sizes from the GP distribution conditional on observed data. We need to make sure that the new
  //           // imputations are jointly determined with the observed data, i.e., are correlated to the observed data.
  //           
  //           current_tumor_size[current_t_missing_measure_idx] = 
  //             exp(gp_pred_rng( // This function I took from the Stan docs that makes it easy to make these draws.
  //               all_measure_t[current_t_missing_measure_idx],
  //               log(current_tumor_size[current_t_measure_idx]), 
  //               all_measure_t[current_t_measure_idx], 
  //               pop_tumor_vcov[current_t_measure_idx, current_t_measure_idx],
  //               pop_tumor_gp_alpha,
  //               pop_tumor_gp_rho,
  //               delta
  //             )); 
  //         }
  //         
  //         all_tumor_size[full_measure_pos:full_measure_end] = current_tumor_size;
  //       }
  //       
  //       if (gen_tumor_sizes) {
  //         rep_tumor_size[obs_measure_pos:obs_measure_end] = exp(multi_normal_rng(
  //           rep_vector(tumor_gp_intercept[tumor_pos], n_measures[tumor_pos]), 
  //           pop_tumor_vcov[current_t_measure_idx, current_t_measure_idx]
  //         ));
  //       }
  //       
  //       full_measure_pos = full_measure_end + 1;
  //       obs_measure_pos = obs_measure_end + 1;
  //       missing_measure_pos = missing_measure_end + 1;
  //       tumor_pos += 1;
  //     }
    // }
  }
  */
}
