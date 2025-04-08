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
  vector[last_predict_visit] pop_pred_tumor_gp;
  matrix[n_trials, last_predict_visit] trial_pred_tumor_gp; 
  matrix[n_patients, last_predict_visit] patient_pred_tumor_gp;
  
  // vector[n_trials] trial_K_obs_condition_number = zeros_vector(n_trials);
  // vector[n_trials] trial_K_obs_min_eigenvalue = zeros_vector(n_trials);
  // vector[n_trials] trial_K_obs_max_eigenvalue = zeros_vector(n_trials);
  // 
  // vector[n_patients] patient_K_obs_condition_number = zeros_vector(n_patients);
  // vector[n_patients] patient_K_obs_min_eigenvalue = zeros_vector(n_patients);
  // vector[n_patients] patient_K_obs_max_eigenvalue = zeros_vector(n_patients);
  // 
  // vector[n_patients] patient_K_cond_condition_number = zeros_vector(n_patients);
  // vector[n_patients] patient_K_cond_min_eigenvalue = zeros_vector(n_patients);
  // vector[n_patients] patient_K_cond_max_eigenvalue = zeros_vector(n_patients);
   
  if (predict_missing_sizes) {
    if (separate_trial_tumor_gp) {
      pop_pred_tumor_gp = zeros_vector(last_predict_visit);
    } else if (patient_gp_only) {
      pop_pred_tumor_gp = rep_vector(tumor_mean[1], last_predict_visit);
    } else {
      // matrix[last_predict_visit, last_predict_visit] pop_K = gp_matern52_cov(
      matrix[last_predict_visit, last_predict_visit] pop_K = gp_exp_quad_cov(
        all_tumor_measure_t[:last_predict_visit], pop_tumor_gp_alpha[1], pop_tumor_gp_rho[1], delta
      );
      
      pop_pred_tumor_gp[pop_unique_visits] = pop_obs_tumor_gp;
      
      // matrix[n_pop_unique_missing_visits, n_pop_unique_visits] pop_K_pred_obs = gp_matern52_cov(
      matrix[n_pop_unique_missing_visits, n_pop_unique_visits] pop_K_pred_obs = gp_exp_quad_cov(
        all_tumor_measure_t[pop_unique_missing_visits], all_tumor_measure_t[pop_unique_visits], pop_tumor_gp_alpha[1], pop_tumor_gp_rho[1]
      );
      
      pop_pred_tumor_gp[pop_unique_missing_visits] = multi_normal_rng(
        tumor_mean[1], pop_obs_tumor_gp, pop_K[pop_unique_visits, pop_unique_visits], pop_K_pred_obs, pop_K[pop_unique_missing_visits, pop_unique_missing_visits]
      );
    }

    for (s in 1:n_trials) {
      int actual_s = min(s, n_tumor_separate_trials);
      
      array[n_trial_unique_visits[s]] int curr_trial_visits = get_int_sub_array(trial_unique_visits, trial_unique_visits_pos, s); 
      array[n_trial_unique_missing_visits[s]] int curr_trial_missing_visits = get_int_sub_array(trial_unique_missing_visits, trial_unique_missing_visits_pos, s); 
      
     trial_pred_tumor_gp[s, curr_trial_visits] = get_sub_row_vector(trial_obs_tumor_gp, trial_unique_visits_pos, s);
      
     if (separate_trial_tumor_gp || add_trial_level_tumor) { 
       matrix [last_predict_visit, last_predict_visit] trial_K;
      
       if (separate_trial_tumor_gp) {
          // trial_K = gp_matern52_cov(
          trial_K = gp_exp_quad_cov(
            all_tumor_measure_t[:last_predict_visit], pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], delta
          );
          
          // matrix[n_trial_unique_missing_visits[s], n_trial_unique_visits[s]] trial_K_pred_obs = gp_matern52_cov(
          matrix[n_trial_unique_missing_visits[s], n_trial_unique_visits[s]] trial_K_pred_obs = gp_exp_quad_cov(
            all_tumor_measure_t[curr_trial_missing_visits], all_tumor_measure_t[curr_trial_visits], pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s]
          );
          
          trial_pred_tumor_gp[s, curr_trial_missing_visits] = multi_normal_rng(
            tumor_mean[s], 
            trial_pred_tumor_gp[s, curr_trial_visits]',
            trial_K[curr_trial_visits, curr_trial_visits], trial_K_pred_obs, trial_K[curr_trial_missing_visits, curr_trial_missing_visits]
          )';
        } else if (patient_gp_only) {
          trial_pred_tumor_gp[s, curr_trial_missing_visits] = pop_pred_tumor_gp[curr_trial_missing_visits]' + trial_tumor_gp_intercept_effect[s];
        } else {
          // trial_K = gp_matern52_cov(
          trial_K = gp_exp_quad_cov(
            all_tumor_measure_t[:last_predict_visit], trial_tumor_gp_alpha, trial_tumor_gp_rho, delta
          );
          
          // matrix[n_trial_unique_missing_visits[s], n_trial_unique_visits[s]] trial_K_pred_obs = gp_matern52_cov(
          matrix[n_trial_unique_missing_visits[s], n_trial_unique_visits[s]] trial_K_pred_obs = gp_exp_quad_cov(
            all_tumor_measure_t[curr_trial_missing_visits], all_tumor_measure_t[curr_trial_visits], trial_tumor_gp_alpha, trial_tumor_gp_rho
          );
          
          trial_pred_tumor_gp[s, curr_trial_missing_visits] = multi_normal_rng(
            pop_pred_tumor_gp[curr_trial_visits] + trial_tumor_gp_intercept_effect[s],
            pop_pred_tumor_gp[curr_trial_missing_visits] + trial_tumor_gp_intercept_effect[s],
            trial_pred_tumor_gp[s, curr_trial_visits]',
            trial_K[curr_trial_visits, curr_trial_visits], trial_K_pred_obs, trial_K[curr_trial_missing_visits, curr_trial_missing_visits]
          )';
        }
        
        // (trial_K_obs_min_eigenvalue[s], trial_K_obs_max_eigenvalue[s], trial_K_obs_condition_number[s]) = summarize_matrix_eigenvalues(trial_K[curr_trial_visits, curr_trial_visits]);
     }
      
      int curr_trial_patient_pos, curr_trial_patient_end;
      (curr_trial_patient_pos, curr_trial_patient_end) = get_pos(trial_patient_pos, s);

      for (i in curr_trial_patient_pos:curr_trial_patient_end) {
        // matrix [last_predict_visit, last_predict_visit] patient_K = gp_matern52_cov(
        matrix [last_predict_visit, last_predict_visit] patient_K = gp_exp_quad_cov(
          // all_tumor_measure_t[:last_predict_visit], patient_tumor_gp_alpha[actual_s], patient_rho[i], tumor_sd[actual_s]^2 + delta
          all_tumor_measure_t[:last_predict_visit], patient_tumor_gp_alpha[actual_s], exp(log_patient_rho[i]), delta
        );
        
        int n_curr_measured_visits = last_predict_visit - n_patient_unique_mn_visits[i];
        int n_curr_mn_visits = n_patient_unique_mn_visits[i];
        
        array[n_curr_measured_visits] int curr_patient_measured_visits = get_int_sub_array(measured_tumor_visits, patient_measured_tumor_visits_pos, i);
        array[n_curr_mn_visits] int curr_patient_mn_visits = get_int_sub_array(patient_unique_mn_visits, patient_unique_mn_visits_pos, i); 
        
        row_vector[n_patient_unique_visits[i]] curr_patient_sld = get_sub_row_vector(post_treat_sld, post_treat_visits_pos, i);
        row_vector[n_patient_unique_visits[i]] curr_patient_gp = get_sub_row_vector(patient_obs_tumor_gp, patient_unique_visits_pos, i);
        array[n_curr_measured_visits] int measured2patient_idx = get_int_sub_array(measured2patient_visits_idx, patient_measured_tumor_visits_pos, i); 
        
        if (n_curr_measured_visits > 0) { 
          patient_pred_tumor_gp[i, curr_patient_measured_visits] = curr_patient_gp[measured2patient_idx];
          
          // matrix[n_curr_mn_visits, n_curr_measured_visits] patient_K_pred_obs = gp_matern52_cov(
          matrix[n_curr_mn_visits, n_curr_measured_visits] patient_K_pred_obs = gp_exp_quad_cov(
            all_tumor_measure_t[curr_patient_mn_visits], all_tumor_measure_t[curr_patient_measured_visits], patient_tumor_gp_alpha[actual_s], exp(log_patient_rho[i])
          );
          
          // (patient_K_obs_min_eigenvalue[i], patient_K_obs_max_eigenvalue[i], patient_K_obs_condition_number[i]) = summarize_matrix_eigenvalues(
          //   patient_K[curr_patient_measured_visits, curr_patient_measured_visits]
          // );
          
          // vector[n_curr_mn_visits] mu_cond;
          // matrix[n_curr_mn_visits, n_curr_mn_visits] Sigma_cond;
          // 
          // (mu_cond, Sigma_cond) = gp_conditional(
          //   pop_pred_tumor_gp[curr_patient_measured_visits] + trial_pred_tumor_gp[s, curr_patient_measured_visits]' + patient_tumor_gp_intercept_effect[i],
          //   pop_pred_tumor_gp[curr_patient_mn_visits] + trial_pred_tumor_gp[s, curr_patient_mn_visits]' + patient_tumor_gp_intercept_effect[i],
          //   measured_sld',
          //   patient_K[curr_patient_measured_visits, curr_patient_measured_visits], patient_K_pred_obs, patient_K[curr_patient_mn_visits, curr_patient_mn_visits],
          //   delta
          // );
          // 
          // Sigma_cond = (Sigma_cond + Sigma_cond') / 2;
          // 
          // (patient_K_cond_min_eigenvalue[i], patient_K_cond_max_eigenvalue[i], patient_K_cond_condition_number[i]) = summarize_matrix_eigenvalues(Sigma_cond);
          
          patient_pred_tumor_gp[i, curr_patient_mn_visits] = multi_normal_rng(
            trial_pred_tumor_gp[s, curr_patient_measured_visits]' + patient_tumor_gp_intercept_effect[i],
            trial_pred_tumor_gp[s, curr_patient_mn_visits]' + patient_tumor_gp_intercept_effect[i],
            patient_pred_tumor_gp[i, curr_patient_measured_visits]',   // measured_sld',
            patient_K[curr_patient_measured_visits, curr_patient_measured_visits], patient_K_pred_obs, patient_K[curr_patient_mn_visits, curr_patient_mn_visits]
          )';
        } else {
          patient_pred_tumor_gp[i, curr_patient_mn_visits] = multi_normal_rng(
            // pop_pred_tumor_gp[curr_patient_mn_visits] + trial_pred_tumor_gp[s, curr_patient_mn_visits]' + patient_tumor_gp_intercept_effect[i],
            trial_pred_tumor_gp[s, curr_patient_mn_visits]' + patient_tumor_gp_intercept_effect[i],
            patient_K[curr_patient_mn_visits, curr_patient_mn_visits]
          )';
        }
      }
    }
  }
  
  vector<lower = 0, upper = 1>[separate_trial_tumor_gp || patient_gp_only ? 0 : max_t_width] pop_rho_corr; 
  array[patient_gp_only ? 0 : n_tumor_separate_trials] vector<lower = 0, upper = 1>[max_t_width] trial_rho_corr; 
  array[n_patients] vector<lower = 0, upper = 1>[max_t_width] patient_rho_corr; 
  
  if (!separate_trial_tumor_gp) {
    if (!patient_gp_only) {
      pop_rho_corr = exp(-0.5 * square(to_vector(all_tumor_measure_t) / pop_tumor_gp_rho[1]));
      
      if (add_trial_level_tumor) {
        trial_rho_corr[1] = exp(-0.5 * square(to_vector(all_tumor_measure_t) / trial_tumor_gp_rho));
      }
    }
  } else {
    for (s in 1:n_trials) {
      trial_rho_corr[s] = exp(-0.5 * square(to_vector(all_tumor_measure_t) / pop_tumor_gp_rho[s]));
    }
  }
 
  for (i in 1:n_patients) { 
    patient_rho_corr[i] = exp(-0.5 * square(to_vector(all_tumor_measure_t) / exp(log_patient_rho[i])));
  }
}
