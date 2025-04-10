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
  array[n_patients] vector[last_predict_visit] patient_pred_latent_sld;
   
  if (predict_missing_sizes) {
    for (i in 1:n_patients) {
      // matrix [last_predict_visit, last_predict_visit] patient_K = gp_matern52_cov(
      matrix [last_predict_visit, last_predict_visit] patient_K = gp_exp_quad_cov(
        all_tumor_measure_t[:last_predict_visit], pop_tumor_gp_alpha, patient_tumor_gp_rho[i], delta
      );
      
      int n_curr_measured_visits = last_predict_visit - n_patient_unique_mn_visits[i];
      int n_curr_mn_visits = n_patient_unique_mn_visits[i];
      
      array[n_curr_measured_visits] int curr_patient_measured_visits = get_int_sub_array(measured_tumor_visits, patient_measured_tumor_visits_pos, i);
      array[n_curr_mn_visits] int curr_patient_mn_visits = get_int_sub_array(patient_unique_mn_visits, patient_unique_mn_visits_pos, i); 
      
      vector[n_patient_unique_visits[i]] curr_patient_sld = get_sub_vector(post_treat_sld, post_treat_visits_pos, i);
      vector[n_patient_unique_visits[i]] curr_patient_gp = get_sub_vector(patient_obs_tumor_gp, patient_unique_visits_pos, i);
      array[n_curr_measured_visits] int measured2patient_idx = get_int_sub_array(measured2patient_visits_idx, patient_measured_tumor_visits_pos, i); 
      
      if (n_curr_measured_visits > 0) { 
        patient_pred_latent_sld[i, curr_patient_measured_visits] = curr_patient_gp[measured2patient_idx];
        
        // matrix[n_curr_mn_visits, n_curr_measured_visits] patient_K_pred_obs = gp_matern52_cov(
        matrix[n_curr_mn_visits, n_curr_measured_visits] patient_K_pred_obs = gp_exp_quad_cov(
          all_tumor_measure_t[curr_patient_mn_visits], all_tumor_measure_t[curr_patient_measured_visits], pop_tumor_gp_alpha, patient_tumor_gp_rho[i]
        );
        
        patient_pred_latent_sld[i, curr_patient_mn_visits] = multi_normal_rng(
          0, 0,
          patient_pred_latent_sld[i, curr_patient_measured_visits],
          patient_K[curr_patient_measured_visits, curr_patient_measured_visits], patient_K_pred_obs, patient_K[curr_patient_mn_visits, curr_patient_mn_visits]
        );
      } else {
        patient_pred_latent_sld[i, curr_patient_mn_visits] = multi_normal_rng(
          zeros_vector(n_curr_mn_visits),
          patient_K[curr_patient_mn_visits, curr_patient_mn_visits]
        );
      }
    }
  }
  
  vector<lower = 0, upper = 1>[max_t_width] pop_rho_corr = gp_exp_quad_cov(all_tumor_measure_t, 1, exp(log_pop_tumor_gp_rho))[, 1]; 
  array[n_trials] vector<lower = 0, upper = 1>[max_t_width] trial_rho_corr; 
  array[n_patients] vector<lower = 0, upper = 1>[max_t_width] patient_rho_corr; 
  
  for (s in 1:n_trials) {
    trial_rho_corr[s] = gp_exp_quad_cov(all_tumor_measure_t, 1, exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[s]))[, 1]; 
  }
 
  for (i in 1:n_patients) { 
    patient_rho_corr[i] = gp_exp_quad_cov(all_tumor_measure_t, 1, patient_tumor_gp_rho[i])[, 1]; 
  }
}
