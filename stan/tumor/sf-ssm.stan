functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "sf-ssm_functions.stan"
}  

data {
  #include "../base_data.stan"
  
  int<lower = 0, upper = 1> fit_tumor_data;
  int<lower = 0> sf_rep_T;
  
  real<lower = 0> decrease_process_sd_sd;
  real<lower = 0> growth_process_sd_sd;
  real<lower = 0> decrease_process_alpha;
  real<lower = 0> decrease_process_beta;
  real<lower = 0> growth_process_alpha;
  real<lower = 0> growth_process_beta;
  real<lower = 0> measure_sd_sd;
  
  real log_decrease_rate_mean;
  real<lower = 0> log_decrease_rate_sd;
  real log_growth_rate_mean;
  real<lower = 0> log_growth_rate_sd;
  real<lower = 0> patient_log_decrease_rate_sd_sd;
  real<lower = 0> patient_log_growth_rate_sd_sd;
  
  int test_patients_pos, test_patients_end;
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
 
  row_vector[2] H = [1, 1];
}

parameters {
  real log_decrease_rate;      // Tumor regression rate (d)
  real log_growth_rate;        // Tumor growth rate (g)
  real<lower=0> decrease_process_sd;
  real<lower=0> growth_process_sd;
  real<lower=0> measure_sd;        // Measurement noise variance
  
  vector[n_patients] raw_patient_log_decrease_rate_effect;      
  real<lower = 0> patient_log_decrease_rate_sd;
  vector[n_patients] raw_patient_log_growth_rate_effect;      
  real<lower = 0> patient_log_growth_rate_sd;
}

transformed parameters {
  vector[n_patients] patient_log_decrease_rate_effect = raw_patient_log_decrease_rate_effect * patient_log_decrease_rate_sd;      
  vector[n_patients] patient_log_growth_rate_effect = raw_patient_log_growth_rate_effect * patient_log_growth_rate_sd;      
  
  vector<lower = 0>[n_patients] patient_decrease_rate = exp(log_decrease_rate + patient_log_decrease_rate_effect); 
  vector<lower = 0>[n_patients] patient_growth_rate = exp(log_growth_rate + patient_log_growth_rate_effect); 
  
  print("exp(log_decrease_rate) = ", exp(log_decrease_rate), ", exp(log_growth_rate) = ", exp(log_growth_rate));
  
  // for (i in linspaced_int_array(n_test_patients, 1, n_patients)) {
  for (i in test_patients_pos:test_patients_end) {
    vector[2] F_diag = diagonal(sf_create_F(patient_decrease_rate[i], patient_growth_rate[i], 1));
     
    print(
      i, ": F_diag = ", F_diag, ", condition number = ", max(F_diag) / min(F_diag), 
      ", patient_decrease_rate[i] = ", patient_decrease_rate[i], 
      ", patient_growth_rate[i] = ", patient_growth_rate[i] 
    );
  }
}

model {
  decrease_process_sd ~ normal(0, decrease_process_sd_sd);
  growth_process_sd ~ normal(0, growth_process_sd_sd);
  // decrease_process_sd ~ inv_gamma(decrease_process_alpha, decrease_process_beta);
  // growth_process_sd ~ inv_gamma(growth_process_alpha, growth_process_beta);
  measure_sd ~ normal(0, measure_sd_sd);
  
  log_decrease_rate ~ normal(log_decrease_rate_mean, log_decrease_rate_sd);
  log_growth_rate ~ normal(log_growth_rate_mean, log_growth_rate_sd);
  
  raw_patient_log_decrease_rate_effect ~ std_normal();
  raw_patient_log_growth_rate_effect ~ std_normal();
  
  patient_log_decrease_rate_sd ~ normal(0, patient_log_decrease_rate_sd_sd);
  patient_log_growth_rate_sd ~ normal(0, patient_log_growth_rate_sd_sd);
  
 
  if (fit_tumor_data) { 
    for (i in 1:n_patients) {
      int visit_pos, visit_end;
      (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
      
      tuple(matrix[2, n_patient_visits[i]], array[n_patient_visits[i]] matrix[2, 2]) forward_res = 
        // sf_forward_filter_sqrt(
        sf_forward_filter(
          normalized_sld[visit_pos:visit_end], t_patient_visits[visit_pos:visit_end], 
          patient_decrease_rate[i], patient_growth_rate[i], 
          decrease_process_sd, growth_process_sd, measure_sd
        );
        
      matrix[2, n_patient_visits[i]] filtered_means = forward_res.1;
        
      row_vector[n_patient_visits[i]] pred_sld = H * filtered_means - 1.0;
     
      // The first observation is ignored. The normalization might be messing things up. 
      for (v in 2:n_patient_visits[i]) {
        int v_idx = visit_pos + v - 1; 
        
        if (normalized_sld[v_idx] > 0) {
          normalized_sld[v_idx] ~ normal(pred_sld[v], measure_sd);
        } else {
          target += normal_lcdf(lod | pred_sld[v], measure_sd);
        }
      }
    }
  }
}

generated quantities {
  array[n_patients] matrix[2, sf_rep_T + 1] rep_state = rep_array(rep_matrix(0, 2, sf_rep_T + 1), n_patients);
  array[n_patients] vector<lower = 0>[sf_rep_T + 1] rep_sld = rep_array(zeros_vector(sf_rep_T + 1), n_patients);

  if (sf_rep_T > 0) {
    // for (i in 1:n_patients) {
    // for (i in linspaced_int_array(10, 1, n_patients)) {
    for (i in test_patients_pos:test_patients_end) {
      int visit_pos, visit_end;
      (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
      
      // print(i, ": -");

      if (t_patient_visits[visit_pos] <= 0) {
        // (rep_state[i], rep_sld[i]) = sf_forecast_sqrt_rng(
        (rep_state[i], rep_sld[i]) = sf_forecast_rng(
          normalized_sld[visit_pos:visit_pos], t_patient_visits[visit_pos:visit_pos], rep_time_points, 
          patient_decrease_rate[i], patient_growth_rate[i], 
          decrease_process_sd, growth_process_sd, measure_sd
        );

        rep_sld[i] *= sum_tumor_size[visit_pos];

        for (t in 1:sf_rep_T) {
          if (rep_sld[i, t + 1] < lod) {
            rep_sld[i, t + 1] = 0.0;
          }
        }
      }
    }
  }
}