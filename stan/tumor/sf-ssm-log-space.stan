functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "sf-ssls_functions.stan"
}  

data {
  #include "../base_data.stan"
  
  int<lower = 0, upper = 1> fit_tumor_data;
  int<lower = 0> sf_rep_T;
  int<lower = 0, upper = 1> debug;
  
  real<lower = 0> decrease_process_sd_sd;
  real<lower = 0> growth_process_sd_sd;
  real<lower = 0> decrease_process_alpha;
  real<lower = 0> decrease_process_beta;
  real<lower = 0> growth_process_alpha;
  real<lower = 0> growth_process_beta;
  real<lower = 0> process_corr_param;
  real<lower = 0> measure_sd_sd;
  
  real log_decrease_rate_mean;
  real<lower = 0> log_decrease_rate_sd;
  real log_growth_rate_mean;
  real<lower = 0> log_growth_rate_sd;
  real log_rate_ratio_mean;
  real<lower = 0> log_rate_ratio_sd;
  real<lower = 0> growth_lag_alpha;
  real<lower = 0> growth_lag_beta;
  real growth_lag_mean;
  real<lower = 0> growth_lag_sd;
  real<lower = 0> patient_log_growth_lag_sd_sd;
  real<lower = 0> growth_transition_rate_sd;
  real<lower = 0> patient_log_decrease_rate_sd_sd;
  real<lower = 0> patient_log_growth_rate_sd_sd;
  real<lower = 0> rate_corr_param;
  real<lower = 0> pop_decrease_prop_logis_sd;
  real<lower = 0> patient_decrease_prop_logis_sd_sd;
  real<lower = 0> log_lod_sd;
  
  int test_patients_pos, test_patients_end;
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
  
  int<lower = 1> n_total_visits_m1 = sum(n_patient_visits) - n_patients;
  array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
  
  real log_lod = log(0.1);
}

parameters {
  real pop_log_decrease_rate;      // Tumor regression rate (d)
  real<lower = 0> log_rate_ratio; // d > g
  // real pop_log_growth_rate;        // Tumor growth rate (g)
  // real<lower=0> log_rate_difference; // Must be positive
  
  real pop_log_growth_lag;
  real<lower = 0> growth_transition_rate;
  
  real<lower = 0> patient_log_growth_lag_sd;
  vector<offset = pop_log_growth_lag, multiplier = patient_log_growth_lag_sd>[n_patients] patient_log_growth_lag;
  
  real<lower=0> decrease_process_sd;
  real<lower=0> growth_process_sd;
  cholesky_factor_corr[2] L_process_corr;
  real<lower=0> measure_sd;        // Measurement noise variance
  
  // vector[n_patients] raw_patient_log_decrease_rate_effect;
  real<lower = 0> patient_log_decrease_rate_sd;
  // vector<offset = pop_log_decrease_rate, multiplier = patient_log_decrease_rate_sd>[n_patients] patient_log_decrease_rate;
  // vector[n_patients] raw_patient_log_growth_rate_effect;
  real<lower = 0> patient_log_growth_rate_sd;
  // vector<offset = pop_log_growth_rate, multiplier = patient_log_growth_rate_sd>[n_patients] patient_log_growth_rate;
  // vector<multiplier = patient_log_growth_rate_sd>[n_patients] patient_log_growth_rate_nc;
  cholesky_factor_corr[2] L_rate_corr;
  matrix[n_patients, 2] raw_patient_log_rate_effect;
  
  real pop_decrease_prop_logis;

  // // vector[n_patients] raw_patient_decrease_prop_logis_effect;
  real<lower = 0> patient_decrease_prop_logis_sd;
  vector<offset = pop_decrease_prop_logis, multiplier = patient_decrease_prop_logis_sd>[n_patients] patient_decrease_prop_logis;
  
  // real log_lod;
  
  matrix[n_total_visits_m1, 2] raw_states; 
}

transformed parameters {
  real pop_log_growth_rate = pop_log_decrease_rate - log_rate_ratio;
  // vector[n_patients] patient_log_growth_rate = pop_log_growth_rate + patient_log_growth_rate_nc;
 
  cholesky_factor_cov[2] L_rate_cov = diag_pre_multiply([patient_log_decrease_rate_sd, patient_log_growth_rate_sd]', L_rate_corr); 
  matrix[n_patients, 2] patient_log_rate = rep_matrix([pop_log_decrease_rate, pop_log_growth_rate], n_patients) + raw_patient_log_rate_effect * L_rate_cov';
  
  // vector[n_patients] patient_decrease_prop_logis = pop_decrease_prop_logis + patient_decrease_prop_logis_effect;      
  // real pop_log_decrease_prop = -log1p_exp(-pop_decrease_prop_logis);
  // real pop_log_growth_prop = pop_log_decrease_prop - pop_decrease_prop_logis;
  vector[n_patients] patient_log_decrease_prop = -log1p_exp(-patient_decrease_prop_logis);
  vector[n_patients] patient_log_growth_prop = patient_log_decrease_prop - patient_decrease_prop_logis;
  
  // real pop_log_growth_rate = pop_log_decrease_rate - log_rate_difference;
  
  // matrix[fit_tumor_data ? sum(n_patient_visits[1:1]) : 0, 2] states; 
  matrix[sum(n_patient_visits[test_patients_pos:test_patients_end]), 2] expected_states; 
  matrix[sum(n_patient_visits[test_patients_pos:test_patients_end]), 2] states; 
  
  for (i in test_patients_pos:test_patients_end) {
    int visit_pos, visit_end;
    (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
    int visit_m1_pos, visit_m1_end;
    (visit_m1_pos, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
    
    row_vector[2] initial_state = [patient_log_decrease_prop[i], patient_log_growth_prop[i]]; //  log([0.99, 0.01]);
    
    (expected_states[visit_pos:visit_end], states[visit_pos:visit_end]) = sf_log_space_trajectory_ncp(
      raw_states[visit_m1_pos:visit_m1_end], initial_state, t_patient_visits[visit_pos:visit_end],
      // exp(patient_log_decrease_rate[i]), exp(patient_log_growth_rate[i]), exp(patient_log_growth_lag[i]), growth_transition_rate,
      exp(patient_log_rate[i, 1]), exp(patient_log_rate[i, 2]), exp(patient_log_growth_lag[i]), growth_transition_rate,
      [decrease_process_sd, growth_process_sd]', L_process_corr
    );
  }
}

model {
  decrease_process_sd ~ normal(0, decrease_process_sd_sd);
  growth_process_sd ~ normal(0, growth_process_sd_sd);
  L_process_corr ~ lkj_corr_cholesky(process_corr_param);
  
  measure_sd ~ normal(0, measure_sd_sd);
  
  pop_log_decrease_rate ~ normal(log_decrease_rate_mean, log_decrease_rate_sd);
  log_rate_ratio ~ normal(log_rate_ratio_mean, log_rate_ratio_sd); 
  // pop_log_growth_rate ~ normal(log_growth_rate_mean, log_growth_rate_sd);
  // log_rate_difference ~ normal(2, 0.5);
  // growth_lag ~ inv_gamma(growth_lag_alpha, growth_lag_beta);
  pop_log_growth_lag ~ normal(growth_lag_mean, growth_lag_sd);
  growth_transition_rate ~ normal(0, growth_transition_rate_sd);
 
  patient_log_growth_lag_sd ~ normal(0, patient_log_growth_lag_sd_sd); 
  patient_log_growth_lag ~ normal(pop_log_growth_lag, patient_log_growth_lag_sd);
  
  patient_log_decrease_rate_sd ~ normal(0, patient_log_decrease_rate_sd_sd);
  // patient_log_decrease_rate ~ normal(pop_log_decrease_rate, patient_log_decrease_rate_sd);
  patient_log_growth_rate_sd ~ normal(0, patient_log_growth_rate_sd_sd);
  // patient_log_growth_rate_nc ~ normal(0, patient_log_growth_rate_sd);
  L_rate_corr ~ lkj_corr_cholesky(rate_corr_param);
  to_vector(raw_patient_log_rate_effect) ~ std_normal();
  
  pop_decrease_prop_logis ~ normal(0, pop_decrease_prop_logis_sd);
  // // raw_patient_decrease_prop_logis_effect ~ std_normal();
  patient_decrease_prop_logis_sd ~ normal(0, patient_decrease_prop_logis_sd_sd);
  patient_decrease_prop_logis ~ normal(pop_decrease_prop_logis, patient_decrease_prop_logis_sd);
  
  // log_lod ~ normal(log(lod), log_lod_sd);
  
  to_vector(raw_states) ~ std_normal();
 
  if (debug) { 
    print("exp(pop_log_decrease_rate) = ", exp(pop_log_decrease_rate), ", exp(pop_log_growth_rate) = ", exp(pop_log_growth_rate));
  }
 
  if (fit_tumor_data) { 
    // for (i in 1:n_patients) {
    for (i in test_patients_pos:test_patients_end) {
      int visit_pos, visit_end;
      (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
      
      normalized_sld[visit_pos:visit_end] ~ sf_log_space_obs(states[visit_pos:visit_end], measure_sd, log_lod - log(sum_tumor_size[1]));
  
      if (debug) {    
        print(i, ": normalized_sld = ", normalized_sld[visit_pos:visit_end], ", exp(states) = ", exp(states));
      }
    }
  }
}

generated quantities {
  corr_matrix[2] process_corr = L_process_corr * L_process_corr';
  corr_matrix[2] rate_corr = L_rate_corr * L_rate_corr';
  
  // array[n_patients] matrix[2, sf_rep_T + 1] rep_state = rep_array(rep_matrix(0, 2, sf_rep_T + 1), n_patients);
  // array[n_patients] vector<lower = 0>[sf_rep_T + 1] rep_sld = rep_array(zeros_vector(sf_rep_T + 1), n_patients);
  // 
  // if (sf_rep_T > 0) {
  //   // for (i in 1:n_patients) {
  //   // for (i in linspaced_int_array(10, 1, n_patients)) {
  //   for (i in test_patients_pos:test_patients_end) {
  //     int visit_pos, visit_end;
  //     (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
  //     
  //     // print(i, ": -");
  // 
  //     if (t_patient_visits[visit_pos] <= 0) {
  //       // (rep_state[i], rep_sld[i]) = sf_forecast_sqrt_rng(
  //       (rep_state[i], rep_sld[i]) = sf_forecast_rng(
  //         normalized_sld[visit_pos:visit_pos], t_patient_visits[visit_pos:visit_pos], rep_time_points, 
  //         patient_decrease_rate[i], patient_growth_rate[i], 
  //         decrease_process_sd, growth_process_sd, measure_sd
  //       );
  // 
  //       rep_sld[i] *= sum_tumor_size[visit_pos];
  // 
  //       for (t in 1:sf_rep_T) {
  //         if (rep_sld[i, t + 1] < lod) {
  //           rep_sld[i, t + 1] = 0.0;
  //         }
  //       }
  //     }
  //   }
  // }
}