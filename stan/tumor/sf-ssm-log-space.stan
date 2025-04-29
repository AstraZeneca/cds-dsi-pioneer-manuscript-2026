functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "sf-ssls_functions.stan"
}  

data {
  #include "../base_data.stan"
  
  int<lower = 0, upper = 1> fit_tumor_data;
  int<lower = 0> sf_rep_T;
  int<lower = 0, upper = 1> debug;
  int<lower = 0, upper = 1> pop_param_only;
  int<lower = 0, upper = 1> independ_long_process_noise;
  int<lower = 0, upper = 1> independ_cross_process_noise;
  int<lower = 1, upper = n_patients> train_patients_pos, train_patients_end;
  
  // GP parameters
  real<lower = 0> pop_tumor_gp_rho_meanlog;
  real<lower = 0> pop_tumor_gp_rho_sdlog;
  real<lower = 0> log_patient_tumor_gp_rho_sd_sd;
  
  // Process noise parameters
  real<lower = 0> pop_decrease_process_sd_sd;
  real<lower = 0> pop_growth_process_sd_sd;
  real<lower = 0> process_corr_param;
  real<lower = 0> measure_sd_sd;
  
  // Population rate parameters
  real pop_log_net_rate_mean;
  real<lower = 0> pop_log_net_rate_sd;
  real pop_log_rate_ratio_mean;
  real<lower = 0> pop_log_rate_ratio_sd;
  real<lower = 0> patient_log_net_rate_sd_sd;
  
  // Growth lag parameters
  real growth_lag_mean;
  real<lower = 0> growth_lag_sd;
  real<lower = 0> patient_log_growth_lag_sd_sd;
  real<lower = 0> log_growth_transition_rate_sd;
  
  // Correlation parameters
  real<lower = 0> rate_corr_param;
  
  // Proportion parameters
  real pop_decrease_prop_logis_mean;
  real<lower = 0> pop_decrease_prop_logis_sd;
  real<lower = 0> patient_decrease_prop_logis_sd_sd;
  real<lower = 0> log_lod_sd;
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
  
  int<lower = 1> n_total_visits_m1 = sum(n_patient_visits) - n_patients;
  array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
  
  int<lower = 1> n_train_patients = train_patients_end - train_patients_pos + 1;
  int<lower = 1> n_total_train_visits_m1 = sum(n_patient_visits[train_patients_pos:train_patients_end]) - n_train_patients;
  
  
  real log_lod = log(0.1);
}

parameters {
  real pop_log_net_rate;            // Population-level net rate (log(d-g))
  real<lower = 0.125> pop_log_rate_ratio; // Population-level ratio (log(d/g))

  // Patient-level variation for net rate only
  real<lower=0> patient_log_net_rate_sd;
  // vector<offset = pop_log_net_rate, multiplier = patient_log_net_rate_sd>[n_patients] patient_log_net_rate;
  vector[pop_param_only ? 0 : n_train_patients] raw_patient_log_net_rate;

  real pop_log_growth_lag;
  real pop_log_growth_transition_rate;

  real<lower = 0> patient_log_growth_lag_sd;
  // vector<offset = pop_log_growth_lag, multiplier = patient_log_growth_lag_sd>[n_patients] patient_log_growth_lag;
  vector[pop_param_only ? 0 : n_train_patients] raw_patient_log_growth_lag;
  
  // Patient-level GP
  // row_vector<lower = 0>[2] pop_tumor_gp_alpha;
  real log_pop_tumor_gp_rho;
  
  // Hierarchical length-scale parameter
  // vector<lower = 0>[2] log_trial_tumor_gp_rho_sd;
  // vector[add_trial_level_tumor_gp_param ? n_trials : 0] raw_log_trial_tumor_gp_rho_effect; 
  
  real<lower = 0> log_patient_tumor_gp_rho_sd;
  vector[independ_long_process_noise ? 0 : n_train_patients] raw_log_patient_tumor_gp_rho_effect; 
  
  matrix[n_total_train_visits_m1, 2] raw_patient_process_noise;

  // real<lower=0> pop_decrease_process_sd;
  // real<lower=0> pop_growth_process_sd;
  vector<lower = 0>[2] pop_process_sd;
  cholesky_factor_corr[independ_cross_process_noise ? 0 : 2] L_process_corr;
  
  real<lower=0> measure_sd;        // Measurement noise variance
  
  // real<lower = 0> patient_log_decrease_rate_sd;
  // real<lower = 0> patient_log_growth_rate_sd;
  
  real pop_decrease_prop_logis;

  real<lower = 0> patient_decrease_prop_logis_sd;
  // vector<offset = pop_decrease_prop_logis, multiplier = patient_decrease_prop_logis_sd>[n_patients] patient_decrease_prop_logis;
  vector[pop_param_only ? 0 : n_train_patients] raw_patient_decrease_prop_logis;
  
  // real log_lod;
  
  matrix[n_total_train_visits_m1, 2] raw_states; 
}

transformed parameters {
  vector[n_train_patients] patient_log_net_rate = rep_vector(pop_log_net_rate, n_train_patients);
  vector[n_train_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_train_patients);
  
  if (!pop_param_only) {
    patient_log_net_rate += patient_log_net_rate_sd * raw_patient_log_net_rate;
    patient_log_growth_lag += patient_log_growth_lag_sd * raw_patient_log_growth_lag;
  }
  
  vector[n_train_patients] patient_log_growth_rate = patient_log_net_rate - log_diff_exp(pop_log_rate_ratio, 0);
  vector[n_train_patients] patient_log_decrease_rate = patient_log_growth_rate + pop_log_rate_ratio;
  
  matrix[n_total_train_visits_m1, 2] patient_process_noise;
  
  vector[n_train_patients] patient_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_train_patients);
  
  if (!pop_param_only) {
    patient_decrease_prop_logis += patient_decrease_prop_logis_sd * raw_patient_decrease_prop_logis;
  }
  
  vector[n_train_patients] patient_log_decrease_prop = -log1p_exp(-patient_decrease_prop_logis);
  vector[n_train_patients] patient_log_growth_prop = patient_log_decrease_prop - patient_decrease_prop_logis;
  
  matrix[sum(n_patient_visits[train_patients_pos:train_patients_end]), 2] expected_states; 
  matrix[sum(n_patient_visits[train_patients_pos:train_patients_end]), 2] states; 
  
  profile("states") {
    vector[independ_long_process_noise ? 0 : n_train_patients] log_patient_tumor_gp_rho_effect; 
    vector[independ_long_process_noise ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
    vector[independ_long_process_noise ? 0 : n_train_patients] patient_tumor_gp_rho;
    
    if (!independ_long_process_noise) {
      log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
      log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
      patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial[train_patients_pos:train_patients_end]] + log_patient_tumor_gp_rho_effect);
    }
    
    for (i in train_patients_pos:train_patients_end) {
      int visit_pos, visit_end;
      (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
      int visit_m1_pos, visit_m1_end;
      (visit_m1_pos, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
      
      // array[n_patient_unique_visits[i]] real time_points = all_tumor_measure_t[get_int_sub_array(patient_unique_visits_idx, patient_unique_visits_pos, i)];
      array[n_patient_visits[i]] real time_points = t_patient_visits[visit_pos:visit_end];
      
      if (!independ_long_process_noise) {
        patient_process_noise[visit_m1_pos:visit_m1_end] = calc_gp_pred(
        // patient_obs_tumor_gp[patient_gp_pos:patient_gp_end] = ncp_gp_matern52(
          // all_tumor_measure_t[get_int_sub_array(patient_unique_visits_idx, patient_unique_visits_pos, i)], 
          time_points,
          patient_tumor_gp_rho[i], delta, 
          pop_process_sd, 
          independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
          raw_patient_process_noise[visit_m1_pos:visit_m1_end],
          1
        );
      } else {
        if (independ_cross_process_noise) {
          patient_process_noise[visit_m1_pos:visit_m1_end] = raw_patient_process_noise[visit_m1_pos:visit_m1_end];
        } else {
          patient_process_noise[visit_m1_pos:visit_m1_end] = raw_patient_process_noise[visit_m1_pos:visit_m1_end] * L_process_corr'; 
        }
        
        matrix[n_patient_visits[i], 2] scaled_process_sd = scale_process_sd(time_points, pop_process_sd);
        patient_process_noise[visit_m1_pos:visit_m1_end] = patient_process_noise[visit_m1_pos:visit_m1_end] .* scaled_process_sd[2:];
      }
  
      row_vector[2] initial_state = [patient_log_decrease_prop[i], patient_log_growth_prop[i]];
  
      (expected_states[visit_pos:visit_end], states[visit_pos:visit_end]) = sf_log_space_trajectory_ncp(
        raw_states[visit_m1_pos:visit_m1_end], initial_state, time_points,
        exp(patient_log_decrease_rate[i]), exp(patient_log_growth_rate[i]), exp(patient_log_growth_lag[i]), exp(pop_log_growth_transition_rate),
        patient_process_noise[visit_m1_pos:visit_m1_end]
      );
    }
  }
}

model {
  pop_log_net_rate ~ normal(pop_log_net_rate_mean, pop_log_net_rate_sd);
  pop_log_rate_ratio ~ normal(pop_log_rate_ratio_mean, pop_log_rate_ratio_sd);
  patient_log_net_rate_sd ~ normal(0, patient_log_net_rate_sd_sd);
  // patient_log_net_rate ~ normal(pop_log_net_rate, patient_log_net_rate_sd);
  raw_patient_log_net_rate ~ std_normal(); 

  pop_log_growth_lag ~ normal(growth_lag_mean, growth_lag_sd);
  pop_log_growth_transition_rate ~ normal(0, log_growth_transition_rate_sd);
  patient_log_growth_lag_sd ~ normal(0, patient_log_growth_lag_sd_sd);
  // patient_log_growth_lag ~ normal(pop_log_growth_lag, patient_log_growth_lag_sd);
  raw_patient_log_growth_lag ~ std_normal(); 
  
  // pop_tumor_gp_alpha ~ normal(0, pop_tumor_gp_alpha_sd);
  log_pop_tumor_gp_rho ~ normal(pop_tumor_gp_rho_meanlog, pop_tumor_gp_rho_sdlog);
  log_patient_tumor_gp_rho_sd ~ normal(0, log_patient_tumor_gp_rho_sd_sd);
  to_vector(raw_log_patient_tumor_gp_rho_effect) ~ std_normal(); 
  
  to_vector(raw_patient_process_noise) ~ std_normal();
  
  pop_process_sd[1] ~ normal(0, pop_decrease_process_sd_sd);
  pop_process_sd[2] ~ normal(0, pop_growth_process_sd_sd);
  
  measure_sd ~ normal(0, measure_sd_sd);

  if (!independ_cross_process_noise) {
    L_process_corr ~ lkj_corr_cholesky(process_corr_param);
  }
  
  pop_decrease_prop_logis ~ normal(pop_decrease_prop_logis_mean, pop_decrease_prop_logis_sd);
  patient_decrease_prop_logis_sd ~ normal(0, patient_decrease_prop_logis_sd_sd);
  // patient_decrease_prop_logis ~ normal(pop_decrease_prop_logis, patient_decrease_prop_logis_sd);
  raw_patient_decrease_prop_logis ~ std_normal();
  
  // log_lod ~ normal(log(lod), log_lod_sd);
  
  to_vector(raw_states) ~ std_normal();

  profile("loglik") { 
    if (fit_tumor_data) { 
      // for (i in 1:n_patients) {
      for (i in train_patients_pos:train_patients_end) {
        int visit_pos, visit_end;
        (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
        
        normalized_sld[visit_pos:visit_end] ~ sf_log_space_obs(states[visit_pos:visit_end], measure_sd, log_lod - log(sum_tumor_size[visit_pos]));
    
        if (debug) {    
          print(i, ": normalized_sld = ", normalized_sld[visit_pos:visit_end], ", exp(states) = ", exp(states));
        }
      }
    }
  }
}

generated quantities {
  corr_matrix[independ_cross_process_noise ? 0 : 2] process_corr;
  
  if (!independ_cross_process_noise) {
    process_corr = L_process_corr * L_process_corr';
  }
  
  vector<lower = 0, upper = 1>[max_t_width] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  matrix[max_t_width, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);
  
  // array[n_patients] matrix[2, sf_rep_T + 1] rep_state = rep_array(rep_matrix(0, 2, sf_rep_T + 1), n_patients);
  // array[n_patients] vector<lower = 0>[sf_rep_T + 1] rep_sld = rep_array(zeros_vector(sf_rep_T + 1), n_patients);
  // 
  // if (sf_rep_T > 0) {
  //   // for (i in 1:n_patients) {
  //   // for (i in linspaced_int_array(10, 1, n_patients)) {
  //   for (i in train_patients_pos:train_patients_end) {
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
  //         pop_decrease_process_sd, growth_process_sd, measure_sd
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