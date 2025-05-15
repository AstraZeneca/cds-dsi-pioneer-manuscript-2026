functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "sf-ssls_functions.stan"
  #include "recist.stanfunctions"
}  

data {
  #include "../base_data.stan"
  
  int<lower = 0, upper = 1> fit_tumor_data;
  int<lower = 0> sf_rep_T;
  int<lower = 0, upper = 1> debug;
  int<lower = 0, upper = 1> pop_growth_lag_param_only;
  int<lower = 0, upper = 1> pop_initial_states_param_only;
  int<lower = 0, upper = 1> pop_rates_param_only;
  int<lower = 0, upper = 1> pop_rho_param_only; 
  int<lower = 0, upper = 1> independ_long_process_noise;
  int<lower = 0, upper = 1> independ_cross_process_noise;
  int<lower = 0, upper = 1> forecast;
  int<lower = 0, upper = 1> run_parallel;
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
  
  // Define RECIST categories as integers
  int CR = 1;  // Complete Response
  int PR = 2;  // Partial Response  
  int SD = 3;  // Stable Disease
  int PD = 4;  // Progressive Disease
}

parameters {
  real pop_log_net_rate;            // Population-level net rate (log(d-g))
  real<lower = 0.125> pop_log_rate_ratio; // Population-level ratio (log(d/g))

  // Patient-level variation for net rate only
  real<lower=0> patient_log_net_rate_sd;
  // vector<offset = pop_log_net_rate, multiplier = patient_log_net_rate_sd>[n_patients] patient_log_net_rate;
  vector[pop_rates_param_only ? 0 : n_train_patients] raw_patient_log_net_rate;

  real pop_log_growth_lag;
  real pop_log_growth_transition_rate;

  real<lower = 0> patient_log_growth_lag_sd;
  // vector<offset = pop_log_growth_lag, multiplier = patient_log_growth_lag_sd>[n_patients] patient_log_growth_lag;
  vector[pop_growth_lag_param_only ? 0 : n_train_patients] raw_patient_log_growth_lag;
  
  // Patient-level GP
  // row_vector<lower = 0>[2] pop_tumor_gp_alpha;
  real log_pop_tumor_gp_rho;
  
  // Hierarchical length-scale parameter
  // vector<lower = 0>[2] log_trial_tumor_gp_rho_sd;
  // vector[add_trial_level_tumor_gp_param ? n_trials : 0] raw_log_trial_tumor_gp_rho_effect; 
  
  real<lower = 0> log_patient_tumor_gp_rho_sd;
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_train_patients] raw_log_patient_tumor_gp_rho_effect; 
  
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
  vector[pop_initial_states_param_only ? 0 : n_train_patients] raw_patient_decrease_prop_logis;
  
  // real log_lod;
}

transformed parameters {
  vector[n_train_patients] patient_log_net_rate_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_net_rate = rep_vector(pop_log_net_rate, n_train_patients);
  vector[n_train_patients] patient_log_growth_lag_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_train_patients);
  
  if (!pop_rates_param_only) {
    patient_log_net_rate_effect = patient_log_net_rate_sd * raw_patient_log_net_rate;
    patient_log_net_rate += patient_log_net_rate_effect;
  }
  
  if (!pop_growth_lag_param_only) {
    patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
    patient_log_growth_lag += patient_log_growth_lag_effect;
  }
  
  vector[n_train_patients] patient_log_growth_rate = patient_log_net_rate - log_diff_exp(pop_log_rate_ratio, 0);
  vector[n_train_patients] patient_log_decrease_rate = patient_log_growth_rate + pop_log_rate_ratio;
  
  vector[n_train_patients] patient_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_train_patients);
  
  if (!pop_initial_states_param_only) {
    patient_decrease_prop_logis += patient_decrease_prop_logis_sd * raw_patient_decrease_prop_logis;
  }
  
  vector[n_train_patients] patient_log_decrease_prop = -log1p_exp(-patient_decrease_prop_logis);
  vector[n_train_patients] patient_log_growth_prop = patient_log_decrease_prop - patient_decrease_prop_logis;
  
  matrix[sum(n_patient_visits[train_patients_pos:train_patients_end]), 2] states; 
  
  vector[independ_long_process_noise ? 0 : n_train_patients] patient_tumor_gp_rho;
  
  profile("states") {
    vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_train_patients] log_patient_tumor_gp_rho_effect; 
    vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
    
    if (!independ_long_process_noise) {
      if (pop_rho_param_only) {
        patient_tumor_gp_rho = rep_vector(exp(log_pop_tumor_gp_rho), n_train_patients);
      } else {
        log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
        log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
        patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial[train_patients_pos:train_patients_end]] + log_patient_tumor_gp_rho_effect);
      }
    }
    
    states = calc_states(
      create_pos(patient_visit_pos, train_patients_pos, train_patients_end),
      get_int_sub_array(t_patient_visits, patient_visit_pos, train_patients_pos, train_patients_end),
      independ_long_process_noise ? zeros_vector(n_train_patients) : patient_tumor_gp_rho,
      delta,
      pop_process_sd, 
      independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr, 
      raw_patient_process_noise,
      independ_long_process_noise, independ_cross_process_noise,
      append_col(patient_log_decrease_prop, patient_log_growth_prop),
      exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
      exp(patient_log_growth_lag), exp(pop_log_growth_transition_rate),
      run_parallel && !debug,
      debug
    ); 
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
  vector[n_train_patients] patient_log_growth_rate_residual = patient_log_growth_rate - (pop_log_net_rate - log_diff_exp(pop_log_rate_ratio, 0));
  real pop_log_growth_rate = pop_log_net_rate - log_diff_exp(pop_log_rate_ratio, 0);
  vector[n_train_patients] patient_log_decrease_rate_residual = patient_log_decrease_rate - (pop_log_growth_rate + pop_log_rate_ratio);
  vector[n_train_patients] patient_decrease_prop_residual = inv_logit(patient_decrease_prop_logis) - inv_logit(pop_decrease_prop_logis);
  
  corr_matrix[independ_cross_process_noise ? 0 : 2] process_corr;
  
  if (!independ_cross_process_noise) {
    process_corr = L_process_corr * L_process_corr';
  }
  
  vector<lower = 0, upper = 1>[max_t_width] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  matrix[max_t_width, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);
  matrix[n_total_train_visits_m1, 2] obs_patient_process_noise;
  matrix[forecast ? get_pos_total_size(forecast_visits_pos) : 0, 2] forecast_patient_process_noise;
  matrix[forecast ? get_pos_total_size(forecast_visits_pos) : 0, 2] forecast_patient_states;
  vector[forecast ? get_pos_total_size(forecast_visits_pos) : 0] forecast_patient_log_sld;
  array[forecast ? get_pos_total_size(forecast_visits_pos) : 0] int<lower = CR, upper = PD> forecast_recist;
  
  for (i in train_patients_pos:train_patients_end) {
    int visit_pos, visit_end;
    (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
    
    int visit_m1_start, visit_m1_end;
    (visit_m1_start, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
    
    int forecast_visit_start, forecast_visit_end;
    (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
    int forecast_size = forecast_visit_end - forecast_visit_start + 1; 
    
    obs_patient_process_noise[visit_m1_start:visit_m1_end] = calc_patient_process_noise(
      raw_patient_process_noise[visit_m1_start:visit_m1_end], get_int_sub_array(t_patient_visits, patient_visit_pos, i), patient_tumor_gp_rho[i], delta, pop_process_sd, L_process_corr,
      independ_long_process_noise, independ_cross_process_noise
    );
    
    array[n_patient_forecast_visits[i] + 1] int forecast_time = linspaced_int_array(n_patient_forecast_visits[i] + 1, patient_last_obs_visit[i], last_predict_visit);
    
    if (forecast && n_patient_forecast_visits[i] > 0) {
      // matrix[n_patient_visits[i], 2] curr_obs_states = sf_log_space_trajectory_ncp(
      //     [ patient_log_decrease_prop[i], patient_log_growth_prop[i] ],
      //     get_int_sub_array(t_patient_visits, patient_visit_pos, i),
      //     exp(patient_log_decrease_rate[i]), exp(patient_log_growth_rate[i]),
      //     exp(patient_log_growth_lag[i]), exp(pop_log_growth_transition_rate),
      //     obs_patient_process_noise[visit_m1_start:visit_m1_end]
      //   ).2;
      //   
      // for (t in 1:n_patient_visits[i]) {
      //   int dec = abs(states[visit_pos:visit_end][t, 1] - curr_obs_states[t, 1]) > 1e-6;
      //   int gro = abs(states[visit_pos:visit_end][t, 2] - curr_obs_states[t, 2]) > 1e-6;
      //   
      //   if (dec || gro) {
      //     print("initial_state = ", [ patient_log_decrease_prop[i], patient_log_growth_prop[i] ], ", time_points = ", get_int_sub_array(t_patient_visits, patient_visit_pos, i),
      //           ", dec rate = ", exp(patient_log_decrease_rate[i]), ", gro rate = ", exp(patient_log_growth_rate[i]), ", lag = ", exp(patient_log_growth_lag[i]), ", transit = ", exp(pop_log_growth_transition_rate));
      //     
      //     fatal_error(i, ": t = ", t, ", dec = ", dec, ", gro = ", gro, ", states[visit_pos:visit_end][t, 1] = ", states[visit_pos:visit_end][t, 1], ", curr_obs_states[t, 1] = ", curr_obs_states[t, 1],
      //     ", states[visit_pos:visit_end][t, 2] = ", states[visit_pos:visit_end][t, 2], ", curr_obs_states[t, 2] = ", curr_obs_states[t, 2]);
      //   }
      // }
      
      forecast_patient_process_noise[forecast_visit_start:forecast_visit_end] = multi_normal_rng(
        obs_patient_process_noise[visit_m1_start:visit_m1_end],
        get_int_sub_array(t_patient_visits, patient_visit_pos, i)[2:],
        forecast_time[2:],
        independ_long_process_noise ? 0 : patient_tumor_gp_rho[i],
        pop_process_sd,
        independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
        delta
      );
      
      forecast_patient_states[forecast_visit_start:forecast_visit_end] = sf_log_space_trajectory_ncp(
        states[patient_visit_pos[i + 1] - 1],
        forecast_time,
        exp(patient_log_decrease_rate[i]), exp(patient_log_growth_rate[i]),
        exp(patient_log_growth_lag[i]), exp(pop_log_growth_transition_rate),
        forecast_patient_process_noise[forecast_visit_start:forecast_visit_end]
      ).2[2:];
    
      forecast_patient_log_sld[forecast_visit_start:forecast_visit_end] = 
      to_vector(normal_rng(
        to_vector(log_sum_exp(
          forecast_patient_states[forecast_visit_start:forecast_visit_end, 1], forecast_patient_states[forecast_visit_start:forecast_visit_end, 2]
        )) + log(sum_tumor_size[visit_pos]),
        rep_vector(measure_sd, forecast_size)
      ));
      
      forecast_recist[forecast_visit_start:forecast_visit_end] = calculate_target_recist(
        append_row(sum_tumor_size[visit_pos], exp(forecast_patient_log_sld[forecast_visit_start:forecast_visit_end]))
      );
    }
  }
}

