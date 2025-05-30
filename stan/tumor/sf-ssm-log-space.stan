functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
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
  
  array[sum(n_patient_visits)] int<lower = 1, upper = 5> recist;
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive. The last week observed with no progression.
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  array[n_patients] int<lower = 0> interval_censored; // The number of weeks after `pfs` that actual progression could have happened. E.g., zero means progression happened the next week.
  
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
  real<lower = 0> trial_log_net_rate_sd_sd;
  real<lower = 0> patient_log_net_rate_sd_sd;
  real<lower = 0> patient_log_rate_ratio_sd_sd;
  
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
  real<lower = 0> trial_decrease_prop_logis_sd_sd;
  real<lower = 0> patient_decrease_prop_logis_sd_sd;
  
  real<lower = 0> log_lod_sd;
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
  
  int<lower = 1> n_total_visits_m1 = sum(n_patient_visits) - n_patients;
  array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
 
  int<lower = 1, upper = n_trials> n_train_trials = max(patient_trial) - min(patient_trial) + 1; 
  int<lower = 1> n_train_patients = train_patients_end - train_patients_pos + 1;
  array[n_train_patients] int<lower = 0> n_train_patient_visits = n_patient_visits[train_patients_pos:train_patients_end];
  int<lower = 1> n_total_train_visits = sum(n_train_patient_visits);
  int<lower = 1> n_total_train_visits_m1 = n_total_train_visits - n_train_patients;
  
  array[n_total_train_visits] int train_patient_visits = get_int_sub_array(t_patient_visits, patient_visit_pos, train_patients_pos, train_patients_end);
  
  array[n_train_patients + 1] int<lower = 1> train_patient_visit_pos = create_pos(n_train_patient_visits);
  array[n_train_patients + 1] int<lower = 1> train_patient_visit_m1_pos = create_pos(n_train_patient_visits, -1);
  array[n_train_patients + 1] int<lower = 1> train_forecast_visits_pos = create_pos(n_patient_forecast_visits[train_patients_pos:train_patients_end]);
  
  int<lower = 1> n_total_train_forecast_visits = get_pos_total_size(train_forecast_visits_pos);
  
  array[n_total_train_visits] int train_obs_recist = get_int_sub_array(recist, patient_visit_pos, train_patients_pos, train_patients_end);
  
  int<lower = 0, upper = n_patients> n_train_right_censored_patients = sum(right_censored[train_patients_pos:train_patients_end]);
  int<lower = 0, upper = n_patients> n_train_right_uncensored_patients = n_patients - n_train_right_censored_patients; 
  array[n_train_right_uncensored_patients] int<lower = 1, upper = n_patients> train_right_uncensored_patients;
 
  {
    int right_uncensored_idx = 1;
    
    for (i in train_patients_pos:train_patients_end) {
      if (right_censored[i]) {
      } else {
        train_right_uncensored_patients[right_uncensored_idx] = i;
        right_uncensored_idx += 1;
      }
    }
  }
  
  array[n_patients] int<lower = 1> ub_pfs_p1;
  
  for (i in 1:n_patients) {
    ub_pfs_p1[i] = pfs[i] + interval_censored[i] + 1; 
  }
  
  real log_lod = log(0.1);
  
  // Define RECIST categories as integers
  int CR = 1;  // Complete Response
  int PR = 2;  // Partial Response  
  int SD = 3;  // Stable Disease
  int PD = 4;  // Progressive Disease
}

parameters {
  real pop_log_net_rate;            // Population-level net rate (log(d-g))
  real<lower = 0> pop_log_rate_ratio; // Population-level ratio (log(d/g))
  
  // Trial-level variation for net rate only
  real<lower=0> trial_log_net_rate_sd;
  vector[pop_rates_param_only ? 0 : n_train_trials] raw_trial_log_net_rate;

  // Patient-level variation for net rate only
  real<lower=0> patient_log_net_rate_sd;
  // vector<offset = pop_log_net_rate, multiplier = patient_log_net_rate_sd>[n_patients] patient_log_net_rate;
  vector[pop_rates_param_only ? 0 : n_train_patients] raw_patient_log_net_rate;
 
  // real<lower=0> patient_log_rate_ratio_sd;
  // vector[pop_rates_param_only ? 0 : n_train_patients] raw_patient_log_rate_ratio; 

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
  
  real<lower = 0> trial_decrease_prop_logis_sd;
  vector[pop_initial_states_param_only ? 0 : n_train_trials] raw_trial_decrease_prop_logis;

  real<lower = 0> patient_decrease_prop_logis_sd;
  // vector<offset = pop_decrease_prop_logis, multiplier = patient_decrease_prop_logis_sd>[n_patients] patient_decrease_prop_logis;
  vector[pop_initial_states_param_only ? 0 : n_train_patients] raw_patient_decrease_prop_logis;
  
  // real log_lod;
}

transformed parameters {
  vector[n_train_trials] trial_log_net_rate_effect = zeros_vector(n_train_trials);
  vector[n_train_patients] patient_log_net_rate_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_net_rate = rep_vector(pop_log_net_rate, n_train_patients);
  vector[n_train_patients] patient_log_rate_ratio_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_rate_ratio = rep_vector(pop_log_rate_ratio, n_train_patients);
  vector[n_train_patients] patient_log_growth_lag_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_train_patients);
  
  if (!pop_rates_param_only) {
    trial_log_net_rate_effect = trial_log_net_rate_sd * raw_trial_log_net_rate;
    patient_log_net_rate_effect = patient_log_net_rate_sd * raw_patient_log_net_rate;
    patient_log_net_rate += trial_log_net_rate_effect[patient_trial[train_patients_pos:train_patients_end]] + patient_log_net_rate_effect;
    
    // patient_log_rate_ratio_effect = patient_log_rate_ratio_sd * raw_patient_log_rate_ratio;
    // patient_log_rate_ratio += patient_log_rate_ratio_effect;
  }
  
  if (!pop_growth_lag_param_only) {
    patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
    patient_log_growth_lag += patient_log_growth_lag_effect;
  }
  
  vector[n_train_patients] patient_log_growth_rate = patient_log_net_rate - log_diff_exp(patient_log_rate_ratio, zeros_vector(n_train_patients));
  vector[n_train_patients] patient_log_decrease_rate = patient_log_growth_rate + patient_log_rate_ratio;
  
  vector[n_trials] trial_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_trials);
  vector[n_train_patients] patient_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_train_patients);
  
  if (!pop_initial_states_param_only) {
    trial_decrease_prop_logis += trial_decrease_prop_logis_sd * raw_trial_decrease_prop_logis;
    patient_decrease_prop_logis += 
      raw_trial_decrease_prop_logis[patient_trial[train_patients_pos:train_patients_end]] +
      patient_decrease_prop_logis_sd * raw_patient_decrease_prop_logis;
  }
  
  vector[n_train_patients] patient_log_decrease_prop = -log1p_exp(-patient_decrease_prop_logis);
  vector[n_train_patients] patient_log_growth_prop = patient_log_decrease_prop - patient_decrease_prop_logis;
  
  vector[n_train_patients] patient_tumor_gp_rho = independ_long_process_noise ? zeros_vector(n_train_patients) : rep_vector(exp(log_pop_tumor_gp_rho), n_train_patients);  
  
  matrix[n_total_train_visits, 2] states; 
  
  profile("states") {
    vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_train_patients] log_patient_tumor_gp_rho_effect; 
    vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
    
    if (!independ_long_process_noise && !pop_rho_param_only) {
      log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
      log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
      patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial[train_patients_pos:train_patients_end]] + log_patient_tumor_gp_rho_effect);
    }

    states = calc_states(
      train_patient_visit_pos,
      train_patient_visits,
      patient_tumor_gp_rho,
      delta,
      pop_process_sd,
      independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
      // raw_patient_process_noise,
      independ_long_process_noise, independ_cross_process_noise,
      append_col(patient_log_decrease_prop, patient_log_growth_prop),
      exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
      rep_vector(0.0001, n_train_patients), // exp(patient_log_growth_lag), 
      0.0001, // exp(pop_log_growth_transition_rate),
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
  // patient_log_rate_ratio_sd ~ normal(0, patient_log_rate_ratio_sd_sd);
  // raw_patient_log_rate_ratio ~ std_normal();
  trial_log_net_rate_sd ~ normal(0, trial_log_net_rate_sd_sd);
  raw_trial_log_net_rate ~ std_normal(); 

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
  trial_decrease_prop_logis_sd ~ normal(0, trial_decrease_prop_logis_sd_sd);
  raw_trial_decrease_prop_logis ~ std_normal();
  patient_decrease_prop_logis_sd ~ normal(0, patient_decrease_prop_logis_sd_sd);
  // patient_decrease_prop_logis ~ normal(pop_decrease_prop_logis, patient_decrease_prop_logis_sd);
  raw_patient_decrease_prop_logis ~ std_normal();
  
  // log_lod ~ normal(log(lod), log_lod_sd);
  
  profile("loglik") { 
    if (fit_tumor_data) { 
      for (i in train_patients_pos:train_patients_end) {
        int train_idx = i - train_patients_pos + 1;
        
        int train_visit_start, train_visit_end;
        (train_visit_start, train_visit_end) = get_pos(train_patient_visit_pos, train_idx);
        
        int visit_pos, visit_end;
        (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
        
        normalized_sld[visit_pos:visit_end] ~ sf_log_space_obs(states[train_visit_start:train_visit_end], measure_sd, log_lod - log(sum_tumor_size[visit_pos]));
      }
    }
  }
}

generated quantities {
  real pop_log_growth_rate = pop_log_net_rate - log_diff_exp(pop_log_rate_ratio, 0);
  vector[n_trials] trial_log_growth_rate = pop_log_growth_rate + trial_log_net_rate_effect;
  vector[n_trials] trial_log_growth_rate_residual = trial_log_growth_rate - pop_log_growth_rate;
  vector[n_train_patients] patient_log_growth_rate_residual = patient_log_growth_rate - trial_log_growth_rate[patient_trial[train_patients_pos:train_patients_end]];
  
  real pop_log_decrease_rate = pop_log_growth_rate + pop_log_rate_ratio;
  vector[n_trials] trial_log_decrease_rate = trial_log_growth_rate + pop_log_rate_ratio; 
  vector[n_trials] trial_log_decrease_rate_residual = trial_log_decrease_rate - pop_log_decrease_rate;
  vector[n_train_patients] patient_log_decrease_rate_residual = patient_log_decrease_rate - trial_log_decrease_rate[patient_trial[train_patients_pos:train_patients_end]];
  
  vector[n_train_patients] patient_decrease_prop_residual = inv_logit(patient_decrease_prop_logis) - inv_logit(trial_decrease_prop_logis[patient_trial[train_patients_pos:train_patients_end]]);
  vector[n_trials] trial_decrease_prop_residual = inv_logit(trial_decrease_prop_logis) - inv_logit(pop_decrease_prop_logis);
  real pop_log_decrease_prop = -log1p_exp(-pop_decrease_prop_logis);  
  real pop_log_growth_prop = pop_log_decrease_prop - pop_decrease_prop_logis;
  vector[n_trials] trial_log_decrease_prop = -log1p_exp(-trial_decrease_prop_logis);
  vector[n_trials] trial_log_growth_prop = trial_log_decrease_prop - trial_decrease_prop_logis;
  
  corr_matrix[independ_cross_process_noise ? 0 : 2] process_corr;
  
  if (!independ_cross_process_noise) {
    process_corr = L_process_corr * L_process_corr';
  }
  
  vector<lower = 0, upper = 1>[max_t_width] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  
  matrix[max_t_width, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);
  
  matrix[n_total_train_visits_m1, 2] obs_patient_process_noise;
  matrix[forecast ? n_total_train_forecast_visits : 0, 2] forecast_patient_process_noise;
  
  matrix[forecast ? n_total_train_forecast_visits : 0, 2] forecast_patient_states;
  
  vector[n_total_train_visits] rep_patient_log_sld;
  vector[forecast ? n_total_train_forecast_visits : 0] forecast_patient_log_sld;
  
  array[n_total_train_visits] int<lower = CR, upper = PD + 1> rep_recist;
  array[forecast ? n_total_train_forecast_visits : 0] int<lower = CR, upper = PD> forecast_recist;
 
  array[n_train_patients] int<lower = 0> spop_pfs; // Zero means right censored
  array[n_train_patients] int<lower = 0, upper = 1> spop_right_censored; 
  // Forecasting for right censored patients 
  array[forecast ? n_train_right_censored_patients : 0] int<lower = 0> forecast_pfs; // Zero means right censored
  array[forecast ? n_train_right_censored_patients : 0] int<lower = 0, upper = 1> forecast_right_censored; 
  vector<lower = 0, upper = 1>[max_all_t + 1] sample_km_est = zeros_vector(max_all_t + 1), spop_km_est = zeros_vector(max_all_t + 1);  
  
  // RECIST prediction accuracy metrics
  int<lower=0> correct_recist_predictions = 0;
  int<lower=0> total_recist_predictions = 0;
  matrix[4, 4] recist_confusion_matrix = rep_matrix(0, PD, PD); // rows = observed, cols = predicted
  real weighted_recist_accuracy_linear = 0;
  real weighted_recist_accuracy_quadratic = 0;
  int<lower=0> correct_recist_response_class = 0;
  int<lower=0> correct_recist_disease_control = 0;
  
  // Per-category metrics
  vector[PD] recist_category_sensitivity = zeros_vector(PD); // true positive rate per category
  vector[PD] recist_category_precision = zeros_vector(PD);   // positive predictive value per category
  vector[PD] recist_category_counts = zeros_vector(PD);      // number of observations per category

  {
    int right_censored_idx = 1;
    
    for (i in train_patients_pos:train_patients_end) {
      int train_idx = i - train_patients_pos + 1;
  
      int train_visit_start, train_visit_end;
      (train_visit_start, train_visit_end) = get_pos(train_patient_visit_pos, train_idx);
  
      int train_visit_m1_start, train_visit_m1_end;
      (train_visit_m1_start, train_visit_m1_end) = get_pos(train_patient_visit_m1_pos, train_idx);
      int train_visit_m1_size = train_visit_m1_end - train_visit_m1_start + 1;
  
      int train_forecast_visit_start, train_forecast_visit_end;
      (train_forecast_visit_start, train_forecast_visit_end) = get_pos(train_forecast_visits_pos, train_idx);
  
      int visit_pos, visit_end;
      (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
  
      int visit_m1_start, visit_m1_end;
      (visit_m1_start, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
  
      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = forecast_visit_end - forecast_visit_start + 1;
      
      obs_patient_process_noise[train_visit_m1_start:train_visit_m1_end] = rep_matrix(0, train_visit_m1_size, 2); 
      // obs_patient_process_noise[train_visit_m1_start:train_visit_m1_end] = calc_patient_process_noise(
      //   raw_patient_process_noise[train_visit_m1_start:train_visit_m1_end],
      //   get_int_sub_array(train_patient_visits, train_patient_visit_pos, train_idx),
      //   exp(log_pop_tumor_gp_rho), delta,
      //   pop_process_sd, L_process_corr,
      //   independ_long_process_noise, independ_cross_process_noise
      // );
      
      rep_patient_log_sld[train_visit_start] = log(sum_tumor_size[visit_pos]);
      rep_patient_log_sld[(train_visit_start + 1):train_visit_end] =
        to_vector(normal_rng(
          to_vector(log_sum_exp(states[(train_visit_start + 1):train_visit_end, 1], states[(train_visit_start + 1):train_visit_end, 2])) + log(sum_tumor_size[visit_pos]),
          rep_vector(measure_sd, n_patient_visits[i] - 1)
        ));
        
      rep_recist[train_visit_start] = PD + 1; 
      rep_recist[(train_visit_start + 1):train_visit_end] = calculate_target_recist(exp(rep_patient_log_sld[train_visit_start:train_visit_end]) * 10);
      
      array[n_patient_forecast_visits[i] + 1] int forecast_time = linspaced_int_array(n_patient_forecast_visits[i] + 1, patient_last_obs_visit[i], last_predict_visit);
      
      if (forecast && n_patient_forecast_visits[i] > 0) {
        // assert_matching_states(
        //   states[visit_pos:visit_end], 
        //   [ patient_log_decrease_prop[train_idx], patient_log_growth_prop[train_idx] ],
        //   get_int_sub_array(t_patient_visits, patient_visit_pos, i),
        //   exp(patient_log_decrease_rate[train_idx]), exp(patient_log_growth_rate[train_idx]),
        //   0.0001, // exp(patient_log_growth_lag[i]), 
        //   0.0001, // exp(pop_log_growth_transition_rate),
        //   rep_matrix(0, get_pos_size(train_patient_visit_m1_pos, train_idx), 2),
        //   debug
        // );
        
        if (independ_long_process_noise) {
          forecast_patient_process_noise[train_forecast_visit_start:train_forecast_visit_end] = multi_normal_rng(
            forecast_size, pop_process_sd, independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr
          );
        } else {
          forecast_patient_process_noise[train_forecast_visit_start:train_forecast_visit_end] = multi_normal_rng(
            obs_patient_process_noise[train_visit_m1_start:train_visit_m1_end],
            get_int_sub_array(t_patient_visits, patient_visit_pos, i)[2:],
            forecast_time[2:],
            exp(log_pop_tumor_gp_rho),
            pop_process_sd,
            independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
            delta
          );
        }
        
        forecast_patient_states[train_forecast_visit_start:train_forecast_visit_end] = sf_log_space_trajectory_ncp(
          states[train_patient_visit_pos[train_idx + 1] - 1],
          forecast_time,
          exp(patient_log_decrease_rate[train_idx]), exp(patient_log_growth_rate[train_idx]),
          0.0001, // exp(patient_log_growth_lag[train_idx]), 
          0.0001, // exp(pop_log_growth_transition_rate),
          rep_matrix(0, forecast_size, 2)
          // forecast_patient_process_noise[train_forecast_visit_start:train_forecast_visit_end]
        ).2[2:];
  
        forecast_patient_log_sld[train_forecast_visit_start:train_forecast_visit_end] =
          to_vector(normal_rng(
            to_vector(log_sum_exp(
              forecast_patient_states[train_forecast_visit_start:train_forecast_visit_end, 1], forecast_patient_states[train_forecast_visit_start:train_forecast_visit_end, 2]
            )) + log(sum_tumor_size[visit_pos]),
            rep_vector(measure_sd, forecast_size)
          ));
  
        forecast_recist[train_forecast_visit_start:train_forecast_visit_end] = calculate_target_recist(
          append_row(sum_tumor_size[visit_pos], exp(forecast_patient_log_sld[train_forecast_visit_start:train_forecast_visit_end])) * 10,
          // min(rep_patient_log_sld[train_visit_start:train_visit_end]) * 10
          min(sum_tumor_size[visit_pos:visit_end]) * 10
        );
        
        if (right_censored[i]) {
          forecast_pfs[right_censored_idx] = find_first(forecast_recist[train_forecast_visit_start:train_forecast_visit_end], PD);
          forecast_right_censored[right_censored_idx] = forecast_pfs[right_censored_idx] == 0;
          
          right_censored_idx += 1;
        }
        
        spop_pfs[train_idx] = find_first(append_array(rep_recist[(train_visit_start + 1):train_visit_end], forecast_recist[train_forecast_visit_start:train_forecast_visit_end]), PD);
      } else {
        spop_pfs[train_idx] = find_first(rep_recist[(train_visit_start + 1):train_visit_end], PD);
      }
      
      spop_right_censored[train_idx] = spop_pfs[train_idx] == 0;
      
      for (t in train_visit_start:train_visit_end) {
        if (train_obs_recist[t] <= PD) {
          // Update all metrics using the function
          (correct_recist_predictions, recist_confusion_matrix, recist_category_counts,
           weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
           correct_recist_response_class, correct_recist_disease_control) = update_recist_metrics(
            train_obs_recist[t], rep_recist[t],
            correct_recist_predictions, recist_confusion_matrix, recist_category_counts,
            weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
            correct_recist_response_class, correct_recist_disease_control
          );
          
          total_recist_predictions += 1;
        } 
      }
    }
   
    if (forecast) { 
      sample_km_est = estimate_kaplan_meier(
        append_array(ub_pfs_p1[train_right_uncensored_patients], forecast_pfs), 
        append_array(right_censored[train_right_uncensored_patients], forecast_right_censored), 
        max_all_t).1; 
    }
    
    spop_km_est = estimate_kaplan_meier(spop_pfs, spop_right_censored, max_all_t).1; 
  }
  
  // Calculate final metrics
  real recist_accuracy;
  real recist_response_accuracy;
  real recist_disease_control_accuracy;
  real recist_response_sensitivity;
  real recist_response_specificity;
  real recist_progression_sensitivity;
  real recist_progression_specificity;
  
  (recist_accuracy, recist_response_accuracy, recist_disease_control_accuracy,
   weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
   recist_category_sensitivity, recist_category_precision,
   recist_response_sensitivity, recist_response_specificity,
   recist_progression_sensitivity, recist_progression_specificity) = calculate_recist_summary_metrics(
    correct_recist_predictions, total_recist_predictions,
    recist_confusion_matrix, recist_category_counts,
    weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
    correct_recist_response_class, correct_recist_disease_control
  );
}

