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
  int<lower = 0, upper = 1> pop_covar_coef_only;
  int<lower = 0, upper = 1> independ_long_process_noise;
  int<lower = 0, upper = 1> independ_cross_process_noise;
  int<lower = 0, upper = 1> run_parallel;
  int<lower = 1, upper = n_patients> train_patients_pos, train_patients_end;

  int<lower = 0, upper = 1> add_trial_level_net_rate; 
  int<lower = 0, upper = 1> add_trial_level_baseline_hazard;
  int<lower = 0, upper = 1> add_trial_level_prop;
  
  array[sum(n_patient_visits)] int<lower = 1, upper = 5> recist;
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive. The last week observed with no progression.
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  array[n_patients] int<lower = 0> interval_censored; // The number of weeks after `pfs` that actual progression could have happened. E.g., zero means progression happened the next week.
  
  array[n_patients] int<lower = 0> target_pfs; // PFS based on target tumor SLD only 
  array[n_patients] int<lower = 0, upper = 1> target_right_censored;
  
  array[n_patients] int<lower = 0> death_week;
 
  int<lower = 0> n_covar; 
  matrix[n_patients, n_covar] covar_design_matrix;
  
  int<lower = 0> n_pfs_timepoints;
  array[n_pfs_timepoints] int<lower = 0> pfs_timepoints; // In months
  
  #include "sf-ssls-hyperparam.stan" 
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
  #include "other_events_transformed_data.stan"
  
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
  
  array[n_trials + 1] int train_trial_right_censored_pos, train_trial_right_uncensored_pos;
 
  array[n_trials + 1] int train_trial_patient_pos = resize_pos(trial_patient_pos, train_patients_pos, train_patients_end);
  print("train_trial_patient_pos = ", train_trial_patient_pos);
  
  array[n_train_patients] int<lower = 1> train_patient_trial = patient_trial[train_patients_pos:train_patients_end];
 
  {
    int right_uncensored_idx = 1;
    array[n_trials] int n_train_trial_censored = zeros_int_array(n_trials), n_train_trial_uncensored = zeros_int_array(n_trials);
    
    for (i in train_patients_pos:train_patients_end) {
      if (right_censored[i]) {
        n_train_trial_censored[patient_trial[i]] += 1;
      } else {
        train_right_uncensored_patients[right_uncensored_idx] = i;
        right_uncensored_idx += 1;
        n_train_trial_uncensored[patient_trial[i]] += 1;
      }
    }
    
    train_trial_right_censored_pos = create_pos(n_train_trial_censored);
    train_trial_right_uncensored_pos = create_pos(n_train_trial_uncensored);
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
  
  // Define non-target status categories
  int NT_CR = 1;
  int NT_STABLE = 2;  // Non-CR/Non-PD
  int NT_PD = 3;
  
  array[n_pfs_timepoints] int<lower = 0> sorted_pfs_timepoints = sort_asc(pfs_timepoints); 
}

parameters {
  #include "other_events_parameters.stan"
  
  real pop_log_net_rate;            // Population-level net rate (log(d-g))
  real<lower = 0> pop_log_rate_ratio; // Population-level ratio (log(d/g))
  
  // Trial-level variation for net rate only
  real<lower=0> trial_log_net_rate_sd;
  vector[pop_rates_param_only || !add_trial_level_net_rate ? 0 : n_train_trials] raw_trial_log_net_rate;

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
  vector[pop_initial_states_param_only || !add_trial_level_prop ? 0 : n_train_trials] raw_trial_decrease_prop_logis;

  real<lower = 0> patient_decrease_prop_logis_sd;
  // vector<offset = pop_decrease_prop_logis, multiplier = patient_decrease_prop_logis_sd>[n_patients] patient_decrease_prop_logis;
  vector[pop_initial_states_param_only ? 0 : n_train_patients] raw_patient_decrease_prop_logis;
  
  // real log_lod;
  
  // Covariate effects on rates
  vector[n_covar] pop_log_net_rate_coef;      // Population-level covariate effects on net rate
  // vector[n_covar] pop_log_rate_ratio_coef;    // Population-level covariate effects on rate ratio
  
  // Optional: hierarchical covariate effects
  row_vector<lower=0>[pop_covar_coef_only ? 0 : n_covar] trial_log_net_rate_coef_sd;
  matrix[pop_covar_coef_only ? 0 : n_train_trials, n_covar] raw_trial_log_net_rate_coef;
  
  // real<lower=0> patient_log_net_rate_coef_sd;
  // matrix[pop_rates_param_only ? 0 : n_train_patients, n_covar] raw_patient_log_net_rate_coef;
}

transformed parameters {
  #include "other_events_transformed_parameters.stan"
  
  vector[n_train_trials] trial_log_net_rate_effect = zeros_vector(n_train_trials);
  vector[n_train_patients] patient_log_net_rate_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_net_rate = rep_vector(pop_log_net_rate, n_train_patients);
  // vector[n_train_patients] patient_log_rate_ratio_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_rate_ratio = rep_vector(pop_log_rate_ratio, n_train_patients);
  vector[n_train_patients] patient_log_growth_lag_effect = zeros_vector(n_train_patients);
  vector[n_train_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_train_patients);
  
  // Calculate linear predictors for rates
  vector[n_train_patients] patient_log_net_rate_linpred = covar_design_matrix[train_patients_pos:train_patients_end] * pop_log_net_rate_coef;
  // vector[n_train_patients] patient_log_rate_ratio_linpred = covar_design_matrix[train_patients_pos:train_patients_end] * pop_log_rate_ratio_coef;
  
  matrix[n_train_trials, n_covar] trial_log_net_rate_coef = rep_matrix(0, n_train_trials, n_covar);
  
  if (!pop_covar_coef_only) {
    // Trial-level covariate effects
    trial_log_net_rate_coef = rep_matrix(trial_log_net_rate_coef_sd, n_train_trials) .* raw_trial_log_net_rate_coef;
    
    patient_log_net_rate_linpred += rows_dot_product(covar_design_matrix[train_patients_pos:train_patients_end],
                                                     trial_log_net_rate_coef[train_patient_trial]);
    
    // Patient-level covariate effects (if you want this level of complexity)
    // matrix[n_train_patients, n_covar] patient_log_net_rate_coef = patient_log_net_rate_coef_sd * raw_patient_log_net_rate_coef;
    // patient_log_net_rate_lp += rows_dot_product(covar_design_matrix[train_patients_pos:train_patients_end], patient_log_net_rate_coef);
  }

  // Add hierarchical covariate effects if needed
  if (!pop_rates_param_only && add_trial_level_net_rate) {
    trial_log_net_rate_effect = trial_log_net_rate_sd * raw_trial_log_net_rate;
    patient_log_net_rate_effect = patient_log_net_rate_sd * raw_patient_log_net_rate;
    patient_log_net_rate += trial_log_net_rate_effect[train_patient_trial] + patient_log_net_rate_effect;
  }

  // Update the rate calculations to include covariate effects
  patient_log_net_rate += patient_log_net_rate_linpred + trial_log_net_rate_effect[train_patient_trial] + patient_log_net_rate_effect;
  // patient_log_rate_ratio = pop_log_rate_ratio + patient_log_rate_ratio_lp + patient_log_rate_ratio_effect;
  
  if (!pop_growth_lag_param_only) {
    patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
    patient_log_growth_lag += patient_log_growth_lag_effect;
  }
  
  vector[n_train_patients] patient_log_growth_rate = patient_log_net_rate - log_diff_exp(patient_log_rate_ratio, zeros_vector(n_train_patients));
  vector[n_train_patients] patient_log_decrease_rate = patient_log_growth_rate + patient_log_rate_ratio;
  
  vector[n_trials] trial_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_trials);
  vector[n_train_patients] patient_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_train_patients);
  
  if (!pop_initial_states_param_only && add_trial_level_prop) {
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
      run_parallel, // && !debug,
      0 // debug 
    );
  }
}

model {
  #include "other_events_priors.stan"
  #include "sf-ssls-priors.stan"

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
      
      for (s in 1:n_trials) {
        for (k in 1:n_causes) {
          target += sum(get_sub_vector(patient_response_lp[, k], train_trial_patient_pos, s));
        }
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
  
  vector<lower = 0, upper = 1>[max_all_t] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  
  matrix[max_all_t, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);
  
  matrix[n_total_train_visits_m1, 2] obs_patient_process_noise;
  matrix[n_total_train_forecast_visits, 2] forecast_patient_process_noise;
  
  matrix[n_total_train_forecast_visits, 2] forecast_patient_states;
  
  vector[n_total_train_visits] rep_patient_log_sld;
  vector[n_total_train_forecast_visits] forecast_patient_log_sld;
  
  array[n_total_train_visits] int<lower = CR, upper = PD + 1> rep_recist = rep_array(PD + 1, n_total_train_visits);
  array[n_total_train_forecast_visits] int<lower = CR, upper = PD> forecast_recist;
 
  array[n_train_patients] int<lower = 0> spop_target_pfs, spop_non_target_pfs, spop_pfs, spop_target_obs_cens_pfs; 
  array[n_train_patients] int<lower = 0, upper = 1> spop_target_right_censored, spop_non_target_right_censored, spop_right_censored, spop_target_obs_cens_right_censored; 
  // Forecasting for right censored patients 
  array[n_train_right_censored_patients] int<lower = 0> forecast_target_pfs, forecast_non_target_pfs, forecast_pfs; // Zero means right censored
  array[n_train_right_censored_patients] int<lower = 0, upper = 1> forecast_target_right_censored, forecast_non_target_right_censored, forecast_right_censored; 
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_km_est; // sample_target_km_est, sample_non_target_km_est,
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] spop_target_km_est, spop_non_target_km_est, spop_km_est, spop_target_obs_cens_km_est;
  
  array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] forecast_target_pfs_n; //, forecast_pfs_n;
    
  array[n_train_patients] int<lower = 0, upper = 1> forecast_confirmed_response;
  vector<lower = 0, upper = 1>[n_trials] forecast_target_orr;
  
  {
    int right_censored_idx = 1;
    
    for (i in train_patients_pos:train_patients_end) {
      int train_idx = i - train_patients_pos + 1;
  
      int train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end;
      (train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end) = get_visit_pos(train_patient_visit_pos, train_idx, n_patient_screening_visits[i]);
      int train_visit_size = get_pos_size(train_patient_visit_pos, train_idx);
  
      int train_visit_m1_start,  train_visit_m1_end;
      (train_visit_m1_start, train_visit_m1_end) = get_pos(train_patient_visit_m1_pos, train_idx);
      int train_visit_m1_size = get_pos_size(train_patient_visit_m1_pos, train_idx);
  
      int train_forecast_visit_start, train_forecast_visit_end;
      (train_forecast_visit_start, train_forecast_visit_end) = get_pos(train_forecast_visits_pos, train_idx);
  
      int visit_pos, visit_screening_end, visit_treat_pos, visit_end;
      (visit_pos, visit_screening_end, visit_treat_pos, visit_end) = get_visit_pos(patient_visit_pos, i, n_patient_screening_visits[i]);
      
      int visit_m1_start, visit_m1_end;
      (visit_m1_start, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
  
      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = get_pos_size(forecast_visits_pos, i);
      
      array[n_patient_visits[i]] int curr_visits = get_int_sub_array(train_patient_visits, train_patient_visit_pos, train_idx);
      
      int n_obs_treat_visits = train_visit_size - n_patient_screening_visits[i];
    
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
        
      array[n_patient_forecast_visits[i] + 1] int forecast_time = linspaced_int_array(n_patient_forecast_visits[i] + 1, patient_last_obs_visit[i], last_predict_visit);
     
      (spop_non_target_pfs[train_idx], spop_non_target_right_censored[train_idx]) = survival_time_rng(log_cond_prob_surv[1, train_idx]);
      
      if (n_patient_forecast_visits[i] > 0) {
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
      }
      
      array[n_obs_treat_visits + forecast_size] int full_predict_recist = calculate_target_recist(
        append_row(exp(rep_patient_log_sld[train_visit_start:train_visit_end]), exp(forecast_patient_log_sld[train_forecast_visit_start:train_forecast_visit_end])) * 10,
        n_patient_screening_visits[i]
      );
      
      rep_recist[train_treat_visit_start:train_visit_end] = full_predict_recist[:n_obs_treat_visits]; 
      forecast_recist[train_forecast_visit_start:train_forecast_visit_end] = full_predict_recist[(n_obs_treat_visits + 1):];
          
      if (n_patient_forecast_visits[i] > 0) {
        forecast_confirmed_response[train_idx] = find_first(
          append_array(recist[visit_pos:visit_end], forecast_recist[train_forecast_visit_start:train_forecast_visit_end]), { PR, CR }, 2) > 0;
         
        if (right_censored[i]) {
          (forecast_non_target_pfs[right_censored_idx], forecast_non_target_right_censored[right_censored_idx]) = survival_time_rng(
            log_cond_prob_surv[1, train_idx], pfs[i] + interval_censored[i], right_censored[i], 0
          );
          
          forecast_target_pfs[right_censored_idx] = find_first(forecast_recist[train_forecast_visit_start:train_forecast_visit_end], PD);
          forecast_target_right_censored[right_censored_idx] = forecast_target_pfs[right_censored_idx] == 0;
          
          if (!forecast_target_right_censored[right_censored_idx]) {
            forecast_target_pfs[right_censored_idx] = forecast_time[forecast_target_pfs[right_censored_idx]];
          } else {
            forecast_target_pfs[right_censored_idx] = max_all_t; 
          }
          
          forecast_pfs[right_censored_idx] = min(forecast_target_pfs[right_censored_idx], forecast_non_target_pfs[right_censored_idx]);
          forecast_right_censored[right_censored_idx] = forecast_target_right_censored[right_censored_idx] && forecast_non_target_right_censored[right_censored_idx];
          
          right_censored_idx += 1;
        }
        
        spop_target_pfs[train_idx] = find_first(append_array(rep_recist[train_treat_visit_start:train_visit_end], 
                                                             forecast_recist[train_forecast_visit_start:train_forecast_visit_end]), 
                                        PD); 
      } else {
        spop_target_pfs[train_idx] = find_first(rep_recist[train_treat_visit_start:train_visit_end], PD); 
      }
      
      // Why add one? We're passing the recist array excluding the first one.
      spop_target_right_censored[train_idx] = spop_target_pfs[train_idx] == 0;
      
      if (!spop_target_right_censored[train_idx]) { 
        // Why add one? We're passing the recist array excluding the screening visits.
        spop_target_pfs[train_idx] += n_patient_screening_visits[i];
        
        spop_target_pfs[train_idx] = spop_target_pfs[train_idx] <= n_train_patient_visits[train_idx] ? 
                              curr_visits[spop_target_pfs[train_idx]] :
                              forecast_time[spop_target_pfs[train_idx] - n_train_patient_visits[train_idx]]; 
      } else {
        spop_target_pfs[train_idx] = max_all_t;
      }
     
      // spop_target_obs_cens_right_censored[train_idx] = (spop_target_pfs[train_idx] > pfs[i] && right_censored[i]) || spop_target_right_censored[train_idx]; 
      spop_target_obs_cens_right_censored[train_idx] = right_censored[i]; 
      spop_target_obs_cens_pfs[train_idx] = right_censored[i] ? min(spop_target_pfs[train_idx], pfs[i]) : spop_target_pfs[train_idx];
      
      spop_pfs[train_idx] = min(spop_non_target_pfs[train_idx] + 1, 
                                max(0, spop_target_pfs[train_idx])); // BUG a couple of patients end up with negative weeks. We need to figure out why.
      spop_right_censored[train_idx] = spop_target_right_censored[train_idx] && spop_non_target_right_censored[train_idx]; 
    }
      
    for (s in 1:n_trials) {
      if (get_pos_size(train_trial_patient_pos, s) > 0) {
        int n_curr_uncensored_obs = get_pos_size(train_trial_right_uncensored_pos, s);
        array[n_curr_uncensored_obs] int curr_uncensored_obs = get_int_sub_array(train_right_uncensored_patients, train_trial_right_uncensored_pos, s);
        
        sample_km_est[s] = estimate_kaplan_meier(
          append_array(ub_pfs_p1[curr_uncensored_obs], get_int_sub_array(forecast_pfs, train_trial_right_censored_pos, s)), 
          append_array(right_censored[curr_uncensored_obs], get_int_sub_array(forecast_right_censored, train_trial_right_censored_pos, s)), 
          max_all_t).1; 
           
        for (n in 1:n_pfs_timepoints) {
          forecast_target_pfs_n[s, n] = sample_km_est[s, pfs_timepoints[n] * 4]; 
        }
                                             
        forecast_target_orr[s] = mean(get_int_sub_array(forecast_confirmed_response, train_trial_patient_pos, s));
        
        spop_target_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_target_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_target_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 0).1; 
                                               
        spop_target_obs_cens_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_target_obs_cens_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_target_obs_cens_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 0).1; 
                                               
        spop_non_target_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_non_target_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_non_target_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 1).1; 
        
        spop_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 0).1; 
      } else {
        spop_target_km_est[s] = zeros_vector(max_all_t + 1);
        spop_target_obs_cens_km_est[s] = zeros_vector(max_all_t + 1);
        spop_non_target_km_est[s] = zeros_vector(max_all_t + 1);
        spop_km_est[s] = zeros_vector(max_all_t + 1);
        sample_km_est[s] = zeros_vector(max_all_t + 1);
      }
    }
  }
  
  #include "sf-ssls-accuracy_gen_quant.stan"
}
