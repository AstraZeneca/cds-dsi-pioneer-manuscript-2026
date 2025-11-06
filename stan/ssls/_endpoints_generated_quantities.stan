// Latent states //////////////////////////////////////////

matrix[n_total_visits_m1, 2] obs_patient_process_noise;
matrix[n_total_forecast_visits, 2] forecast_patient_process_noise;

matrix[n_total_forecast_visits, 2] forecast_patient_states;

vector[sum(n_patient_visits)] rep_mean_patient_log_sld, rep_patient_log_sld;
vector[n_total_forecast_visits] forecast_mean_patient_log_sld, forecast_patient_log_sld;

array[sum(n_patient_visits)] int<lower = CR, upper = PD + 1> rep_recist = rep_array(PD + 1, sum(n_patient_visits));
array[n_total_forecast_visits] int<lower = CR, upper = PD> forecast_recist;
  
// Endpoints (PFS, ORR, Median PFS, PFSn, ...) ////////////

array[n_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, spop_target_obs_cens_pfs; // spop_other_events_pfs, spop_pfs,  
array[n_patients] int<lower = 0, upper = 1> 
  sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored; // spop_other_events_right_censored, spop_right_censored; 

// Forecasting for right censored patients 
array[n_right_censored_patients] int<lower = 0> forecast_target_pfs; //, forecast_other_events_pfs, forecast_pfs; 
array[n_right_censored_patients] int<lower = 0, upper = 1> forecast_target_right_censored; //, forecast_other_events_right_censored, forecast_right_censored; 

array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_target_km_est, // sample_other_events_km_est,
                                                            spop_target_km_est, spop_target_obs_cens_km_est; // spop_other_events_km_est, spop_km_est, 

array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_target_km_est, 
  cond_spop_target_km_est, cond_spop_target_obs_cens_km_est; // cond_spop_other_events_km_est, cond_spop_km_est, 

array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n; 
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n;

array[n_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
array[n_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;

vector<lower = 0, upper = 1>[n_trials] sample_target_orr, spop_target_orr;
vector<lower = 0, upper = 1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group), cond_spop_target_orr = zeros_vector(n_cond_group);

// PFS quantiles (replaces median-only variablesHere are some guidelines for where to put different kinds of code:)
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_trials, n_pfs_quantiles] int sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

profile("gen_quant") {
  int right_censored_idx = 1;
  
  // Generate states for all patients at once
  (forecast_patient_states, rep_patient_log_sld, rep_mean_patient_log_sld,
   forecast_patient_log_sld, forecast_mean_patient_log_sld,
   obs_patient_process_noise, forecast_patient_process_noise) = 
    generate_all_patients_states_with_means_rng(
      states,
      patient_visit_pos,
      patient_visit_m1_pos,
      forecast_visits_pos,
      n_patient_visits,
      n_patient_screening_visits,
      n_patient_forecast_visits,
      patient_last_obs_visit,
      last_predict_visit,
      t_patient_visits,
      patient_log_decrease_rate,
      patient_log_growth_rate,
      sum_tumor_size,
      0.0001, 0.0001, // forecast_growth_lag, forecast_growth_transition
      measure_sd,
      independ_long_process_noise,
      independ_cross_process_noise,
      pop_process_sd,
      L_process_corr,
      log_pop_tumor_gp_rho,
      delta
    );
  
  // Calculate RECIST classifications for all patients at once
  (rep_recist, forecast_recist) = calculate_all_patients_recist(
    rep_mean_patient_log_sld,
    forecast_mean_patient_log_sld,
    patient_visit_pos,
    forecast_visits_pos,
    n_patient_visits,
    n_patient_screening_visits,
    n_patient_forecast_visits
  );
  
  // Calculate patient-level endpoints (PFS, response) for all patients at once
  (sample_target_pfs, sample_target_right_censored,
   spop_target_pfs, spop_target_right_censored,
   spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
   sample_target_confirmed_response, sample_target_unconfirmed_response,
   spop_target_confirmed_response, spop_target_unconfirmed_response,
   forecast_target_pfs, forecast_target_right_censored) =
    calculate_all_patients_endpoints(
      recist,
      rep_recist,
      forecast_recist,
      pfs,
      interval_censored,
      right_censored,
      patient_visit_pos,
      forecast_visits_pos,
      n_patient_visits,
      n_patient_screening_visits,
      n_patient_forecast_visits,
      patient_last_obs_visit,
      last_predict_visit,
      t_patient_visits,
      max_all_t
    );
    
  // Aggregate to trial-level metrics
  (sample_target_orr, spop_target_orr,
   sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
   sample_target_quant_pfs, spop_target_quant_pfs,
   sample_target_quant_pfs_exceeds_max, spop_target_quant_pfs_exceeds_max,
   sample_target_pfs_n, spop_target_pfs_n) =
    aggregate_trial_metrics(
      sample_target_confirmed_response,
      spop_target_confirmed_response,
      sample_target_pfs,
      sample_target_right_censored,
      spop_target_pfs,
      spop_target_right_censored,
      spop_target_obs_cens_pfs,
      spop_target_obs_cens_right_censored,
      trial_patient_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );
  
  // Aggregate to conditional group-level metrics
  (cond_sample_target_orr, cond_spop_target_orr,
   cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
   cond_sample_target_quant_pfs, cond_spop_target_quant_pfs,
   cond_sample_target_quant_pfs_exceeds_max, cond_spop_target_quant_pfs_exceeds_max,
   cond_sample_target_pfs_n, cond_spop_target_pfs_n) =
    aggregate_conditional_group_metrics(
      sample_target_confirmed_response,
      spop_target_confirmed_response,
      sample_target_pfs,
      sample_target_right_censored,
      spop_target_pfs,
      spop_target_right_censored,
      spop_target_obs_cens_pfs,
      spop_target_obs_cens_right_censored,
      cond_group,
      cond_group_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );
}