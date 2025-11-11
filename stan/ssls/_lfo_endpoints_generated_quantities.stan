// LFO Endpoints Generated Quantities
// This version uses cutoff-aware COMPACT data (only n_cutoff_observed_patients)

// Latent states for cutoff-observed patients only
matrix[n_cutoff_visits_m1, 2] cutoff_obs_patient_process_noise;
matrix[n_cutoff_total_forecast_visits, 2] cutoff_forecast_patient_process_noise;

matrix[n_cutoff_total_forecast_visits, 2] cutoff_forecast_patient_states;

vector[n_cutoff_visits] cutoff_rep_mean_patient_log_sld, cutoff_rep_patient_log_sld;
vector[n_cutoff_total_forecast_visits] cutoff_forecast_mean_patient_log_sld, cutoff_forecast_patient_log_sld;

array[n_cutoff_visits] int<lower = CR, upper = PD + 1> cutoff_rep_recist = rep_array(PD + 1, n_cutoff_visits);
array[n_cutoff_total_forecast_visits] int<lower = CR, upper = PD> cutoff_forecast_recist;
  
// Endpoints (PFS, ORR, Median PFS, PFSn, ...) - sized for cutoff-observed patients
array[n_cutoff_observed_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, spop_target_obs_cens_pfs;
array[n_cutoff_observed_patients] int<lower = 0, upper = 1> sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored;

// Forecasting for right-censored patients at cutoff
array[n_cutoff_right_censored_patients] int<lower = 0> forecast_target_pfs;
array[n_cutoff_right_censored_patients] int<lower = 0, upper = 1> forecast_target_right_censored;

array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est;
array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est;

array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n,
                                                               sample_other_events_pfs_n, spop_other_events_pfs_n,
                                                               sample_pfs_n, spop_pfs_n;
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n,
                                                                   cond_sample_other_events_pfs_n, cond_spop_other_events_pfs_n,
                                                                   cond_sample_pfs_n, cond_spop_pfs_n;

array[n_cutoff_observed_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
array[n_cutoff_observed_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;

vector<lower = 0, upper = 1>[n_trials] sample_target_orr, spop_target_orr;
vector<lower = 0, upper = 1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group), cond_spop_target_orr = zeros_vector(n_cond_group);

// PFS quantiles - target events
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_trials, n_pfs_quantiles] int sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// PFS quantiles - other events
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_other_events_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_other_events_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_other_events_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_other_events_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_trials, n_pfs_quantiles] int sample_other_events_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_other_events_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_other_events_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_other_events_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// PFS quantiles - combined
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_trials, n_pfs_quantiles] int sample_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

vector<lower = 0>[n_trials] sample_target_median_pfs = zeros_vector(n_trials), spop_target_median_pfs = zeros_vector(n_trials);
vector<lower = 0>[n_cond_group] cond_sample_target_median_pfs = zeros_vector(n_cond_group), cond_spop_target_median_pfs = zeros_vector(n_cond_group);

array[n_trials] int sample_target_median_pfs_exceeds_max = zeros_int_array(n_trials), spop_target_median_pfs_exceeds_max = zeros_int_array(n_trials);
array[n_cond_group] int cond_sample_target_median_pfs_exceeds_max = zeros_int_array(n_cond_group), cond_spop_target_median_pfs_exceeds_max = zeros_int_array(n_cond_group);

profile("gen_quant") {
  // Subset states and patient-level parameters using fancy indexing
  // Note: cutoff_state_indices extracts ALL visits (including first) for cutoff-observed patients
  matrix[n_cutoff_visits, 2] cutoff_states = states[cutoff_state_indices, ];
  vector[n_cutoff_observed_patients] cutoff_patient_log_decrease_rate = patient_log_decrease_rate[cutoff_observed_patients];
  vector[n_cutoff_observed_patients] cutoff_patient_log_growth_rate = patient_log_growth_rate[cutoff_observed_patients];
  vector[n_cutoff_visits] cutoff_sum_tumor_size = sum_tumor_size[cutoff_state_indices];
  
  // Generate states for cutoff-observed patients (including forecasts for censored patients)
  (cutoff_forecast_patient_states, cutoff_rep_patient_log_sld, cutoff_rep_mean_patient_log_sld,
   cutoff_forecast_patient_log_sld, cutoff_forecast_mean_patient_log_sld,
   cutoff_obs_patient_process_noise, cutoff_forecast_patient_process_noise) = 
    generate_all_patients_states_with_means_rng(
      cutoff_states,
      cutoff_patient_visit_pos,
      cutoff_patient_visit_m1_pos,
      cutoff_forecast_visits_pos,
      cutoff_patient_last_obs_visit,
      max_all_t,  // Forecast up to max time
      cutoff_t_patient_visits,
      cutoff_patient_log_decrease_rate,
      cutoff_patient_log_growth_rate,
      cutoff_sum_tumor_size,
      0.0001, 0.0001, // forecast_growth_lag, forecast_growth_transition
      measure_sd,
      independ_long_process_noise,
      independ_cross_process_noise,
      pop_process_sd,
      L_process_corr,
      log_pop_tumor_gp_rho,
      delta,
      cutoff_n_patient_screening_visits
    );
  
  // Calculate RECIST classifications for cutoff-observed patients (including forecasts)
  (cutoff_rep_recist, cutoff_forecast_recist) = calculate_all_patients_recist(
    cutoff_rep_mean_patient_log_sld,
    cutoff_forecast_mean_patient_log_sld,
    cutoff_patient_visit_pos,
    cutoff_forecast_visits_pos,
    cutoff_n_patient_screening_visits
  );
  
  // Calculate patient-level endpoints for cutoff-observed patients (including forecasts)
  array[n_cutoff_observed_patients] int sample_other_events_pfs, sample_other_events_right_censored;
  array[n_cutoff_observed_patients] int spop_other_events_pfs, spop_other_events_right_censored;
  array[n_cutoff_observed_patients] int sample_pfs, sample_right_censored;
  array[n_cutoff_observed_patients] int spop_pfs, spop_right_censored;
  
  (sample_target_pfs, sample_target_right_censored,
   spop_target_pfs, spop_target_right_censored,
   spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
   sample_other_events_pfs, sample_other_events_right_censored,
   spop_other_events_pfs, spop_other_events_right_censored,
   sample_pfs, sample_right_censored,
   spop_pfs, spop_right_censored,
   sample_target_confirmed_response, sample_target_unconfirmed_response,
   spop_target_confirmed_response, spop_target_unconfirmed_response,
   forecast_target_pfs, forecast_target_right_censored) =
    calculate_all_patients_endpoints_rng(
      cutoff_recist,
      cutoff_rep_recist,
      cutoff_forecast_recist,
      log_cond_prob_surv[cutoff_observed_patients, ],  // Subset to cutoff-observed patients
      cutoff_pfs,
      cutoff_interval_censored,
      cutoff_right_censored,
      cutoff_other_events_right_censored,
      cutoff_patient_visit_pos,
      cutoff_forecast_visits_pos,
      cutoff_patient_last_obs_visit,
      max_all_t,
      cutoff_t_patient_visits,
      max_all_t,
      cutoff_n_patient_screening_visits
    );
    
  // Aggregate to trial-level metrics (using cutoff-observed patients only)
  array[n_trials] vector[max_all_t + 1] sample_other_events_km_est, spop_other_events_km_est;
  array[n_trials] vector[max_all_t + 1] sample_km_est, spop_km_est;
  
  (sample_target_orr, spop_target_orr,
   sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
   sample_other_events_km_est, spop_other_events_km_est,
   sample_km_est, spop_km_est,
   sample_target_quant_pfs, spop_target_quant_pfs,
   sample_target_quant_pfs_exceeds_max, spop_target_quant_pfs_exceeds_max,
   sample_other_events_quant_pfs, spop_other_events_quant_pfs,
   sample_other_events_quant_pfs_exceeds_max, spop_other_events_quant_pfs_exceeds_max,
   sample_quant_pfs, spop_quant_pfs,
   sample_quant_pfs_exceeds_max, spop_quant_pfs_exceeds_max,
   sample_target_pfs_n, spop_target_pfs_n,
   sample_other_events_pfs_n, spop_other_events_pfs_n,
   sample_pfs_n, spop_pfs_n) =
    aggregate_trial_metrics(
      sample_target_confirmed_response,
      spop_target_confirmed_response,
      sample_target_pfs,
      sample_target_right_censored,
      spop_target_pfs,
      spop_target_right_censored,
      spop_target_obs_cens_pfs,
      spop_target_obs_cens_right_censored,
      sample_other_events_pfs,
      sample_other_events_right_censored,
      spop_other_events_pfs,
      spop_other_events_right_censored,
      sample_pfs,
      sample_right_censored,
      spop_pfs,
      spop_right_censored,
      cutoff_trial_patient_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );
  
  // Aggregate to conditional group-level metrics (using cutoff-observed patients only)
  array[n_cond_group] vector[max_all_t + 1] cond_sample_other_events_km_est, cond_spop_other_events_km_est;
  array[n_cond_group] vector[max_all_t + 1] cond_sample_km_est, cond_spop_km_est;
  
  (cond_sample_target_orr, cond_spop_target_orr,
   cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
   cond_sample_other_events_km_est, cond_spop_other_events_km_est,
   cond_sample_km_est, cond_spop_km_est,
   cond_sample_target_quant_pfs, cond_spop_target_quant_pfs,
   cond_sample_target_quant_pfs_exceeds_max, cond_spop_target_quant_pfs_exceeds_max,
   cond_sample_other_events_quant_pfs, cond_spop_other_events_quant_pfs,
   cond_sample_other_events_quant_pfs_exceeds_max, cond_spop_other_events_quant_pfs_exceeds_max,
   cond_sample_quant_pfs, cond_spop_quant_pfs,
   cond_sample_quant_pfs_exceeds_max, cond_spop_quant_pfs_exceeds_max,
   cond_sample_target_pfs_n, cond_spop_target_pfs_n,
   cond_sample_other_events_pfs_n, cond_spop_other_events_pfs_n,
   cond_sample_pfs_n, cond_spop_pfs_n) =
    aggregate_conditional_group_metrics(
      sample_target_confirmed_response,
      spop_target_confirmed_response,
      sample_target_pfs,
      sample_target_right_censored,
      spop_target_pfs,
      spop_target_right_censored,
      spop_target_obs_cens_pfs,
      spop_target_obs_cens_right_censored,
      sample_other_events_pfs,
      sample_other_events_right_censored,
      spop_other_events_pfs,
      spop_other_events_right_censored,
      sample_pfs,
      sample_right_censored,
      spop_pfs,
      spop_right_censored,
      cutoff_cond_group,
      cutoff_cond_group_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );
}
