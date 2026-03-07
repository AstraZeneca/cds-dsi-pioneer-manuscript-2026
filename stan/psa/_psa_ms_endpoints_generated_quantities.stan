// ============================================================================
// PSA MULTISTATE ENDPOINT GENERATED QUANTITIES
// ============================================================================
// Computes PFS, KM curves, quantiles, and OS endpoints using PCWG3 response
// categories and multistate hazards.  Mirrors _tumor_endpoints_generated_quantities.stan
// for the PSA biomarker (pioneer model).
//
// Prerequisites (already computed earlier in the GQ block):
//   rep_pcwg3, forecast_pcwg3   — from _psa_endpoints_generated_quantities.stan
//   forecast_patient_states     — from modules/state_space/generated_quantities.stan
//   log_cond_surv_*             — from modules/multistate/transformed_parameters.stan
//
// Variable mapping vs. tumor model:
//   recist / rep_recist / forecast_recist  →  pcwg3_category / rep_pcwg3 / forecast_pcwg3
//   pfs / interval_censored / right_censored  →  psa_pfs / psa_interval_censored / psa_right_censored
//   target_pfs / target_right_censored  →  psa_pfs / psa_right_censored
//     (no "target vs. non-target" distinction in the PSA context)

// Patient-level PFS arrays
array[n_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs,
                                 spop_target_obs_cens_pfs, sample_pfs, spop_pfs;
array[n_patients] int<lower = 0, upper = 1>
  sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
  sample_ms_right_censored, spop_ms_right_censored,
  sample_right_censored, spop_right_censored;

// OS endpoints
array[n_patients] int<lower = 0> sample_os, spop_os;
array[n_patients] int<lower = 0, upper = 1> sample_os_censored, spop_os_censored;

// Forecast for right-censored patients
array[sum(psa_right_censored)] int<lower = 0> forecast_target_pfs;
array[sum(psa_right_censored)] int<lower = 0, upper = 1> forecast_target_right_censored;

// Trial-level KM curves
array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_target_km_est, spop_target_km_est,
                                                            spop_target_obs_cens_km_est,
                                                            sample_ms_km_est, spop_ms_km_est,
                                                            sample_km_est, spop_km_est;
array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_target_km_est, cond_spop_target_km_est,
                                                                cond_spop_target_obs_cens_km_est,
                                                                cond_sample_ms_km_est, cond_spop_ms_km_est,
                                                                cond_sample_km_est, cond_spop_km_est;

// PFS at fixed timepoints
array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n,
                                                               sample_ms_pfs_n, spop_ms_pfs_n,
                                                               sample_pfs_n, spop_pfs_n;
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n,
                                                                   cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
                                                                   cond_sample_pfs_n, cond_spop_pfs_n;

// Response (PSA50 = unconfirmed, undetectable = confirmed)
array[n_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
array[n_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;
vector<lower = 0, upper = 1>[n_trials] sample_target_orr, spop_target_orr;
vector<lower = 0, upper = 1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group),
                                            cond_spop_target_orr = zeros_vector(n_cond_group);

// PFS quantiles — target (PSA-based) events
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_trials, n_pfs_quantiles] int sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// PFS quantiles — multistate events
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_trials, n_pfs_quantiles] int sample_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// PFS quantiles — combined (PSA + multistate)
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_trials, n_pfs_quantiles] int sample_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// OS KM curves
array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_os_km_est, spop_os_km_est;
array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_os_km_est, cond_spop_os_km_est;

// OS quantiles
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// OS at fixed timepoints
array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_os_n, spop_os_n;
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_os_n, cond_spop_os_n;

profile("psa_ms_endpoints") {
  // Calculate patient-level endpoints using PCWG3 categories as the biomarker signal.
  // pcwg3_category uses the same 1–4 progression code as RECIST (4 = PSA-PD),
  // so calculate_all_patients_endpoints_rng is fully compatible.
  (sample_target_pfs, sample_target_right_censored,
   spop_target_pfs, spop_target_right_censored,
   spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
   sample_ms_pfs, sample_ms_right_censored,
   spop_ms_pfs, spop_ms_right_censored,
   sample_pfs, sample_right_censored,
   spop_pfs, spop_right_censored,
   sample_target_confirmed_response, sample_target_unconfirmed_response,
   spop_target_confirmed_response, spop_target_unconfirmed_response,
   forecast_target_pfs, forecast_target_right_censored,
   sample_os, sample_os_censored,
   spop_os, spop_os_censored
  ) = calculate_all_patients_endpoints_rng(
    pcwg3_category,
    rep_pcwg3,
    forecast_pcwg3,
    log_cond_surv_01,
    enable_ms_02,
    enable_ms_12,
    ms_time_scale_12,
    log_cond_surv_02,
    log_cond_surv_12_s,
    log_cond_surv_12_t,
    psa_pfs,
    psa_interval_censored,
    psa_right_censored,
    psa_pfs,              // target_pfs = psa_pfs (no target vs. non-target distinction)
    psa_right_censored,   // target_right_censored = psa_right_censored
    ms_censored_01,       // ms_right_censored: censoring status for 0→1 transition
    ms_final_state,
    ms_time_02,
    ms_censored_02,
    ms_censored_12,
    ms_time_12,
    patient_visit_pos,
    forecast_visits_pos,
    patient_last_obs_visit,
    last_predict_visit,
    t_patient_visits,
    max_all_t,
    n_patient_screening_visits
  );

  // Aggregate to trial-level metrics
  (sample_target_orr, spop_target_orr,
   sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
   sample_ms_km_est, spop_ms_km_est,
   sample_km_est, spop_km_est,
   sample_target_quant_pfs, spop_target_quant_pfs,
   sample_target_quant_pfs_exceeds_max, spop_target_quant_pfs_exceeds_max,
   sample_ms_quant_pfs, spop_ms_quant_pfs,
   sample_ms_quant_pfs_exceeds_max, spop_ms_quant_pfs_exceeds_max,
   sample_quant_pfs, spop_quant_pfs,
   sample_quant_pfs_exceeds_max, spop_quant_pfs_exceeds_max,
   sample_target_pfs_n, spop_target_pfs_n,
   sample_ms_pfs_n, spop_ms_pfs_n,
   sample_pfs_n, spop_pfs_n,
   sample_os_km_est, spop_os_km_est,
   sample_os_quant, spop_os_quant,
   sample_os_quant_exceeds_max, spop_os_quant_exceeds_max,
   sample_os_n, spop_os_n) =
    aggregate_trial_metrics(
      sample_target_confirmed_response,
      spop_target_confirmed_response,
      sample_target_pfs,
      sample_target_right_censored,
      spop_target_pfs,
      spop_target_right_censored,
      spop_target_obs_cens_pfs,
      spop_target_obs_cens_right_censored,
      sample_ms_pfs,
      sample_ms_right_censored,
      spop_ms_pfs,
      spop_ms_right_censored,
      sample_pfs,
      sample_right_censored,
      spop_pfs,
      spop_right_censored,
      sample_os,
      sample_os_censored,
      spop_os,
      spop_os_censored,
      trial_patient_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );

  // Aggregate to conditional group-level metrics
  (cond_sample_target_orr, cond_spop_target_orr,
   cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
   cond_sample_ms_km_est, cond_spop_ms_km_est,
   cond_sample_km_est, cond_spop_km_est,
   cond_sample_target_quant_pfs, cond_spop_target_quant_pfs,
   cond_sample_target_quant_pfs_exceeds_max, cond_spop_target_quant_pfs_exceeds_max,
   cond_sample_ms_quant_pfs, cond_spop_ms_quant_pfs,
   cond_sample_ms_quant_pfs_exceeds_max, cond_spop_ms_quant_pfs_exceeds_max,
   cond_sample_quant_pfs, cond_spop_quant_pfs,
   cond_sample_quant_pfs_exceeds_max, cond_spop_quant_pfs_exceeds_max,
   cond_sample_target_pfs_n, cond_spop_target_pfs_n,
   cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
   cond_sample_pfs_n, cond_spop_pfs_n,
   cond_sample_os_km_est, cond_spop_os_km_est,
   cond_sample_os_quant, cond_spop_os_quant,
   cond_sample_os_quant_exceeds_max, cond_spop_os_quant_exceeds_max,
   cond_sample_os_n, cond_spop_os_n) =
    aggregate_conditional_group_metrics(
      sample_target_confirmed_response,
      spop_target_confirmed_response,
      sample_target_pfs,
      sample_target_right_censored,
      spop_target_pfs,
      spop_target_right_censored,
      spop_target_obs_cens_pfs,
      spop_target_obs_cens_right_censored,
      sample_ms_pfs,
      sample_ms_right_censored,
      spop_ms_pfs,
      spop_ms_right_censored,
      sample_pfs,
      sample_right_censored,
      spop_pfs,
      spop_right_censored,
      sample_os,
      sample_os_censored,
      spop_os,
      spop_os_censored,
      cond_group,
      cond_group_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );
}
