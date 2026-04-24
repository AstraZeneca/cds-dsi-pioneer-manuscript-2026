// ============================================================================
// MS-Standalone LFO: Posterior-predictive KM endpoints from re-censored events
// ============================================================================
// Mirrors the regular ms-standalone GQ (tumor/_ms_standalone_generated_quantities.stan)
// — calculate_all_patients_endpoints_rng + aggregate_trial_metrics +
// aggregate_conditional_group_metrics — but operating on the COMPACT
// cutoff-observed cohort with re-censored (lfo_ms_*) event arrays.
//
// Trial-level KMs are expanded back to all n_patients (non-cutoff-observed
// treated as censored at t=0) so `sample_os_km_est[n_trials]` matches the
// regular model's convention. Cond_group-level KMs stay on the compact
// cohort.
//
// Requires in scope (from transformed_data / other GQ includes):
//   n_cutoff_observed_patients, cutoff_observed_patients[],
//   patient_to_cutoff_idx[], cutoff_trial_patient_pos[],
//   cutoff_cond_group[], cutoff_cond_group_pos[]
//   lfo_ms_final_state, lfo_ms_time_01/02/12/03/32,
//   lfo_ms_censored_01, lfo_ms_prog_deterministic
//   ms_os_event_12  (un-recensored; gating by lfo_ms_censored_12 below)
//   log_cond_surv_01/02/03/12_s/12_t/32
//   enable_ms_02/03/32/12, ms_time_scale_12
//   patient_visit_pos, forecast_visits_pos, forecast_obs_visits_pos,
//   patient_last_obs_visit, last_predict_visit, t_patient_visits,
//   n_patient_screening_visits, n_pfs_quantiles, pfs_quantiles, pfs_timepoints,
//   forecast_observation_interval

// --- Re-derive 02/12/32 censoring indicators on the re-censored arrays ------
array[n_patients] int lfo_ms_censored_02;
array[n_patients] int lfo_ms_censored_12;
array[n_patients] int lfo_ms_censored_32;
(lfo_ms_censored_02, lfo_ms_censored_12, lfo_ms_censored_32) =
  derive_ms_censoring_indicators(
    lfo_ms_final_state, lfo_ms_time_01, lfo_ms_time_32
  );

// --- Patient-level endpoint arrays (compact) --------------------------------
array[n_cutoff_observed_patients] int<lower=0>
  sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs,
  spop_target_obs_cens_pfs, sample_pfs, spop_pfs;
array[n_cutoff_observed_patients] int<lower=0, upper=1>
  sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
  sample_ms_right_censored, spop_ms_right_censored,
  sample_right_censored, spop_right_censored;
array[n_cutoff_observed_patients] int<lower=0> sample_os, spop_os;
array[n_cutoff_observed_patients] int<lower=0, upper=1>
  sample_os_censored, spop_os_censored,
  spop_is_dropout, sample_is_dropout,
  sample_target_confirmed_response, spop_target_confirmed_response,
  sample_target_unconfirmed_response, spop_target_unconfirmed_response;

// --- Trial-level KMs (expanded back to all n_patients) ----------------------
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1]
  sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
  sample_ms_pfs_km_est, spop_ms_pfs_km_est,
  sample_pfs_km_est,    spop_pfs_km_est,
  sample_os_km_est,     spop_os_km_est;
array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints]
  sample_target_pfs_n, spop_target_pfs_n,
  sample_ms_pfs_n,     spop_ms_pfs_n,
  sample_pfs_n,        spop_pfs_n,
  sample_os_n,         spop_os_n;
array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_target_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_target_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  sample_ms_pfs_quant     = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_ms_pfs_quant       = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  sample_pfs_quant        = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_pfs_quant          = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  sample_os_quant         = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_os_quant           = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int
  sample_target_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_target_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  sample_ms_pfs_quant_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_ms_pfs_quant_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  sample_pfs_quant_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_pfs_quant_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  sample_os_quant_exceeds_max         = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_os_quant_exceeds_max           = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
vector<lower=0, upper=1>[n_trials] sample_target_orr, spop_target_orr;

// --- Cond-group KMs (compact cohort only) -----------------------------------
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
  cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
  cond_sample_pfs_km_est,    cond_spop_pfs_km_est,
  cond_sample_os_km_est,     cond_spop_os_km_est;
array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints]
  cond_sample_target_pfs_n, cond_spop_target_pfs_n,
  cond_sample_ms_pfs_n,     cond_spop_ms_pfs_n,
  cond_sample_pfs_n,        cond_spop_pfs_n,
  cond_sample_os_n,         cond_spop_os_n;
array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_target_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_target_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_sample_ms_pfs_quant     = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_pfs_quant       = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_sample_pfs_quant        = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_pfs_quant          = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_sample_os_quant         = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_os_quant           = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int
  cond_sample_target_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_target_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_sample_ms_pfs_quant_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_pfs_quant_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_sample_pfs_quant_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_pfs_quant_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_sample_os_quant_exceeds_max         = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_os_quant_exceeds_max           = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
vector<lower=0, upper=1>[n_cond_group]
  cond_sample_target_orr = zeros_vector(n_cond_group),
  cond_spop_target_orr   = zeros_vector(n_cond_group);

// --- Compute patient-level endpoints on the compact cohort ------------------
// Dummy biomarker categories: standalone model has no biomarker dynamics, so
// target ORR is suppressed and km_est == ms_km_est.
{
  array[n_total_visits] int obs_biomarker_cat = rep_array(0, n_total_visits);
  array[n_total_visits] int rep_biomarker_cat = rep_array(0, n_total_visits);
  array[n_total_forecast_obs_visits] int forecast_obs_biomarker_cat =
    rep_array(0, n_total_forecast_obs_visits);
  // Burden/target surrogate: PFS = progression (0→1) OR death (0→2);
  // dropout (0→3) is censored. Same pattern as the regular standalone GQ.
  array[n_patients] int burden_pfs;
  array[n_patients] int burden_right_censored;
  for (i in 1:n_patients) {
    if (lfo_ms_censored_01[i] == 0) {
      burden_pfs[i] = lfo_ms_time_01[i];
      burden_right_censored[i] = 0;
    } else if (lfo_ms_censored_02[i] == 0) {
      burden_pfs[i] = lfo_ms_time_02[i];
      burden_right_censored[i] = 0;
    } else {
      burden_pfs[i] = min(lfo_ms_time_01[i], lfo_ms_time_02[i]);
      burden_right_censored[i] = 1;
    }
  }
  array[n_patients] int burden_target_pfs            = rep_array(max_all_t + 1, n_patients);
  array[n_patients] int burden_target_right_censored = rep_array(1, n_patients);

  array[n_cutoff_observed_patients] int forecast_target_pfs_local;
  array[n_cutoff_observed_patients] int forecast_target_right_censored_local;

  (sample_target_pfs, sample_target_right_censored,
   spop_target_pfs, spop_target_right_censored,
   spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
   sample_ms_pfs, sample_ms_right_censored,
   spop_ms_pfs, spop_ms_right_censored,
   sample_pfs, sample_right_censored,
   spop_pfs, spop_right_censored,
   sample_target_confirmed_response, sample_target_unconfirmed_response,
   spop_target_confirmed_response, spop_target_unconfirmed_response,
   forecast_target_pfs_local, forecast_target_right_censored_local,
   sample_os, sample_os_censored,
   spop_os, spop_os_censored,
   spop_is_dropout, sample_is_dropout
  ) = calculate_all_patients_endpoints_rng(
    cutoff_observed_patients,
    obs_biomarker_cat,
    rep_biomarker_cat,
    forecast_obs_biomarker_cat,
    forecast_obs_visits_pos,
    forecast_observation_interval,
    log_cond_surv_01,
    enable_ms_02, enable_ms_03, enable_ms_32, enable_ms_12,
    ms_time_scale_12,
    log_cond_surv_02, log_cond_surv_03, log_cond_surv_32,
    log_cond_surv_12_s, log_cond_surv_12_t,
    burden_pfs,
    burden_right_censored,
    burden_target_pfs,
    burden_target_right_censored,
    lfo_ms_censored_01,
    lfo_ms_final_state,
    lfo_ms_time_02,
    lfo_ms_censored_02,
    lfo_ms_censored_12,
    lfo_ms_time_12,
    ms_os_event_12,
    lfo_ms_time_03,
    lfo_ms_time_32,
    lfo_ms_censored_32,
    patient_visit_pos,
    forecast_visits_pos,
    patient_last_obs_visit,
    last_predict_visit,
    t_patient_visits,
    max_all_t,
    n_patient_screening_visits,
    0,                // enable_ms_visit_gated_01 = 0 (no PSA gating in ms_standalone)
    0.0,              // tv_coef_01_val
    zeros_vector(0),  // forecast_obs_log_psa
    0.0,              // median_log_psa_obs
    0.0               // iqr_log_psa_obs
  );
}

// --- Expand per-patient arrays from compact → n_patients for trial KM -------
// Non-cutoff-observed patients are treated as censored at t=0 so they
// contribute a trivial row to the KM (matches the regular model's convention
// where every enrolled patient participates in the denominator).
array[n_patients] int all_sample_target_pfs              = zeros_int_array(n_patients);
array[n_patients] int all_sample_target_right_censored   = rep_array(1, n_patients);
array[n_patients] int all_spop_target_pfs                = zeros_int_array(n_patients);
array[n_patients] int all_spop_target_right_censored     = rep_array(1, n_patients);
array[n_patients] int all_spop_target_obs_cens_pfs       = zeros_int_array(n_patients);
array[n_patients] int all_spop_target_obs_cens_right_censored = rep_array(1, n_patients);
array[n_patients] int all_sample_ms_pfs                  = zeros_int_array(n_patients);
array[n_patients] int all_sample_ms_right_censored       = rep_array(1, n_patients);
array[n_patients] int all_spop_ms_pfs                    = zeros_int_array(n_patients);
array[n_patients] int all_spop_ms_right_censored         = rep_array(1, n_patients);
array[n_patients] int all_sample_pfs                     = zeros_int_array(n_patients);
array[n_patients] int all_sample_right_censored          = rep_array(1, n_patients);
array[n_patients] int all_spop_pfs                       = zeros_int_array(n_patients);
array[n_patients] int all_spop_right_censored            = rep_array(1, n_patients);
array[n_patients] int all_sample_target_confirmed_response = zeros_int_array(n_patients);
array[n_patients] int all_spop_target_confirmed_response   = zeros_int_array(n_patients);
array[n_patients] int all_sample_os                      = zeros_int_array(n_patients);
array[n_patients] int all_sample_os_censored             = rep_array(1, n_patients);
array[n_patients] int all_spop_os                        = zeros_int_array(n_patients);
array[n_patients] int all_spop_os_censored               = rep_array(1, n_patients);
for (obs_idx in 1:n_cutoff_observed_patients) {
  int i = cutoff_observed_patients[obs_idx];
  all_sample_target_pfs[i]                       = sample_target_pfs[obs_idx];
  all_sample_target_right_censored[i]            = sample_target_right_censored[obs_idx];
  all_spop_target_pfs[i]                         = spop_target_pfs[obs_idx];
  all_spop_target_right_censored[i]              = spop_target_right_censored[obs_idx];
  all_spop_target_obs_cens_pfs[i]                = spop_target_obs_cens_pfs[obs_idx];
  all_spop_target_obs_cens_right_censored[i]     = spop_target_obs_cens_right_censored[obs_idx];
  all_sample_ms_pfs[i]                           = sample_ms_pfs[obs_idx];
  all_sample_ms_right_censored[i]                = sample_ms_right_censored[obs_idx];
  all_spop_ms_pfs[i]                             = spop_ms_pfs[obs_idx];
  all_spop_ms_right_censored[i]                  = spop_ms_right_censored[obs_idx];
  all_sample_pfs[i]                              = sample_pfs[obs_idx];
  all_sample_right_censored[i]                   = sample_right_censored[obs_idx];
  all_spop_pfs[i]                                = spop_pfs[obs_idx];
  all_spop_right_censored[i]                     = spop_right_censored[obs_idx];
  all_sample_target_confirmed_response[i]        = sample_target_confirmed_response[obs_idx];
  all_spop_target_confirmed_response[i]          = spop_target_confirmed_response[obs_idx];
  all_sample_os[i]                               = sample_os[obs_idx];
  all_sample_os_censored[i]                      = sample_os_censored[obs_idx];
  all_spop_os[i]                                 = spop_os[obs_idx];
  all_spop_os_censored[i]                        = spop_os_censored[obs_idx];
}

// --- Trial-level aggregation (all enrolled patients, standard convention) ---
(sample_target_orr, spop_target_orr,
 sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
 sample_ms_pfs_km_est, spop_ms_pfs_km_est,
 sample_pfs_km_est,    spop_pfs_km_est,
 sample_target_pfs_quant, spop_target_pfs_quant,
 sample_target_pfs_quant_exceeds_max, spop_target_pfs_quant_exceeds_max,
 sample_ms_pfs_quant, spop_ms_pfs_quant,
 sample_ms_pfs_quant_exceeds_max, spop_ms_pfs_quant_exceeds_max,
 sample_pfs_quant, spop_pfs_quant,
 sample_pfs_quant_exceeds_max, spop_pfs_quant_exceeds_max,
 sample_target_pfs_n, spop_target_pfs_n,
 sample_ms_pfs_n,     spop_ms_pfs_n,
 sample_pfs_n,        spop_pfs_n,
 sample_os_km_est,    spop_os_km_est,
 sample_os_quant,     spop_os_quant,
 sample_os_quant_exceeds_max, spop_os_quant_exceeds_max,
 sample_os_n,         spop_os_n) =
  aggregate_trial_metrics(
    all_sample_target_confirmed_response, all_spop_target_confirmed_response,
    all_sample_target_pfs, all_sample_target_right_censored,
    all_spop_target_pfs, all_spop_target_right_censored,
    all_spop_target_obs_cens_pfs, all_spop_target_obs_cens_right_censored,
    all_sample_ms_pfs, all_sample_ms_right_censored,
    all_spop_ms_pfs,   all_spop_ms_right_censored,
    all_sample_pfs,    all_sample_right_censored,
    all_spop_pfs,      all_spop_right_censored,
    all_sample_os,     all_sample_os_censored,
    all_spop_os,       all_spop_os_censored,
    trial_patient_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );

// --- Cond_group aggregation (compact cohort, already filtered) --------------
(cond_sample_target_orr, cond_spop_target_orr,
 cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
 cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
 cond_sample_pfs_km_est,    cond_spop_pfs_km_est,
 cond_sample_target_pfs_quant, cond_spop_target_pfs_quant,
 cond_sample_target_pfs_quant_exceeds_max, cond_spop_target_pfs_quant_exceeds_max,
 cond_sample_ms_pfs_quant, cond_spop_ms_pfs_quant,
 cond_sample_ms_pfs_quant_exceeds_max, cond_spop_ms_pfs_quant_exceeds_max,
 cond_sample_pfs_quant, cond_spop_pfs_quant,
 cond_sample_pfs_quant_exceeds_max, cond_spop_pfs_quant_exceeds_max,
 cond_sample_target_pfs_n, cond_spop_target_pfs_n,
 cond_sample_ms_pfs_n,     cond_spop_ms_pfs_n,
 cond_sample_pfs_n,        cond_spop_pfs_n,
 cond_sample_os_km_est,    cond_spop_os_km_est,
 cond_sample_os_quant,     cond_spop_os_quant,
 cond_sample_os_quant_exceeds_max, cond_spop_os_quant_exceeds_max,
 cond_sample_os_n,         cond_spop_os_n) =
  aggregate_conditional_group_metrics(
    sample_target_confirmed_response, spop_target_confirmed_response,
    sample_target_pfs, sample_target_right_censored,
    spop_target_pfs, spop_target_right_censored,
    spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
    sample_ms_pfs, sample_ms_right_censored,
    spop_ms_pfs,   spop_ms_right_censored,
    sample_pfs,    sample_right_censored,
    spop_pfs,      spop_right_censored,
    sample_os,     sample_os_censored,
    spop_os,       spop_os_censored,
    cutoff_cond_group,
    cutoff_cond_group_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );
