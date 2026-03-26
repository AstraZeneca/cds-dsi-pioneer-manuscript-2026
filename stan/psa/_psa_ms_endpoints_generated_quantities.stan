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
array[n_hmc_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs,
                                 spop_target_obs_cens_pfs, sample_pfs, spop_pfs;
array[n_hmc_patients] int<lower = 0, upper = 1>
  sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
  sample_ms_right_censored, spop_ms_right_censored,
  sample_right_censored, spop_right_censored;

// OS endpoints
array[n_hmc_patients] int<lower = 0> sample_os, spop_os;
array[n_hmc_patients] int<lower = 0, upper = 1> sample_os_censored, spop_os_censored;

// Forecast for right-censored patients
array[sum(psa_right_censored)] int<lower = 0> forecast_target_pfs;
array[sum(psa_right_censored)] int<lower = 0, upper = 1> forecast_target_right_censored;

// Trial-level KM curves
array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_target_km_est, spop_target_km_est,
                                                            spop_target_obs_cens_km_est,
                                                            sample_ms_pfs_km_est, spop_ms_pfs_km_est,
                                                            sample_pfs_km_est, spop_pfs_km_est;
array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_target_km_est, cond_spop_target_km_est,
                                                                cond_spop_target_obs_cens_km_est,
                                                                cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
                                                                cond_sample_pfs_km_est, cond_spop_pfs_km_est;

// PFS at fixed timepoints
array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n,
                                                               sample_ms_pfs_n, spop_ms_pfs_n,
                                                               sample_pfs_n, spop_pfs_n;
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n,
                                                                   cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
                                                                   cond_sample_pfs_n, cond_spop_pfs_n;

// Response (PSA50 = unconfirmed, undetectable = confirmed)
array[n_hmc_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
array[n_hmc_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;
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
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_trials, n_pfs_quantiles] int sample_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// PFS quantiles — combined (PSA + multistate)
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_trials, n_pfs_quantiles] int sample_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

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

// ── Radiographic PFS (rPFS = 0→1 radiographic OR 0→2 direct death) ──────
// Meaningful only for full pioneer model where 0→1 = radiographic hazard.
// Does not include PSA-PD events (use pfs_km_est for those).
array[n_hmc_patients] int<lower=0> sample_rpfs, spop_rpfs;
array[n_hmc_patients] int<lower=0, upper=1> sample_rpfs_censored, spop_rpfs_censored;

// Trial-level rPFS
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_rpfs_km_est, spop_rpfs_km_est;
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_rpfs_km_est, cond_spop_rpfs_km_est;

array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_rpfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_rpfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int
  sample_rpfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_rpfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);

array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_rpfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_rpfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int
  cond_sample_rpfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_rpfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints] sample_rpfs_n, spop_rpfs_n;
array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints] cond_sample_rpfs_n, cond_spop_rpfs_n;

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
    hmc_patient_idx,
    pcwg3_category,
    rep_pcwg3,
    forecast_pcwg3,
    log_cond_surv_01,
    enable_ms_02,
    enable_ms_03,
    enable_ms_32,
    enable_ms_12,
    ms_time_scale_12,
    log_cond_surv_02,
    log_cond_surv_03,
    log_cond_surv_32,
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
    ms_time_03,
    ms_time_32,
    ms_censored_32,
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
   sample_ms_pfs_km_est, spop_ms_pfs_km_est,
   sample_pfs_km_est, spop_pfs_km_est,
   sample_target_quant_pfs, spop_target_quant_pfs,
   sample_target_quant_pfs_exceeds_max, spop_target_quant_pfs_exceeds_max,
   sample_ms_pfs_quant, spop_ms_pfs_quant,
   sample_ms_pfs_quant_exceeds_max, spop_ms_pfs_quant_exceeds_max,
   sample_pfs_quant, spop_pfs_quant,
   sample_pfs_quant_exceeds_max, spop_pfs_quant_exceeds_max,
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
      hmc_trial_patient_pos,
      max_all_t,
      pfs_quantiles,
      pfs_timepoints
    );

  // Aggregate to conditional group-level metrics
  (cond_sample_target_orr, cond_spop_target_orr,
   cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
   cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
   cond_sample_pfs_km_est, cond_spop_pfs_km_est,
   cond_sample_target_quant_pfs, cond_spop_target_quant_pfs,
   cond_sample_target_quant_pfs_exceeds_max, cond_spop_target_quant_pfs_exceeds_max,
   cond_sample_ms_pfs_quant, cond_spop_ms_pfs_quant,
   cond_sample_ms_pfs_quant_exceeds_max, cond_spop_ms_pfs_quant_exceeds_max,
   cond_sample_pfs_quant, cond_spop_pfs_quant,
   cond_sample_pfs_quant_exceeds_max, cond_spop_pfs_quant_exceeds_max,
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

// ── Radiographic PFS patient-level derivation ────────────────────────────
// If ms event (radiographic progression): rPFS event at ms_pfs time.
// If died without prior progression (os event, ms censored): rPFS event at os time.
// Otherwise: censored at ms_pfs censoring time.
for (i in 1:n_hmc_patients) {
  if (!spop_ms_right_censored[i]) {
    spop_rpfs[i] = spop_ms_pfs[i]; spop_rpfs_censored[i] = 0;
  } else if (!spop_os_censored[i]) {
    spop_rpfs[i] = spop_os[i]; spop_rpfs_censored[i] = 0;
  } else {
    spop_rpfs[i] = spop_ms_pfs[i]; spop_rpfs_censored[i] = 1;
  }

  if (!sample_ms_right_censored[i]) {
    sample_rpfs[i] = sample_ms_pfs[i]; sample_rpfs_censored[i] = 0;
  } else if (!sample_os_censored[i]) {
    sample_rpfs[i] = sample_os[i]; sample_rpfs_censored[i] = 0;
  } else {
    sample_rpfs[i] = sample_ms_pfs[i]; sample_rpfs_censored[i] = 1;
  }
}

// ── Radiographic PFS trial-level aggregation ─────────────────────────────
for (s in 1:n_trials) {
  if (get_pos_size(hmc_trial_patient_pos, s) > 0) {
    sample_rpfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(sample_rpfs, hmc_trial_patient_pos, s),
      get_int_sub_array(sample_rpfs_censored, hmc_trial_patient_pos, s),
      max_all_t, 0).1;
    spop_rpfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(spop_rpfs, hmc_trial_patient_pos, s),
      get_int_sub_array(spop_rpfs_censored, hmc_trial_patient_pos, s),
      max_all_t, 0).1;

    (sample_rpfs_quant[s], sample_rpfs_quant_exceeds_max[s]) =
      km_quantiles(sample_rpfs_km_est[s], pfs_quantiles);
    (spop_rpfs_quant[s], spop_rpfs_quant_exceeds_max[s]) =
      km_quantiles(spop_rpfs_km_est[s], pfs_quantiles);

    for (n in 1:n_pfs_timepoints) {
      sample_rpfs_n[s, n] = calc_km_pfs_n(sample_rpfs_km_est[s], months_to_weeks(pfs_timepoints[n]));
      spop_rpfs_n[s, n]   = calc_km_pfs_n(spop_rpfs_km_est[s],   months_to_weeks(pfs_timepoints[n]));
    }
  } else {
    sample_rpfs_km_est[s] = zeros_vector(max_all_t + 1);
    spop_rpfs_km_est[s]   = zeros_vector(max_all_t + 1);
    sample_rpfs_n[s]      = zeros_vector(n_pfs_timepoints);
    spop_rpfs_n[s]        = zeros_vector(n_pfs_timepoints);
  }
}

// ── Radiographic PFS conditional group aggregation ───────────────────────
for (c in 1:n_cond_group) {
  int curr_group_size = get_pos_size(cond_group_pos, c);
  if (curr_group_size > 0) {
    array[curr_group_size] int curr_group_patients =
      get_int_sub_array(cond_group, cond_group_pos, c);

    cond_sample_rpfs_km_est[c] = estimate_kaplan_meier(
      sample_rpfs[curr_group_patients], sample_rpfs_censored[curr_group_patients],
      max_all_t, 0).1;
    cond_spop_rpfs_km_est[c] = estimate_kaplan_meier(
      spop_rpfs[curr_group_patients], spop_rpfs_censored[curr_group_patients],
      max_all_t, 0).1;

    (cond_sample_rpfs_quant[c], cond_sample_rpfs_quant_exceeds_max[c]) =
      km_quantiles(cond_sample_rpfs_km_est[c], pfs_quantiles);
    (cond_spop_rpfs_quant[c], cond_spop_rpfs_quant_exceeds_max[c]) =
      km_quantiles(cond_spop_rpfs_km_est[c], pfs_quantiles);

    for (n in 1:n_pfs_timepoints) {
      cond_sample_rpfs_n[c, n] = calc_km_pfs_n(cond_sample_rpfs_km_est[c], months_to_weeks(pfs_timepoints[n]));
      cond_spop_rpfs_n[c, n]   = calc_km_pfs_n(cond_spop_rpfs_km_est[c],   months_to_weeks(pfs_timepoints[n]));
    }
  } else {
    cond_sample_rpfs_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_rpfs_km_est[c]   = zeros_vector(max_all_t + 1);
    cond_sample_rpfs_n[c]      = zeros_vector(n_pfs_timepoints);
    cond_spop_rpfs_n[c]        = zeros_vector(n_pfs_timepoints);
  }
}
