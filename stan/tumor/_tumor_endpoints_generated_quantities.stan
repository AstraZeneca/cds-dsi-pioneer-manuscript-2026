// Tumor-specific trajectory variables (states are in state_space module)

vector[sum(n_patient_visits)] rep_mean_patient_log_sld, rep_patient_log_sld;
vector[n_total_forecast_visits] forecast_mean_patient_log_sld, forecast_patient_log_sld;

array[sum(n_patient_visits)] int<lower = CR, upper = PD + 1> rep_recist = rep_array(PD + 1, sum(n_patient_visits));
array[n_total_forecast_visits] int<lower = CR, upper = PD> forecast_recist;
  
// Endpoints (PFS, ORR, Median PFS, PFSn, ...) ////////////

array[n_forecast_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs, spop_target_obs_cens_pfs,
                                 sample_pfs, spop_pfs;
array[n_forecast_patients] int<lower = 0, upper = 1>
  sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
  sample_ms_right_censored, spop_ms_right_censored,
  sample_right_censored, spop_right_censored;

// OS endpoints
array[n_forecast_patients] int<lower = 0> sample_os, spop_os;
array[n_forecast_patients] int<lower = 0, upper = 1> sample_os_censored, spop_os_censored;

// Dropout flags — 1 if patient exited via cause 3 in this draw (for CIF computation)
array[n_patients] int<lower = 0, upper = 1> spop_is_dropout, sample_is_dropout;

// Forecasting for right censored patients 
array[sum(target_right_censored)] int<lower = 0> forecast_target_pfs;
array[sum(target_right_censored)] int<lower = 0, upper = 1> forecast_target_right_censored; 

array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_target_km_est, spop_target_km_est, 
                                                            spop_target_obs_cens_km_est, 
                                                            sample_ms_km_est, spop_ms_km_est,
                                                            sample_km_est, spop_km_est; 

array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_target_km_est, cond_spop_target_km_est, 
                                                                cond_spop_target_obs_cens_km_est, 
                                                                cond_sample_ms_km_est, cond_spop_ms_km_est,
                                                                cond_sample_km_est, cond_spop_km_est;

array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n, 
                                                               sample_ms_pfs_n, spop_ms_pfs_n,
                                                               sample_pfs_n, spop_pfs_n; 
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n,
                                                                   cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
                                                                   cond_sample_pfs_n, cond_spop_pfs_n;

array[n_forecast_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
array[n_forecast_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;

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

// PFS quantiles - multistate events
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_trials, n_pfs_quantiles] int sample_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// PFS quantiles - combined
array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials] vector<lower = 0>[n_pfs_quantiles] spop_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_spop_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_trials, n_pfs_quantiles] int sample_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int spop_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int cond_sample_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int cond_spop_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// OS KM estimates
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

// OS at timepoints
array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_os_n, spop_os_n;
array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_os_n, cond_spop_os_n;

profile("gen_quant") {
  int right_censored_idx = 1;

  // Generate states for all patients
  if (enable_patient_process_noise_tr) {
    // Process noise ON: Use states_full_grid (dense grid computed in transformed_parameters)
    (forecast_patient_states, rep_patient_log_sld, rep_mean_patient_log_sld,
     forecast_patient_log_sld, forecast_mean_patient_log_sld) =
      generate_all_patients_states_with_means_rng(
        states_full_grid,
        forecast_patient_idx,
        patient_visit_pos,
        patient_visit_m1_pos,
        forecast_visits_pos,
        patient_last_obs_visit,
        last_predict_visit,
        t_patient_visits,
        t_patient_visit_idx,
        sum_tumor_size,
        measure_sd_sld,
        measure_nu,
        n_patient_screening_visits
      );
  } else {
    // Process noise OFF: Compute states on-the-fly using constant rates (original approach)
    // This avoids needing states_full_grid which is not computed when process noise is off
    for (i in 1:n_patients) {
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, i);
      int visit_size = get_pos_size(patient_visit_pos, i);

      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = get_pos_size(forecast_visits_pos, i);

      // Build forecast time WITH anchor (duplicate last observed week as element 1)
      array[forecast_size + 1] real forecast_time = linspaced_array(
        forecast_size + 1,
        patient_last_obs_visit[i],
        last_predict_visit);

      // Generate states using constant rates (computed on-the-fly, not from grid)
      matrix[forecast_size, 2] temp_forecast_patient_states;
      vector[visit_size] temp_rep_patient_log_sld;
      vector[visit_size] temp_rep_mean_patient_log_sld;
      vector[forecast_size] temp_forecast_patient_log_sld;
      vector[forecast_size] temp_forecast_mean_patient_log_sld;
      matrix[visit_size - 1, 2] temp_obs_process_noise;  // Unused but required by function

      (temp_forecast_patient_states, temp_rep_patient_log_sld, temp_rep_mean_patient_log_sld,
       temp_forecast_patient_log_sld, temp_forecast_mean_patient_log_sld, temp_obs_process_noise) =
        generate_patient_states_with_means_rng(
          states[visit_start:visit_end],  // Use states computed in transformed_parameters
          forecast_time,
          patient_log_decrease_rate[i, 1],  // Scalar rate (constant)
          patient_log_growth_rate[i, 1],    // Scalar rate (constant)
          sum_tumor_size[visit_start],
          negative_infinity(),  // growth lag (disabled)
          1.0,                  // growth transition
          rep_matrix(0.0, forecast_size, 2),  // No forecast process noise
          measure_sd_sld,
          measure_nu
        );

      // Store results
      forecast_patient_states[forecast_visit_start:forecast_visit_end] = temp_forecast_patient_states;
      rep_patient_log_sld[visit_start:visit_end] = temp_rep_patient_log_sld;
      rep_mean_patient_log_sld[visit_start:visit_end] = temp_rep_mean_patient_log_sld;
      forecast_patient_log_sld[forecast_visit_start:forecast_visit_end] = temp_forecast_patient_log_sld;
      forecast_mean_patient_log_sld[forecast_visit_start:forecast_visit_end] = temp_forecast_mean_patient_log_sld;
    }
  }
  
  // Weekly RECIST from noisy SLD (for visualization/diagnostics)
  (rep_recist, forecast_recist) = calculate_all_patients_recist(
    rep_patient_log_sld,
    forecast_patient_log_sld,
    patient_visit_pos,
    forecast_visits_pos,
    n_patient_screening_visits
  );

  // Assessment-visit noisy SLD for RECIST endpoint computation
  vector[n_total_forecast_obs_visits] forecast_obs_log_sld;
  array[n_total_forecast_obs_visits] int forecast_obs_recist;
  for (i in 1:n_patients) {
    int n_assessment = get_pos_size(forecast_obs_visits_pos, i);
    if (n_assessment > 0) {
      int assess_start, assess_end;
      (assess_start, assess_end) = get_pos(forecast_obs_visits_pos, i);

      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = forecast_visit_end - forecast_visit_start + 1;

      // Extract mean SLD at assessment weeks (subsample from weekly mean grid)
      vector[n_assessment] obs_visit_mean;
      for (a in 1:n_assessment) {
        int forecast_idx = min(a * forecast_observation_interval, forecast_size);
        obs_visit_mean[a] = forecast_mean_patient_log_sld[forecast_visit_start + forecast_idx - 1];
      }

      // Apply measurement noise at assessment visits only
      forecast_obs_log_sld[assess_start:assess_end] =
        to_vector(student_t_rng(measure_nu, obs_visit_mean, measure_sd_sld));
    }
  }

  // Assessment-visit RECIST from noisy SLD (for endpoint computation)
  {
    array[sum(n_patient_visits)] int unused_rep_recist;
    (unused_rep_recist, forecast_obs_recist) = calculate_all_patients_recist(
        rep_patient_log_sld,
        forecast_obs_log_sld,
        patient_visit_pos,
        forecast_obs_visits_pos,
        n_patient_screening_visits
    );
  }
  
  // Calculate patient-level endpoints (PFS, response) for all patients at once
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
   spop_os, spop_os_censored,
   spop_is_dropout, sample_is_dropout
  ) = calculate_all_patients_endpoints_rng(
    forecast_patient_idx,
    recist,
    rep_recist,
    forecast_recist,
    forecast_obs_recist,
    forecast_obs_visits_pos,
    forecast_observation_interval,
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
    pfs,
    right_censored,
    target_pfs,
    target_right_censored,
    ms_censored_01,
    ms_final_state,
    ms_time_02,
    ms_censored_02,
    ms_censored_12,
    ms_time_12,
    ms_os_event_12,
    ms_time_03,
    ms_time_32,
    ms_censored_32,
    patient_visit_pos,
    forecast_visits_pos,
    patient_last_obs_visit,
    last_predict_visit,
    t_patient_visits,
    max_all_t,
    n_patient_screening_visits,
    0,               // enable_ms_visit_gated_01 = 0 (tumor model does not use visit-gated 0->1)
    0.0,             // tv_coef_01_val
    zeros_vector(0), // forecast_obs_log_psa (unused)
    0.0,             // median_log_psa_obs (unused)
    0.0              // iqr_log_psa_obs (unused)
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
      forecast_trial_patient_pos,
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

// ── Competing Risks CIF (per-trial, empirical subdistribution) ───────────────
// Uses combined PFS (spop_pfs = min(target_recist, ms_hazard)) so that
// RECIST-detected progressions are included in the 0→1 cause.

array[n_trials] vector<lower=0, upper=1>[max_all_t + 1]
  spop_cif_01   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  spop_cif_02   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  spop_cif_03   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_01 = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_02 = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_03 = rep_array(zeros_vector(max_all_t + 1), n_trials);

for (s in 1:n_trials) {
  int n_tr = get_pos_size(trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(trial_patient_pos, s);

    (spop_cif_01[s], spop_cif_02[s], spop_cif_03[s]) = compute_trial_cif(
      spop_pfs[tr_start:tr_end], spop_right_censored[tr_start:tr_end],
      spop_is_dropout[tr_start:tr_end],
      spop_os[tr_start:tr_end],  spop_os_censored[tr_start:tr_end],
      max_all_t
    );

    (sample_cif_01[s], sample_cif_02[s], sample_cif_03[s]) = compute_trial_cif(
      sample_pfs[tr_start:tr_end], sample_right_censored[tr_start:tr_end],
      sample_is_dropout[tr_start:tr_end],
      sample_os[tr_start:tr_end],  sample_os_censored[tr_start:tr_end],
      max_all_t
    );
  }
}
