// ============================================================================
// Radiographic PFS (rPFS) Generated Quantities
//
// Prerequisites (in scope from pioneer.stan GQ block):
//   sample_ms_pfs, sample_ms_right_censored   — computed by burden_endpoints.stan
//   spop_ms_pfs,   spop_ms_right_censored     — computed by burden_endpoints.stan
//   sample_os,     sample_os_censored         — computed by burden_endpoints.stan
//   spop_os,       spop_os_censored           — computed by burden_endpoints.stan
//   forecast_trial_patient_pos, cond_group_pos, cond_group — from transformed data
//   max_all_t, n_pfs_timepoints, n_pfs_quantiles           — from data/flags
//   pfs_quantiles, pfs_timepoints                          — from data
// ============================================================================

// Patient-level rPFS arrays
array[n_forecast_patients] int<lower=0> sample_rpfs, spop_rpfs;
array[n_forecast_patients] int<lower=0, upper=1> sample_rpfs_censored, spop_rpfs_censored;

// Trial-level rPFS KM curves
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_rpfs_km_est, spop_rpfs_km_est;
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_rpfs_km_est, cond_spop_rpfs_km_est;

// rPFS quantiles — trial-level
array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_rpfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_rpfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int
  sample_rpfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_rpfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);

// rPFS quantiles — conditional group-level
array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_rpfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_rpfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int
  cond_sample_rpfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_rpfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// rPFS at fixed timepoints
array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints] sample_rpfs_n, spop_rpfs_n;
array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints] cond_sample_rpfs_n, cond_spop_rpfs_n;

// Patient-level: rPFS = min(radiographic PFS, OS), censored otherwise
for (i in 1:n_forecast_patients) {
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

// Trial-level aggregation
for (s in 1:n_trials) {
  if (get_pos_size(forecast_trial_patient_pos, s) > 0) {
    sample_rpfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(sample_rpfs, forecast_trial_patient_pos, s),
      get_int_sub_array(sample_rpfs_censored, forecast_trial_patient_pos, s),
      max_all_t, 0).1;
    spop_rpfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(spop_rpfs, forecast_trial_patient_pos, s),
      get_int_sub_array(spop_rpfs_censored, forecast_trial_patient_pos, s),
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

// Conditional group aggregation
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
