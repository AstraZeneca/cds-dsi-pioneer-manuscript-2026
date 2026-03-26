// ============================================================================
// Standalone Multistate Endpoint Simulation
//
// Simulates PFS and OS using only the illness-death multistate model,
// without tumor dynamics. PFS = 0→1 transition (competing with 0→2).
// OS = illness-death routing (0→2 direct death, or 0→1→2 via progression).
//
// "sample" = conditional on observed data, forecasting only censored patients
// "spop"   = unconditional posterior predictive from time 0
// ============================================================================

// ── Patient-level PFS and OS ──────────────────────────────────────────────
array[n_patients] int<lower=0> sample_ms_pfs, spop_ms_pfs;
array[n_patients] int<lower=0, upper=1> sample_ms_right_censored, spop_ms_right_censored;

array[n_patients] int<lower=0> sample_os, spop_os;
array[n_patients] int<lower=0, upper=1> sample_os_censored, spop_os_censored;

// Dropout flags — 1 if patient exited via cause 3 in this draw (for CIF computation)
array[n_patients] int<lower=0, upper=1> spop_is_dropout, sample_is_dropout;

// ── Composite PFS (progression OR direct death) ──────────────────────────
array[n_patients] int<lower=0> sample_pfs, spop_pfs;
array[n_patients] int<lower=0, upper=1> sample_pfs_right_censored, spop_pfs_right_censored;

// ── Trial-level PFS KM ───────────────────────────────────────────────────
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_ms_pfs_km_est, spop_ms_pfs_km_est;

array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);

array[n_trials, n_pfs_quantiles] int
  sample_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);

array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints]
  sample_ms_pfs_n, spop_ms_pfs_n;

// ── Trial-level OS KM ────────────────────────────────────────────────────
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_os_km_est, spop_os_km_est;

array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials);

array[n_trials, n_pfs_quantiles] int
  sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);

array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints]
  sample_os_n, spop_os_n;

// ── Composite PFS KM ────────────────────────────────────────────────────
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_pfs_km_est, spop_pfs_km_est;

array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials);

array[n_trials, n_pfs_quantiles] int
  sample_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);

array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints]
  sample_pfs_n, spop_pfs_n;

// ── Conditional group PFS KM ─────────────────────────────────────────────
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est;

array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_cond_group, n_pfs_quantiles] int
  cond_sample_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints]
  cond_sample_ms_pfs_n, cond_spop_ms_pfs_n;

// ── Conditional group OS KM ──────────────────────────────────────────────
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_os_km_est, cond_spop_os_km_est;

array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_cond_group, n_pfs_quantiles] int
  cond_sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints]
  cond_sample_os_n, cond_spop_os_n;

// ── Conditional group composite PFS KM ──────────────────────────────────
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_pfs_km_est, cond_spop_pfs_km_est;

array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_cond_group, n_pfs_quantiles] int
  cond_sample_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints]
  cond_sample_pfs_n, cond_spop_pfs_n;

// ══════════════════════════════════════════════════════════════════════════
// PATIENT-LEVEL SIMULATION
// ══════════════════════════════════════════════════════════════════════════

for (i in 1:n_patients) {

  // Guard: clock-forward matrix is [0,0] when ms_time_scale_12==1 (semi-Markov)
  row_vector[cols(log_cond_surv_12_t)] surv_12_t_i =
    cols(log_cond_surv_12_t) > 0 ? log_cond_surv_12_t[i] : rep_row_vector(0, 0);

  // ── Build visit schedule (shared by spop and sample 0→3) ────────────────
  int v_start = patient_visit_pos[i];
  int v_end = patient_visit_pos[i + 1] - 1;
  int n_obs_v = v_end - v_start + 1;
  int last_obs_wk = t_patient_visits[v_end];
  int n_fc_v = max(0, (max_all_t - last_obs_wk) %/% forecast_observation_interval);
  array[n_obs_v + n_fc_v] int patient_visits = build_spop_visit_schedule(
      t_patient_visits[v_start:v_end], forecast_observation_interval, max_all_t);

  // ── Unconditional (spop) pathway ─────────────────────────────────────────

  // Draw individual competing times
  int spop_t01; int spop_c01;
  if (enable_ms_visit_gated_01) {
    // Visit-gated: one Bernoulli per assessment visit (matches likelihood's
    // sum_at_visits_below). Events can only be detected AT visits.
    (spop_t01, spop_c01) = visit_only_survival_time_rng(
        log_cond_surv_01[i], patient_visits, max_all_t);
  } else {
    // Continuous: per-week Bernoulli (raw 0-indexed → +1 for 1-based week)
    int spop_t01_raw;
    (spop_t01_raw, spop_c01) = survival_time_rng(log_cond_surv_01[i]);
    spop_t01 = spop_t01_raw + 1;
  }

  int spop_t02 = max_all_t + 1; int spop_c02 = 1;
  if (enable_ms_02) {
    int t02_raw; int c02_raw;
    (t02_raw, c02_raw) = survival_time_rng(log_cond_surv_02[i]);
    spop_t02 = t02_raw + 1; spop_c02 = c02_raw;
  }

  int spop_t03 = max_all_t + 1; int spop_c03 = 1;
  if (enable_ms_03) {
    // Visit-only sampling: dropout can only occur AT visits, not between them.
    (spop_t03, spop_c03) = visit_only_survival_time_rng(
        log_cond_surv_03[i], patient_visits, max_all_t);
  }

  int spop_cause; int spop_exit;
  (spop_cause, spop_exit) = classify_spop_exit(
    spop_t01, spop_c01, spop_t02, spop_c02, spop_t03, spop_c03,
    enable_ms_02, enable_ms_03);
  spop_is_dropout[i] = (spop_cause == 3);

  // PFS: no RECIST in standalone — pass dummy target (max_all_t+1, censored=1)
  { int unused_pfs; int unused_cens;
    (unused_pfs, unused_cens, spop_ms_pfs[i], spop_ms_right_censored[i]) =
      derive_spop_pfs(spop_cause, spop_exit,
        max_all_t + 1, 1,    // no RECIST component
        spop_t01, spop_c01, spop_t03, enable_ms_02, enable_ms_03);
  }

  (spop_os[i], spop_os_censored[i]) = derive_spop_os_rng(
    spop_cause, spop_exit,
    enable_ms_12, enable_ms_32, ms_time_scale_12,
    log_cond_surv_12_s[i], surv_12_t_i, log_cond_surv_32[i],
    spop_t01, spop_t02, spop_t03);

  // ── Conditional (sample) pathway ─────────────────────────────────────────

  // 0→1: respect observed data; forecast only censored patients
  int sample_t01; int sample_c01;
  if (ms_censored_01[i]) {
    if (enable_ms_visit_gated_01) {
      // Visit-gated: forecast from observed censoring time, checking only future visits
      (sample_t01, sample_c01) = visit_only_survival_time_rng(
          log_cond_surv_01[i], patient_visits, max_all_t,
          ms_time_01[i], 1);
    } else {
      int t01_raw; int c01_raw;
      (t01_raw, c01_raw) = survival_time_rng(log_cond_surv_01[i], ms_time_01[i], 1, 0);
      sample_t01 = t01_raw + 1; sample_c01 = c01_raw;
    }
  } else {
    sample_t01 = ms_time_01[i]; sample_c01 = 0;
  }

  // 0→2: use observed event or forecast from observed censoring time
  int sample_t02 = max_all_t + 1; int sample_c02 = 1;
  if (enable_ms_02) {
    if (!ms_censored_02[i]) {
      sample_t02 = ms_time_02[i]; sample_c02 = 0;
    } else {
      int t02_raw; int c02_raw;
      (t02_raw, c02_raw) = survival_time_rng(log_cond_surv_02[i], ms_time_01[i], 1, 0);
      sample_t02 = t02_raw + 1; sample_c02 = c02_raw;
    }
  }

  // 0→3: use observed data or forecast for admin-censored patients
  int sample_t03 = max_all_t + 1; int sample_c03 = 1;
  if (enable_ms_03) {
    if (ms_final_state[i] == 3) {
      // Observed dropout — use ground truth
      sample_t03 = ms_time_03[i]; sample_c03 = 0;
    } else if (ms_final_state[i] == 0) {
      // Admin-censored — forecast 0→3 from last observation
      (sample_t03, sample_c03) = visit_only_survival_time_rng(
          log_cond_surv_03[i], patient_visits, max_all_t,
          ms_time_03[i], 1);
    }
    // States 1, 2: event observed → dropout didn't happen; leave as (max_all_t+1, 1)
  }

  int sample_cause; int sample_exit;
  (sample_cause, sample_exit) = classify_sample_exit(
    sample_t01, sample_c01, sample_t02, sample_c02,
    sample_t03, sample_c03,
    ms_censored_02[i],
    enable_ms_02, enable_ms_03);
  sample_is_dropout[i] = (sample_cause == 3);

  // pfs_km_est: composite PFS — cause 1 (progression) OR cause 2 (direct death) = event
  // cause 3 (dropout) or admin-censored = censored
  spop_pfs[i]                  = spop_exit;
  spop_pfs_right_censored[i]   = (spop_cause == 1 || spop_cause == 2) ? 0 : 1;
  sample_pfs[i]                = sample_exit;
  sample_pfs_right_censored[i] = (sample_cause == 1 || sample_cause == 2) ? 0 : 1;

  // PFS: no RECIST in standalone — pass dummy target (max_all_t+1, censored=1)
  { int unused_pfs; int unused_cens;
    (unused_pfs, unused_cens, sample_ms_pfs[i], sample_ms_right_censored[i]) =
      derive_sample_pfs(sample_cause, sample_exit,
        max_all_t + 1, 1,    // no RECIST component
        sample_t01, sample_c01, sample_t03, enable_ms_02, enable_ms_03);
  }

  (sample_os[i], sample_os_censored[i]) = derive_sample_os_rng(
    sample_cause, sample_exit,
    enable_ms_12, enable_ms_32, ms_time_scale_12,
    log_cond_surv_12_s[i], surv_12_t_i, log_cond_surv_32[i],
    sample_t01, sample_t02,
    ms_time_01[i],
    ms_censored_12[i], ms_time_12[i], ms_os_event_12[i],
    ms_censored_32[i], sample_t03, ms_time_32[i]);
}

// ══════════════════════════════════════════════════════════════════════════
// TRIAL-LEVEL AGGREGATION
// ══════════════════════════════════════════════════════════════════════════

for (s in 1:n_trials) {
  if (get_pos_size(trial_patient_pos, s) > 0) {
    // PFS KM curves
    sample_ms_pfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(sample_ms_pfs, trial_patient_pos, s),
      get_int_sub_array(sample_ms_right_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    spop_ms_pfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(spop_ms_pfs, trial_patient_pos, s),
      get_int_sub_array(spop_ms_right_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    // PFS quantiles
    (sample_ms_pfs_quant[s], sample_ms_pfs_quant_exceeds_max[s]) =
      km_quantiles(sample_ms_pfs_km_est[s], pfs_quantiles);
    (spop_ms_pfs_quant[s], spop_ms_pfs_quant_exceeds_max[s]) =
      km_quantiles(spop_ms_pfs_km_est[s], pfs_quantiles);

    // PFS-n at specified timepoints
    for (n in 1:n_pfs_timepoints) {
      sample_ms_pfs_n[s, n] = calc_km_pfs_n(sample_ms_pfs_km_est[s], months_to_weeks(pfs_timepoints[n]));
      spop_ms_pfs_n[s, n] = calc_km_pfs_n(spop_ms_pfs_km_est[s], months_to_weeks(pfs_timepoints[n]));
    }

    // OS KM curves
    sample_os_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(sample_os, trial_patient_pos, s),
      get_int_sub_array(sample_os_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    spop_os_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(spop_os, trial_patient_pos, s),
      get_int_sub_array(spop_os_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    // OS quantiles
    (sample_os_quant[s], sample_os_quant_exceeds_max[s]) =
      km_quantiles(sample_os_km_est[s], pfs_quantiles);
    (spop_os_quant[s], spop_os_quant_exceeds_max[s]) =
      km_quantiles(spop_os_km_est[s], pfs_quantiles);

    // OS-n at specified timepoints
    for (n in 1:n_pfs_timepoints) {
      sample_os_n[s, n] = calc_km_pfs_n(sample_os_km_est[s], months_to_weeks(pfs_timepoints[n]));
      spop_os_n[s, n] = calc_km_pfs_n(spop_os_km_est[s], months_to_weeks(pfs_timepoints[n]));
    }

    // Composite PFS KM
    sample_pfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(sample_pfs, trial_patient_pos, s),
      get_int_sub_array(sample_pfs_right_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    spop_pfs_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(spop_pfs, trial_patient_pos, s),
      get_int_sub_array(spop_pfs_right_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    (sample_pfs_quant[s], sample_pfs_quant_exceeds_max[s]) =
      km_quantiles(sample_pfs_km_est[s], pfs_quantiles);
    (spop_pfs_quant[s], spop_pfs_quant_exceeds_max[s]) =
      km_quantiles(spop_pfs_km_est[s], pfs_quantiles);

    for (n in 1:n_pfs_timepoints) {
      sample_pfs_n[s, n] = calc_km_pfs_n(sample_pfs_km_est[s], months_to_weeks(pfs_timepoints[n]));
      spop_pfs_n[s, n]   = calc_km_pfs_n(spop_pfs_km_est[s],   months_to_weeks(pfs_timepoints[n]));
    }
  } else {
    // Empty trial
    sample_ms_pfs_km_est[s] = zeros_vector(max_all_t + 1);
    spop_ms_pfs_km_est[s] = zeros_vector(max_all_t + 1);
    sample_ms_pfs_n[s] = zeros_vector(n_pfs_timepoints);
    spop_ms_pfs_n[s] = zeros_vector(n_pfs_timepoints);
    sample_os_km_est[s] = zeros_vector(max_all_t + 1);
    spop_os_km_est[s] = zeros_vector(max_all_t + 1);
    sample_os_n[s] = zeros_vector(n_pfs_timepoints);
    spop_os_n[s] = zeros_vector(n_pfs_timepoints);
    sample_pfs_km_est[s] = zeros_vector(max_all_t + 1);
    spop_pfs_km_est[s]   = zeros_vector(max_all_t + 1);
    sample_pfs_n[s]      = zeros_vector(n_pfs_timepoints);
    spop_pfs_n[s]        = zeros_vector(n_pfs_timepoints);
  }
}

// ══════════════════════════════════════════════════════════════════════════
// CONDITIONAL GROUP AGGREGATION
// ══════════════════════════════════════════════════════════════════════════

for (c in 1:n_cond_group) {
  int curr_group_size = get_pos_size(cond_group_pos, c);
  if (curr_group_size > 0) {
    array[curr_group_size] int curr_group_patients = get_int_sub_array(cond_group, cond_group_pos, c);

    // PFS KM
    cond_sample_ms_pfs_km_est[c] = estimate_kaplan_meier(
      sample_ms_pfs[curr_group_patients],
      sample_ms_right_censored[curr_group_patients],
      max_all_t, 0).1;

    cond_spop_ms_pfs_km_est[c] = estimate_kaplan_meier(
      spop_ms_pfs[curr_group_patients],
      spop_ms_right_censored[curr_group_patients],
      max_all_t, 0).1;

    // PFS quantiles
    (cond_sample_ms_pfs_quant[c], cond_sample_ms_pfs_quant_exceeds_max[c]) =
      km_quantiles(cond_sample_ms_pfs_km_est[c], pfs_quantiles);
    (cond_spop_ms_pfs_quant[c], cond_spop_ms_pfs_quant_exceeds_max[c]) =
      km_quantiles(cond_spop_ms_pfs_km_est[c], pfs_quantiles);

    // PFS-n
    for (n in 1:n_pfs_timepoints) {
      cond_sample_ms_pfs_n[c, n] = calc_km_pfs_n(cond_sample_ms_pfs_km_est[c], months_to_weeks(pfs_timepoints[n]));
      cond_spop_ms_pfs_n[c, n] = calc_km_pfs_n(cond_spop_ms_pfs_km_est[c], months_to_weeks(pfs_timepoints[n]));
    }

    // OS KM
    cond_sample_os_km_est[c] = estimate_kaplan_meier(
      sample_os[curr_group_patients],
      sample_os_censored[curr_group_patients],
      max_all_t, 0).1;

    cond_spop_os_km_est[c] = estimate_kaplan_meier(
      spop_os[curr_group_patients],
      spop_os_censored[curr_group_patients],
      max_all_t, 0).1;

    // OS quantiles
    (cond_sample_os_quant[c], cond_sample_os_quant_exceeds_max[c]) =
      km_quantiles(cond_sample_os_km_est[c], pfs_quantiles);
    (cond_spop_os_quant[c], cond_spop_os_quant_exceeds_max[c]) =
      km_quantiles(cond_spop_os_km_est[c], pfs_quantiles);

    // OS-n
    for (n in 1:n_pfs_timepoints) {
      cond_sample_os_n[c, n] = calc_km_pfs_n(cond_sample_os_km_est[c], months_to_weeks(pfs_timepoints[n]));
      cond_spop_os_n[c, n] = calc_km_pfs_n(cond_spop_os_km_est[c], months_to_weeks(pfs_timepoints[n]));
    }

    // Conditional group composite PFS KM
    cond_sample_pfs_km_est[c] = estimate_kaplan_meier(
      sample_pfs[curr_group_patients],
      sample_pfs_right_censored[curr_group_patients],
      max_all_t, 0).1;

    cond_spop_pfs_km_est[c] = estimate_kaplan_meier(
      spop_pfs[curr_group_patients],
      spop_pfs_right_censored[curr_group_patients],
      max_all_t, 0).1;

    (cond_sample_pfs_quant[c], cond_sample_pfs_quant_exceeds_max[c]) =
      km_quantiles(cond_sample_pfs_km_est[c], pfs_quantiles);
    (cond_spop_pfs_quant[c], cond_spop_pfs_quant_exceeds_max[c]) =
      km_quantiles(cond_spop_pfs_km_est[c], pfs_quantiles);

    for (n in 1:n_pfs_timepoints) {
      cond_sample_pfs_n[c, n] = calc_km_pfs_n(cond_sample_pfs_km_est[c], months_to_weeks(pfs_timepoints[n]));
      cond_spop_pfs_n[c, n]   = calc_km_pfs_n(cond_spop_pfs_km_est[c],   months_to_weeks(pfs_timepoints[n]));
    }
  } else {
    // Empty group
    cond_sample_ms_pfs_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_ms_pfs_km_est[c] = zeros_vector(max_all_t + 1);
    cond_sample_ms_pfs_n[c] = zeros_vector(n_pfs_timepoints);
    cond_spop_ms_pfs_n[c] = zeros_vector(n_pfs_timepoints);
    cond_sample_os_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_os_km_est[c] = zeros_vector(max_all_t + 1);
    cond_sample_os_n[c] = zeros_vector(n_pfs_timepoints);
    cond_spop_os_n[c] = zeros_vector(n_pfs_timepoints);
    cond_sample_pfs_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_pfs_km_est[c]   = zeros_vector(max_all_t + 1);
    cond_sample_pfs_n[c]      = zeros_vector(n_pfs_timepoints);
    cond_spop_pfs_n[c]        = zeros_vector(n_pfs_timepoints);
  }
}

// ── Competing Risks CIF (per-trial, empirical subdistribution) ───────────────
// Standalone model only has multistate PFS (no SLD/RECIST pathway).

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
      spop_ms_pfs[tr_start:tr_end], spop_ms_right_censored[tr_start:tr_end],
      spop_is_dropout[tr_start:tr_end],
      spop_os[tr_start:tr_end],  spop_os_censored[tr_start:tr_end],
      max_all_t
    );

    (sample_cif_01[s], sample_cif_02[s], sample_cif_03[s]) = compute_trial_cif(
      sample_ms_pfs[tr_start:tr_end], sample_ms_right_censored[tr_start:tr_end],
      sample_is_dropout[tr_start:tr_end],
      sample_os[tr_start:tr_end],  sample_os_censored[tr_start:tr_end],
      max_all_t
    );
  }
}
