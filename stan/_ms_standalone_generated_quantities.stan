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

// ── Trial-level PFS KM ───────────────────────────────────────────────────
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_ms_km_est, spop_ms_km_est;

array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials);

array[n_trials, n_pfs_quantiles] int
  sample_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);

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

// ── Conditional group PFS KM ─────────────────────────────────────────────
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_ms_km_est, cond_spop_ms_km_est;

array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);

array[n_cond_group, n_pfs_quantiles] int
  cond_sample_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

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

// ══════════════════════════════════════════════════════════════════════════
// PATIENT-LEVEL SIMULATION
// ══════════════════════════════════════════════════════════════════════════

for (i in 1:n_patients) {

  // ── 1. PFS via 0→1 transition ────────────────────────────────────────

  // Unconditional: full posterior predictive from time 0
  (spop_ms_pfs[i], spop_ms_right_censored[i]) = survival_time_rng(log_cond_surv_01[i]);
  spop_ms_pfs[i] += 1;  // Convert to 1-based detection week

  // Conditional: respect observed data, forecast only censored patients
  if (ms_censored_01[i]) {
    // Censored for 0→1: forecast from observed time
    (sample_ms_pfs[i], sample_ms_right_censored[i]) =
      survival_time_rng(log_cond_surv_01[i], ms_time_01[i], 1, 0);
    sample_ms_pfs[i] += 1;
  } else {
    // Event observed: use data (pfs + interval_censored + 1 = detection week)
    sample_ms_pfs[i] = ms_time_01[i] + interval_censored[i] + 1;
    sample_ms_right_censored[i] = 0;
  }

  // ── 2. OS via competing risks (0→1 vs 0→2 vs 0→3) + illness-death routing ──

  // Guard: clock-forward matrix is [0,0] when ms_time_scale_12==1 (semi-Markov)
  row_vector[cols(log_cond_surv_12_t)] surv_12_t_i =
    cols(log_cond_surv_12_t) > 0 ? log_cond_surv_12_t[i] : rep_row_vector(0, 0);

  // --- Unconditional (spop) OS ---
  {
    int spop_time_01 = spop_ms_pfs[i];
    int spop_cens_01 = spop_ms_right_censored[i];

    int spop_time_02 = max_all_t + 1;
    int spop_cens_02 = 1;
    if (enable_ms_02) {
      int t02_raw; int c02_raw;
      (t02_raw, c02_raw) = survival_time_rng(log_cond_surv_02[i]);
      spop_time_02 = t02_raw + 1;
      spop_cens_02 = c02_raw;
    }

    int spop_time_03 = max_all_t + 1;
    int spop_cens_03 = 1;
    if (enable_ms_03) {
      int t03_raw; int c03_raw;
      (t03_raw, c03_raw) = survival_time_rng(log_cond_surv_03[i]);
      spop_time_03 = t03_raw + 1;
      spop_cens_03 = c03_raw;
    }

    // Competing risks truth table
    // Priority on ties: progression (0→1) > direct death (0→2) > dropout (0→3)
    int spop_progressed_first = !spop_cens_01
      && (spop_cens_02 || spop_time_01 <= spop_time_02)
      && (spop_cens_03 || spop_time_01 <= spop_time_03);
    int spop_died_directly = enable_ms_02 && !spop_cens_02
      && (spop_cens_01 || spop_time_02 < spop_time_01)
      && (spop_cens_03 || spop_time_02 <= spop_time_03);
    int spop_dropped_out = enable_ms_03 && !spop_cens_03
      && (spop_cens_01 || spop_time_03 < spop_time_01)
      && (spop_cens_02 || spop_time_03 < spop_time_02);

    if (spop_died_directly) {
      spop_os[i] = spop_time_02;
      spop_os_censored[i] = 0;
      // Death without progression is a PFS event
      spop_ms_pfs[i] = spop_time_02;
      spop_ms_right_censored[i] = 0;
    } else if (spop_dropped_out && enable_ms_32) {
      // Dropped out: sample post-dropout death (3→2) unconditionally
      (spop_os[i], spop_os_censored[i]) = sample_dropout_death_rng(
        log_cond_surv_32[i], spop_time_03, 0);  // sojourn_obs=0: unconditional
      // Dropout censors PFS at dropout time
      spop_ms_pfs[i] = spop_time_03;
      spop_ms_right_censored[i] = 1;
    } else if (spop_dropped_out) {
      // Dropped out but 3→2 not modeled: censor OS at dropout
      spop_os[i] = spop_time_03;
      spop_os_censored[i] = 1;
      spop_ms_pfs[i] = spop_time_03;
      spop_ms_right_censored[i] = 1;
    } else if (spop_progressed_first && enable_ms_12) {
      (spop_os[i], spop_os_censored[i]) = sample_post_progression_death_rng(
        ms_time_scale_12, log_cond_surv_12_s[i], surv_12_t_i,
        spop_time_01, 0);  // sojourn_obs=0: unconditional
    } else {
      // All censored
      spop_os[i] = max({spop_time_01, spop_time_02, spop_time_03});
      spop_os_censored[i] = 1;
    }
  }

  // --- Conditional (sample) OS ---
  {
    int sample_time_01 = sample_ms_pfs[i];
    int sample_cens_01 = sample_ms_right_censored[i];

    int sample_time_02 = max_all_t + 1;
    int sample_cens_02 = 1;
    if (enable_ms_02) {
      if (!ms_censored_02[i]) {
        // Death without progression was observed
        sample_time_02 = ms_time_02[i];
        sample_cens_02 = 0;
      } else {
        // Censored for 0→2: forecast from observed time
        int t02_raw; int c02_raw;
        (t02_raw, c02_raw) = survival_time_rng(log_cond_surv_02[i], ms_time_01[i], 1, 0);
        sample_time_02 = t02_raw + 1;
        sample_cens_02 = c02_raw;
      }
    }

    // State 3: observed dropout is ground truth — takes priority over forecasts
    int sample_dropped_out = enable_ms_03 && ms_final_state[i] == 3;

    // Routing truth table
    int sample_died_directly = !sample_dropped_out && enable_ms_02 && !sample_cens_02
      && (!ms_censored_02[i] || sample_cens_01 || sample_time_02 < sample_time_01);
    int sample_progressed_first = !sample_dropped_out && !sample_died_directly && !sample_cens_01
      && (sample_cens_02 || sample_time_01 <= sample_time_02);

    if (sample_dropped_out) {
      if (!ms_censored_32[i]) {
        // Off-trial death was observed: use exact calendar time
        sample_os[i] = ms_time_03[i] + ms_time_32[i];
        sample_os_censored[i] = 0;
      } else if (enable_ms_32) {
        // Forecast post-dropout death conditioning on observed sojourn survival
        (sample_os[i], sample_os_censored[i]) = sample_dropout_death_rng(
          log_cond_surv_32[i], ms_time_03[i], ms_time_32[i]);
      } else {
        // 3→2 not modeled: censor at last off-trial observation
        sample_os[i] = ms_time_03[i] + ms_time_32[i];
        sample_os_censored[i] = 1;
      }
      // Dropout censors PFS at dropout time
      sample_ms_pfs[i] = ms_time_03[i];
      sample_ms_right_censored[i] = 1;
    } else if (sample_died_directly) {
      sample_os[i] = sample_time_02;
      sample_os_censored[i] = 0;
      // Death without progression is a PFS event
      sample_ms_pfs[i] = sample_time_02;
      sample_ms_right_censored[i] = 0;
    } else if (sample_progressed_first && enable_ms_12) {
      if (!ms_censored_12[i]) {
        // Observed post-progression death: use exact time
        sample_os[i] = ms_time_01[i] + ms_time_12[i];
        sample_os_censored[i] = 0;
      } else {
        // Forecast post-progression death, conditioning on observed sojourn survival
        (sample_os[i], sample_os_censored[i]) = sample_post_progression_death_rng(
          ms_time_scale_12, log_cond_surv_12_s[i], surv_12_t_i,
          sample_time_01, ms_time_12[i]);
      }
    } else {
      // All censored
      sample_os[i] = max(sample_time_01, sample_time_02);
      sample_os_censored[i] = 1;
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════
// TRIAL-LEVEL AGGREGATION
// ══════════════════════════════════════════════════════════════════════════

for (s in 1:n_trials) {
  if (get_pos_size(trial_patient_pos, s) > 0) {
    // PFS KM curves
    sample_ms_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(sample_ms_pfs, trial_patient_pos, s),
      get_int_sub_array(sample_ms_right_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    spop_ms_km_est[s] = estimate_kaplan_meier(
      get_int_sub_array(spop_ms_pfs, trial_patient_pos, s),
      get_int_sub_array(spop_ms_right_censored, trial_patient_pos, s),
      max_all_t, 0).1;

    // PFS quantiles
    (sample_ms_quant_pfs[s], sample_ms_quant_pfs_exceeds_max[s]) =
      km_quantiles(sample_ms_km_est[s], pfs_quantiles);
    (spop_ms_quant_pfs[s], spop_ms_quant_pfs_exceeds_max[s]) =
      km_quantiles(spop_ms_km_est[s], pfs_quantiles);

    // PFS-n at specified timepoints
    for (n in 1:n_pfs_timepoints) {
      sample_ms_pfs_n[s, n] = calc_km_pfs_n(sample_ms_km_est[s], months_to_weeks(pfs_timepoints[n]));
      spop_ms_pfs_n[s, n] = calc_km_pfs_n(spop_ms_km_est[s], months_to_weeks(pfs_timepoints[n]));
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
  } else {
    // Empty trial
    sample_ms_km_est[s] = zeros_vector(max_all_t + 1);
    spop_ms_km_est[s] = zeros_vector(max_all_t + 1);
    sample_ms_pfs_n[s] = zeros_vector(n_pfs_timepoints);
    spop_ms_pfs_n[s] = zeros_vector(n_pfs_timepoints);
    sample_os_km_est[s] = zeros_vector(max_all_t + 1);
    spop_os_km_est[s] = zeros_vector(max_all_t + 1);
    sample_os_n[s] = zeros_vector(n_pfs_timepoints);
    spop_os_n[s] = zeros_vector(n_pfs_timepoints);
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
    cond_sample_ms_km_est[c] = estimate_kaplan_meier(
      sample_ms_pfs[curr_group_patients],
      sample_ms_right_censored[curr_group_patients],
      max_all_t, 0).1;

    cond_spop_ms_km_est[c] = estimate_kaplan_meier(
      spop_ms_pfs[curr_group_patients],
      spop_ms_right_censored[curr_group_patients],
      max_all_t, 0).1;

    // PFS quantiles
    (cond_sample_ms_quant_pfs[c], cond_sample_ms_quant_pfs_exceeds_max[c]) =
      km_quantiles(cond_sample_ms_km_est[c], pfs_quantiles);
    (cond_spop_ms_quant_pfs[c], cond_spop_ms_quant_pfs_exceeds_max[c]) =
      km_quantiles(cond_spop_ms_km_est[c], pfs_quantiles);

    // PFS-n
    for (n in 1:n_pfs_timepoints) {
      cond_sample_ms_pfs_n[c, n] = calc_km_pfs_n(cond_sample_ms_km_est[c], months_to_weeks(pfs_timepoints[n]));
      cond_spop_ms_pfs_n[c, n] = calc_km_pfs_n(cond_spop_ms_km_est[c], months_to_weeks(pfs_timepoints[n]));
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
  } else {
    // Empty group
    cond_sample_ms_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_ms_km_est[c] = zeros_vector(max_all_t + 1);
    cond_sample_ms_pfs_n[c] = zeros_vector(n_pfs_timepoints);
    cond_spop_ms_pfs_n[c] = zeros_vector(n_pfs_timepoints);
    cond_sample_os_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_os_km_est[c] = zeros_vector(max_all_t + 1);
    cond_sample_os_n[c] = zeros_vector(n_pfs_timepoints);
    cond_spop_os_n[c] = zeros_vector(n_pfs_timepoints);
  }
}
