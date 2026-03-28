// ============================================================================
// STANDALONE MULTISTATE ENDPOINT SIMULATION
// ============================================================================
// Uses burden_endpoints.stan (shared with pioneer and tumor models).
// No biomarker dynamics: dummy all-zero categories ensure target_* outputs
// are all-censored. The meaningful outputs are ms_*, pfs_*, and os_* families.
// ============================================================================

// ── Patient-level endpoint arrays ───────────────────────────────────────────
array[n_hmc_patients] int<lower=0> sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs,
                                   spop_target_obs_cens_pfs, sample_pfs, spop_pfs;
array[n_hmc_patients] int<lower=0, upper=1>
  sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
  sample_ms_right_censored, spop_ms_right_censored,
  sample_right_censored, spop_right_censored;
array[n_hmc_patients] int<lower=0> sample_os, spop_os;
array[n_hmc_patients] int<lower=0, upper=1> sample_os_censored, spop_os_censored;
array[n_hmc_patients] int<lower=0, upper=1> spop_is_dropout, sample_is_dropout;
// In standalone all patients are target-censored (no biomarker-based target PFS).
array[n_hmc_patients] int<lower=0> forecast_target_pfs;
array[n_hmc_patients] int<lower=0, upper=1> forecast_target_right_censored;
array[n_hmc_patients] int<lower=0, upper=1> sample_target_confirmed_response, spop_target_confirmed_response;
array[n_hmc_patients] int<lower=0, upper=1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;

// ── Trial-level KM curves ───────────────────────────────────────────────────
vector<lower=0, upper=1>[n_trials] sample_target_orr, spop_target_orr;
vector<lower=0, upper=1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group),
                                       cond_spop_target_orr   = zeros_vector(n_cond_group);
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_target_km_est, spop_target_km_est,
                                                         spop_target_obs_cens_km_est,
                                                         sample_ms_km_est, spop_ms_km_est,
                                                         sample_km_est, spop_km_est;
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1] cond_sample_target_km_est, cond_spop_target_km_est,
                                                             cond_spop_target_obs_cens_km_est,
                                                             cond_sample_ms_km_est, cond_spop_ms_km_est,
                                                             cond_sample_km_est, cond_spop_km_est;
array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n,
                                                            sample_ms_pfs_n, spop_ms_pfs_n,
                                                            sample_pfs_n, spop_pfs_n;
array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n,
                                                                cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
                                                                cond_sample_pfs_n, cond_spop_pfs_n;
array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_target_quant_pfs   = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  sample_ms_quant_pfs     = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_ms_quant_pfs       = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  sample_quant_pfs        = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_quant_pfs          = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_target_quant_pfs = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_target_quant_pfs   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_sample_ms_quant_pfs     = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_quant_pfs       = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_sample_quant_pfs        = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_quant_pfs          = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_trials, n_pfs_quantiles] int
  sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_target_quant_pfs_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  sample_ms_quant_pfs_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_ms_quant_pfs_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  sample_quant_pfs_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_quant_pfs_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group, n_pfs_quantiles] int
  cond_sample_target_quant_pfs_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_target_quant_pfs_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_sample_ms_quant_pfs_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_ms_quant_pfs_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_sample_quant_pfs_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_quant_pfs_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// ── OS KM ───────────────────────────────────────────────────────────────────
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_os_km_est, spop_os_km_est;
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1] cond_sample_os_km_est, cond_spop_os_km_est;
array[n_trials] vector<lower=0>[n_pfs_quantiles]
  sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
  spop_os_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
array[n_trials, n_pfs_quantiles] int
  sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
  spop_os_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_os_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int
  cond_sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_os_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
array[n_trials] vector<lower=0, upper=1>[n_pfs_timepoints] sample_os_n, spop_os_n;
array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints] cond_sample_os_n, cond_spop_os_n;

// ── Sojourn KMs (1→2 post-progression, 3→2 off-trial) ───────────────────────
// sample: conditional on observed progression/dropout (patients who progressed/
//         dropped out in the observed data contribute observed sojourn times;
//         those still in state 0/1 at follow-up contribute censored sojourn times
//         drawn from the posterior).
// spop:   unconditional — every patient contributes a simulated sojourn time
//         drawn from the posterior (whether or not they actually progressed/dropped
//         out in the data).
array[n_trials] vector<lower=0, upper=1>[ms_max_sojourn_t + 1]
  sample_km_12 = rep_array(ones_vector(ms_max_sojourn_t + 1), n_trials),
  spop_km_12   = rep_array(ones_vector(ms_max_sojourn_t + 1), n_trials);
array[n_trials] vector<lower=0, upper=1>[ms_max_sojourn_t_32 + 1]
  sample_km_32 = rep_array(ones_vector(ms_max_sojourn_t_32 + 1), n_trials),
  spop_km_32   = rep_array(ones_vector(ms_max_sojourn_t_32 + 1), n_trials);

// ── Competing Risks CIF ──────────────────────────────────────────────────────
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1]
  spop_cif_01   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  spop_cif_02   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  spop_cif_03   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_01 = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_02 = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_03 = rep_array(zeros_vector(max_all_t + 1), n_trials);

// ── Shared endpoint computation ──────────────────────────────────────────────
// Contract vars are local (not saved to CSV). No biomarker: dummy categories
// never trigger target PD or response detection. burden_pfs combines
// progression (0→1) and death (0→2) as PFS events; all patients are
// target-censored so km_est == ms_km_est in the standalone.
{
  array[n_total_visits] int obs_biomarker_cat          = rep_array(0, n_total_visits);
  array[n_total_visits] int rep_biomarker_cat          = rep_array(0, n_total_visits);
  array[n_total_forecast_obs_visits] int forecast_obs_biomarker_cat = rep_array(0, n_total_forecast_obs_visits);
  // PFS = progression (0→1) OR death (0→2); dropout (0→3) is censored
  array[n_patients] int burden_pfs;
  array[n_patients] int burden_right_censored;
  for (i in 1:n_patients) {
    if (ms_censored_01[i] == 0) {
      burden_pfs[i] = ms_time_01[i];
      burden_right_censored[i] = 0;
    } else if (ms_censored_02[i] == 0) {
      burden_pfs[i] = ms_time_02[i];
      burden_right_censored[i] = 0;
    } else {
      burden_pfs[i] = min(ms_time_01[i], ms_time_02[i]);
      burden_right_censored[i] = 1;
    }
  }
  array[n_patients] int burden_target_pfs              = rep_array(max_all_t + 1, n_patients);
  array[n_patients] int burden_target_right_censored   = rep_array(1, n_patients);
  int  burden_enable_ms_visit_gated_01 = 0;
  real burden_tv_coef_01_val           = 0.0;
  real burden_median_log_psa_obs       = 0.0;
  real burden_iqr_log_psa_obs          = 0.0;
  vector[0] burden_forecast_obs_log_psa;

  #include "modules/state_space/burden_endpoints.stan"
}

// ── Sojourn KM computation ───────────────────────────────────────────────────
// Runs after the local block so GQ-scope endpoint arrays are populated.
for (s in 1:n_trials) {
  int n_tr = get_pos_size(hmc_trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(hmc_trial_patient_pos, s);

    // ── 1→2: post-progression sojourn KM ─────────────────────────────────────
    // sample: patients who progressed (sample_ms_right_censored == 0)
    {
      int n_prog = n_tr - sum(sample_ms_right_censored[tr_start:tr_end]);
      if (n_prog > 0) {
        array[n_prog] int soj; array[n_prog] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (!sample_ms_right_censored[i]) {
            soj[idx]  = max(1, sample_os[i] - sample_ms_pfs[i]);
            cens[idx] = sample_os_censored[i];
            idx += 1;
          }
        }
        sample_km_12[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
      }
    }

    // spop: all patients get a simulated post-progression sojourn
    {
      int n_prog = n_tr - sum(spop_ms_right_censored[tr_start:tr_end]);
      if (n_prog > 0) {
        array[n_prog] int soj; array[n_prog] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (!spop_ms_right_censored[i]) {
            soj[idx]  = max(1, spop_os[i] - spop_ms_pfs[i]);
            cens[idx] = spop_os_censored[i];
            idx += 1;
          }
        }
        spop_km_12[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
      }
    }

    // ── 3→2: off-trial sojourn KM ─────────────────────────────────────────────
    // sample: patients who dropped out (sample_is_dropout == 1)
    {
      int n_drop = sum(to_array_1d(sample_is_dropout[tr_start:tr_end]));
      if (n_drop > 0) {
        array[n_drop] int soj; array[n_drop] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (sample_is_dropout[i]) {
            soj[idx]  = max(1, sample_os[i] - sample_pfs[i]);
            cens[idx] = sample_os_censored[i];
            idx += 1;
          }
        }
        sample_km_32[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
      }
    }

    // spop: all patients get a simulated off-trial sojourn
    {
      int n_drop = sum(spop_is_dropout[tr_start:tr_end]);
      if (n_drop > 0) {
        array[n_drop] int soj; array[n_drop] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (spop_is_dropout[i]) {
            soj[idx]  = max(1, spop_os[i] - spop_pfs[i]);
            cens[idx] = spop_os_censored[i];
            idx += 1;
          }
        }
        spop_km_32[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
      }
    }
  }
}
