// --- Laplace surrogate: constants + validity guard --------------------------
// Constant inverse Vandermonde for the 3 fixed anchors (first anchor must be 0).
matrix[3, 3] surrogate_V;
for (k in 1:3) {
  surrogate_V[k, 1] = 1.0;
  surrogate_V[k, 2] = surrogate_anchor_times[k];
  surrogate_V[k, 3] = surrogate_anchor_times[k] * surrogate_anchor_times[k];
}
matrix[3, 3] surrogate_Vinv = inverse(surrogate_V);

// 3-point Gauss-Hermite nodes/weights for the standard normal (weights sum 1).
vector[3] surrogate_gh_x = [-sqrt(3.0), 0.0, sqrt(3.0)]';
vector[3] surrogate_gh_w = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]';

// Laplace inner-solver knobs (linear-Gaussian surrogate => solver 1; may fall
// back to solver 2 on the collinear (t,t^2) basis — harmless, validated).
real surrogate_tolerance = surrogate_tolerance_in > 0 ? surrogate_tolerance_in : 1e-8;
int surrogate_max_num_steps = surrogate_max_num_steps_in > 0 ? surrogate_max_num_steps_in : 100;
// d-parametric latent dim: 2 burden REs + one frailty latent per frailty slot.
// Detected from the MS intercept-mode config (patient level = last hierarchy col).
int surrogate_patient_lv = n_levels;
// Coupled slots = transitions whose hazard reads the burden coupling. Publication
// = {01, 03}; both read ms_time_varying_covar_01 (transformed_parameters.stan:176,780).
// We detect frailty slots as coupled slots with a patient-level RE-NCP intercept.
int surrogate_n_frailty_slots = 0;
array[2] int surrogate_frailty_slot = {0, 0};   // ascending slot order; 0 = unused
if (enable_background_surrogate == 1) {
  // MS_SLOT_01 / MS_SLOT_03 are in scope (modules/multistate/transformed_data.stan).
  if (enable_ms_01 && ms_level_intercept_mode[MS_SLOT_01, surrogate_patient_lv] == 2) {
    surrogate_n_frailty_slots += 1;
    surrogate_frailty_slot[surrogate_n_frailty_slots] = MS_SLOT_01;
  }
  if (enable_ms_03 && ms_level_intercept_mode[MS_SLOT_03, surrogate_patient_lv] == 2) {
    surrogate_n_frailty_slots += 1;
    surrogate_frailty_slot[surrogate_n_frailty_slots] = MS_SLOT_03;
  }
}
int surrogate_d = 2 + surrogate_n_frailty_slots;
int surrogate_hessian_block_size = surrogate_d;
int surrogate_solver = surrogate_solver_in > 0 ? surrogate_solver_in : 1;
int surrogate_max_steps_line_search = 0;
int surrogate_allow_fallback = 1;
real surrogate_jitter = 1e-10;

// GUARD: the bridge SD assumes the simple summed-variance patient SD, which
// holds for frac/init always and for tr UNLESS the tr SD sub-hierarchy is
// active (then the per-patient SD varies and the single bridge is wrong).
if (enable_background_surrogate == 1 && n_subhier_active_tr_intercept > 0)
  fatal_error("enable_background_surrogate=1 is incompatible with an active tr ",
              "SD sub-hierarchy (n_subhier_active_tr_intercept=",
              n_subhier_active_tr_intercept, "): the surrogate bridge assumes a ",
              "single patient-level tr SD. Disable the surrogate or the tr SD ",
              "sub-hierarchy.");

// First anchor must be 0 (the pinned intercept depends on g(0)=0).
if (enable_background_surrogate == 1 && surrogate_anchor_times[1] != 0.0)
  fatal_error("surrogate_anchor_times[1] must be 0 (baseline); got ",
              surrogate_anchor_times[1]);

// --- Joint-surrogate guards (fire only when the surrogate is on) -------------
int surrogate_frailty_block = 0;   // patient-level MS correlation block index (set by guard 4)
if (enable_background_surrogate == 1) {
  // Guard 6: split must be at a trial (non-patient) level — the baseline gather
  // depends on it (a patient-level split falls into the ->1 fallback).
  if (!(forecast_split_level > 0 && forecast_split_level < n_levels))
    fatal_error("enable_background_surrogate=1 requires 0 < forecast_split_level < ",
                "n_levels (got ", forecast_split_level, "); the surrogate baseline ",
                "gather is only valid for a trial-level split.");

  // Guards 1 & 2: every coupled/frailty slot's patient-level baseline must be a
  // marginalizable Gaussian intercept (mode 0 or 2) with NO patient-level GP.
  array[2] int surrogate_coupled_slot = {MS_SLOT_01, MS_SLOT_03};
  for (s_idx in 1:2) {
    int slot = surrogate_coupled_slot[s_idx];
    int slot_on = (slot == MS_SLOT_01) ? enable_ms_01 : enable_ms_03;
    if (slot_on) {
      if (enable_ms_level_gp[slot, surrogate_patient_lv] != 0)
        fatal_error("enable_background_surrogate=1 forbids a patient-level GP on ",
                    "coupled MS slot ", slot, " (not a marginalizable scalar).");
      int mode = ms_level_intercept_mode[slot, surrogate_patient_lv];
      if (!(mode == 0 || mode == 2))
        fatal_error("enable_background_surrogate=1: coupled MS slot ", slot,
                    " patient-level intercept mode must be 0 (off) or 2 (RE-NCP); ",
                    "got ", mode, ". Patient-level baseline beyond marginalized ",
                    "Gaussian frailty is unsupported by the surrogate.");
    }
  }

  // Guard 7: the joint functor accesses tv_coef_01[2] / tv_coef_03[2] (level +
  // velocity). time_varying_coef_01 is sized 1 when
  // (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01), else
  // n_time_varying_covar (parameters.stan:44). The surrogate requires the 2-feature
  // (level, velocity) basis => assert n_time_varying_covar == 2 AND not the size-1
  // observed-gating case. Publication: enable_ms_visit_gated_latent_01 = TRUE => size 2 (OK).
  if (n_time_varying_covar != 2)
    fatal_error("enable_background_surrogate=1 requires the 2-feature (level, ",
                "velocity) coupling basis (n_time_varying_covar=2); got ",
                n_time_varying_covar, ".");
  if (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01)
    fatal_error("enable_background_surrogate=1 is incompatible with OBSERVED ",
                "visit-gated 0->1 (enable_ms_visit_gated_latent_01=0), which sizes ",
                "time_varying_coef_01 to 1; the surrogate needs the latent 2-feature basis.");
  // Guard 7b: the functor unconditionally reads tv_coef_01[1..2] / tv_coef_03[1..2]
  // whenever a slot is a frailty slot. time_varying_coef_01 is size 0 unless
  // enable_ms_pop_time_varying_cov (parameters.stan:43); time_varying_coef_03 size 0
  // unless enable_ms_pop_time_varying_cov && enable_ms_03_time_varying_cov
  // (parameters.stan:149). Require coupling ON for every detected frailty slot
  // (frailty_slots subset of coupled-TV slots, per spec §2).
  if (surrogate_n_frailty_slots > 0 && !enable_ms_pop_time_varying_cov)
    fatal_error("enable_background_surrogate=1 with frailty slots requires ",
                "enable_ms_pop_time_varying_cov=1 (else time_varying_coef_01 is size 0 ",
                "and the functor reads tv_coef_01[2] out of bounds).");
  for (m in 1:surrogate_n_frailty_slots)
    if (surrogate_frailty_slot[m] == MS_SLOT_03 && !enable_ms_03_time_varying_cov)
      fatal_error("enable_background_surrogate=1 with 0->3 frailty requires ",
                  "enable_ms_03_time_varying_cov=1 (else time_varying_coef_03 is size 0).");

  // Guard 4: the patient-level correlation block's member slots must equal the
  // detected frailty_slots in the same (ascending) order, so Sigma_u rows align
  // with the theta frailty entries. Also record the block index for likelihood.stan.
  if (surrogate_n_frailty_slots >= 2) {
    for (b in 1:n_ms_corr_blocks) {
      if (ms_corr_block_level[b] == surrogate_patient_lv) {
        surrogate_frailty_block = b;
        for (m in 1:surrogate_n_frailty_slots)
          if (ms_corr_block_member_slot[b, m] != surrogate_frailty_slot[m])
            fatal_error("enable_background_surrogate=1: patient-level correlation ",
                        "block member slot order does not match detected frailty ",
                        "slots; Sigma_u would misalign with the latent vector.");
      }
    }
    if (surrogate_frailty_block == 0)
      fatal_error("enable_background_surrogate=1: expected a patient-level MS ",
                  "correlation block for ", surrogate_n_frailty_slots,
                  " frailty slots, found none.");
  }

  // Guard 8: the joint surrogate is currently implemented and validated ONLY for
  // the d=4 case frailty_slots == {01, 03} (the publication config). The d=3
  // single-frailty path would index L_ms_intercept_corr[surrogate_frailty_block]
  // with surrogate_frailty_block==0 (no corr block forms for a single member,
  // r/multistate.R:556), and an 01-only config would leave 0->3 functor inputs
  // (log_pop_lambda_03, tv_coef_03) degenerate while the 0->3 loop still runs.
  // Reject anything but {01,03} rather than silently mis-index/leak; revisit
  // deliberately when a new config needs it.
  if (surrogate_n_frailty_slots > 0
      && !(surrogate_n_frailty_slots == 2
           && surrogate_frailty_slot[1] == MS_SLOT_01
           && surrogate_frailty_slot[2] == MS_SLOT_03))
    fatal_error("enable_background_surrogate=1: the joint surrogate currently ",
                "supports only frailty_slots == {0->1, 0->3} (d=4). Detected ",
                surrogate_n_frailty_slots, " frailty slot(s) [",
                surrogate_frailty_slot[1], ", ", surrogate_frailty_slot[2],
                "]. Extend + re-validate the functor before enabling other sets.");

  // Guard 9 (review wf_e88e7d91 uncertain item -> hard guard): the surrogate
  // hazard anchors burden at the BASELINE (last-screening) week, while the forecast
  // covariate frame anchors at the FIRST visit (states_start_col = 2 - first_visit,
  // _ms_burden_tv_covar.stan:49). These coincide only when each background patient
  // has exactly ONE screening visit. The general visit guard only enforces >= 1
  // (_visit_transformed_data.stan:23), so assert == 1 for background patients to
  // keep the shared tv_coef applied at matched trajectory points.
  if (n_background_patients > 0) {
    for (j in 1:n_background_patients) {
      int p = background_patient_idx[j];
      if (n_patient_screening_visits[p] != 1)
        fatal_error("enable_background_surrogate=1: background patient ", p,
                    " has ", n_patient_screening_visits[p], " screening visits; ",
                    "the surrogate requires exactly 1 so its baseline-anchored ",
                    "burden frame matches the forecast first-visit frame.");
    }
  }
}

vector[enable_background_surrogate == 1 ? n_background_patients * surrogate_d : 0]
  surrogate_theta_0 = rep_vector(0.0,
    enable_background_surrogate == 1 ? n_background_patients * surrogate_d : 0);

// Background-only compact views, assembled from transformed data only (so they
// satisfy the data-only argument qualifiers of surrogate_ll). Per-patient LOD
// offset (burden normalized to each patient's baseline => LOD shifts by
// log_baseline_sld[p]).
array[enable_background_surrogate == 1 ? n_background_patients + 1 : 0] int surrogate_bg_pos;
int surrogate_n_bg_visits = 0;
if (enable_background_surrogate == 1 && n_background_patients > 0) {
  surrogate_bg_pos[1] = 1;
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    surrogate_bg_pos[j + 1] = surrogate_bg_pos[j] + (ve - vs + 1);
  }
  surrogate_n_bg_visits = surrogate_bg_pos[n_background_patients + 1] - 1;
}
vector[surrogate_n_bg_visits] surrogate_bg_obs;
array[surrogate_n_bg_visits] int surrogate_bg_time;
vector[surrogate_n_bg_visits] surrogate_bg_log_lod;
if (enable_background_surrogate == 1 && n_background_patients > 0) {
  int w = 1;
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    for (v in vs:ve) {
      surrogate_bg_obs[w]     = normalized_sld[v];
      surrogate_bg_time[w]    = t_patient_visit_idx[v];
      surrogate_bg_log_lod[w] = log_lod - log_baseline_sld[p];
      w += 1;
    }
  }
}

// --- bg-compact MS arrays for the joint surrogate hazard terms ---------------
// Built only when frailty/coupling is active; otherwise zero-sized (d=2 path).
int surrogate_ms_active = (enable_background_surrogate == 1 && surrogate_n_frailty_slots > 0) ? 1 : 0;
int surrogate_n_wk = surrogate_ms_active ? max_all_t : 0;

// Per-bg event week + censor flag for 0->1 and 0->3, and the per-bg baseline week
// (calendar week of the patient's baseline visit) used to anchor g(tau).
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_event_wk_01;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_censored_01;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_event_wk_03;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_censored_03;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_baseline_week;
// 0->1 visit-gating mask: [bg, week] == 1 at observed visit weeks.
array[surrogate_ms_active ? n_background_patients : 0,
      surrogate_ms_active ? max_all_t : 0] int surrogate_bg_visit_wk_01
  = rep_array(0, surrogate_ms_active ? n_background_patients : 0,
                 surrogate_ms_active ? max_all_t : 0);

if (surrogate_ms_active) {
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    // 0->1 event/censor: ms_censored_01 is a raw data field.
    surrogate_bg_censored_01[j] = ms_censored_01[p];
    surrogate_bg_event_wk_01[j] = ms_censored_01[p] ? max_all_t : ms_time_01[p];
    // 0->3 event/censor: derived. ms_final_state == 3 => dropout event at ms_time_03;
    // otherwise censored at ms_time_03 (= patient_max_t).
    surrogate_bg_censored_03[j] = (ms_final_state[p] == 3) ? 0 : 1;
    surrogate_bg_event_wk_03[j] = ms_time_03[p];
    // Baseline week = calendar week of the patient's BASELINE visit = the LAST
    // screening visit, NOT the first visit. Matches _full_model_transformed_data.stan:20-21
    // (baseline_visit_idx = curr_patient_visit_pos + n_patient_screening_visits[i] - 1),
    // which is the anchor t_patient_visit_idx uses for the SLD term. Using the first
    // visit here would mis-anchor g(tau) by n_screening-1 weeks and break the
    // shared-frame identity between the SLD and hazard terms. n_patient_screening_visits
    // is in scope from _visit_transformed_data.stan:16 (included before this module).
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    surrogate_bg_baseline_week[j] = t_patient_visits[vs + n_patient_screening_visits[p] - 1];
    // Visit-gating mask: mark each observed visit's calendar week.
    for (v in vs:ve) {
      int wk = t_patient_visits[v];
      if (wk >= 1 && wk <= max_all_t) surrogate_bg_visit_wk_01[j, wk] = 1;
    }
  }
}
