// ============================================================================
// STANDALONE MULTISTATE MODEL — LEAVE-FUTURE-OUT CROSS-VALIDATION
// ============================================================================
// Same illness-death structure as ms-standalone.stan, but with temporal
// cross-validation: trains on data re-censored at a calendar cutoff,
// evaluates held-out event log-likelihoods in generated quantities.
//
// Used by the lfo() function in r/accuracy.R, which manages the outer
// refit loop and PSIS approximation.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
}

data {
  #include "_hierarchy_data.stan"

  int<lower=0, upper=1> fit_multistate_data;
  int<lower=1> max_all_t;

  #include "_visit_data.stan"

  // Calendar timing (needed for LFO cutoff computation)
  array[sum(n_patient_visits)] int t_patient_visits_day;
  array[n_patients] int<lower=1> calendar_day;

  #include "modules/visits/data.stan"
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/endpoints/data.stan"
  #include "modules/propensity/flags.stan"
  #include "modules/propensity/hyperparams.stan"

  // LFO cutoff schedule
  #include "modules/state_space/lfo_data.stan"
}

transformed data {
  int n_forecast_patients = n_patients;

  #include "_hierarchy_transformed_data.stan"
  #include "_forecast_routing_transformed_data.stan"

  array[n_forecast_patients] int forecast_patient_idx =
    linspaced_int_array(n_forecast_patients, 1, n_forecast_patients);
  array[n_trials + 1] int forecast_trial_patient_pos = trial_patient_pos;

  #include "_visit_transformed_data.stan"
  #include "modules/visits/transformed_data.stan"

  array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

  #include "_qr_decomposition.stan"
  // Standalone never uses inline burden (no state-space model)
  int ms_needs_inline_burden = 0;
  #include "modules/multistate/transformed_data.stan"
  #include "modules/endpoints/transformed_data.stan"
  #include "modules/propensity/transformed_data.stan"

  // LFO: re-censor events at cutoff, identify training/testing patients
  #include "_ms_standalone_lfo_transformed_data.stan"

  // LFO: compact cutoff-observed cohort mapping (for GQ endpoint aggregation)
  #include "_ms_standalone_lfo_compact_transformed_data.stan"
}

parameters {
  #include "modules/multistate/parameters.stan"
  #include "modules/propensity/parameters.stan"
}

transformed parameters {
  // No state-space model in standalone — always zero-sized placeholder.
  // Latent-PSA variants (visit_gated_latent_01=TRUE) have no trajectory source
  // here; the TP guard `(!latent || size > 0)` evaluates to FALSE, skipping
  // that block without NaN.  Non-latent visit-gated variants still use the
  // observed-PSA path because the same guard becomes `!FALSE || FALSE` = TRUE.
  array[0] matrix[n_patients, max_all_t] ms_time_varying_covar_01;

  #include "modules/multistate/transformed_parameters.stan"
  #include "modules/multistate/cond_surv_transform.stan"
  #include "modules/propensity/transformed_parameters.stan"
}

model {
  #include "modules/propensity/priors.stan"
  #include "modules/multistate/priors.stan"

  // Train on re-censored data. Patients whose last visit is after the cutoff
  // have their lfo_ms_* fields zeroed by recensor_ms_at_cutoff, so the
  // time_*[i] > 0 guards inside multistate_lpmf naturally exclude them —
  // no extra weight gate needed.
  if (fit_multistate_data) {
    profile("multistate loglik") {
      lfo_ms_final_state[forecast_patient_idx] ~ multistate(
          likelihood_weight[forecast_patient_idx],
          enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
          enable_ms_03, enable_ms_32,
          lfo_ms_time_01[forecast_patient_idx], lfo_ms_time_02[forecast_patient_idx],
          lfo_ms_time_12[forecast_patient_idx],
          lfo_ms_time_03[forecast_patient_idx], lfo_ms_time_32[forecast_patient_idx],
          lfo_ms_censored_01[forecast_patient_idx],
          lfo_ms_prog_deterministic[forecast_patient_idx],
          lfo_ms_ic_gap_01[forecast_patient_idx],
          t_patient_visits,
          patient_visit_pos,
          log_cond_surv_01,
          log_cond_surv_02,
          log_cond_surv_12_s,
          log_cond_surv_12_t,
          log_cond_surv_03,
          log_cond_surv_32,
          enable_ms_visit_gated_01
        );
    }
  }
}

generated quantities {
  // LFO: held-out log-likelihoods per transition for temporal cross-validation.
  #include "_ms_standalone_lfo_generated_quantities.stan"

  // Posterior-predictive OS KM from re-censored events — pre-cutoff events
  // are preserved, post-cutoff events are hidden so the KM represents the
  // model's out-of-sample forecast at this cutoff.
  #include "_ms_standalone_lfo_os_km_generated_quantities.stan"
}
