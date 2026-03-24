// Standalone Multistate Illness-Death Model
//
// A self-contained illness-death model with GP baseline hazards and optional
// time-invariant covariates. No tumor dynamics or time-varying covariates.
//
// This model reuses all multistate module includes from the full joint model
// (sf-ssm-log-space.stan) but excludes tumor regression, growth fraction,
// initial state, state-space, and time-varying covariate dependencies.
//
// Transitions:
//   0→1: Progression / PFS event
//   0→2: Death without progression
//   1→2: Post-progression death (semi-Markov, Markov, or extended)

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}

data {
  // =========================================================================
  // BASE HIERARCHY DATA (shared with full model)
  // =========================================================================
  #include "_base_hierarchy_data.stan"

  // Fit control — set to 0 for prior predictive, 1 for posterior
  int<lower=0, upper=1> fit_multistate_data;

  // Time grid length — in the joint model this is computed from visit times;
  // here it is passed directly since there are no tumor visits.
  int<lower=1> max_all_t;

  // Assessment visit arrays — needed for visit-conditioned 0→3 dropout hazard.
  // In the joint model these come from _base_data.stan/_base_transformed_data.stan;
  // here they are passed directly.
  int<lower=0> n_total_visits;
  array[n_total_visits] int<lower=1> t_patient_visits;     // flat visit-week array
  array[n_patients + 1] int<lower=1> patient_visit_pos;   // position array (n_patients+1)
  int<lower=1> forecast_observation_interval;              // weeks between forecast assessments (typically 6)

  // =========================================================================
  // MULTISTATE MODULE DATA
  // =========================================================================
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"

  // =========================================================================
  // ENDPOINT COMPUTATION DATA (shared with full model)
  // =========================================================================
  #include "modules/endpoints/data.stan"
}

transformed data {
  // =========================================================================
  // HIERARCHY VALIDATION + TRIAL POSITION ARRAYS (shared with full model)
  // =========================================================================
  // Standalone: all patients are HMC (no Laplace)
  int n_hmc_patients = n_patients;

  #include "_base_hierarchy_transformed_data.stan"
  #include "_hmc_routing_transformed_data.stan"

  array[n_hmc_patients] int hmc_patient_idx = linspaced_int_array(n_hmc_patients, 1, n_hmc_patients);
  array[n_trials + 1] int hmc_trial_patient_pos = trial_patient_pos;

  // =========================================================================
  // TIME GRID (needed by multistate GP; in the full model this lives in
  // modules/tumor/transformed_data.stan)
  // =========================================================================
  array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

  // =========================================================================
  // QR DECOMPOSITION (shared with full model)
  // =========================================================================
  #include "_qr_decomposition.stan"

  // =========================================================================
  // ENDPOINT POSITION ARRAYS (shared with full model)
  // =========================================================================
  #include "modules/endpoints/transformed_data.stan"

  // =========================================================================
  // MULTISTATE MODULE TRANSFORMED DATA
  // =========================================================================
  #include "modules/multistate/transformed_data.stan"
}

parameters {
  #include "modules/multistate/parameters.stan"
}

transformed parameters {
  // No _ms_time_varying_covar.stan — time-varying covariates disabled.
  // However, the multistate module references ms_time_varying_covar_01 inside
  // a runtime guard. Stan requires identifiers in scope at compile time, so we
  // declare a zero-sized placeholder. The guard
  //   `if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0)`
  // ensures the array is never actually accessed.
  array[enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 ? n_time_varying_covar : 0]
    matrix[n_patients, max_all_t] ms_time_varying_covar_01;

  #include "modules/multistate/transformed_parameters.stan"
}

model {
  #include "modules/multistate/priors.stan"

  // Multistate likelihood (skipped when fit_multistate_data = 0 for prior predictive)
  if (fit_multistate_data) {
    ms_final_state ~ multistate(
      ones_vector(n_patients),   // no propensity weighting in standalone
      enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
      enable_ms_03, enable_ms_32,
      ms_time_01, ms_time_02, ms_time_12,
      ms_time_03, ms_time_32,
      ms_censored_01, ms_censored_02, ms_censored_12,
      ms_censored_32,
      ms_prog_deterministic,
      ms_ic_gap_01,
      t_patient_visits,
      patient_visit_pos,
      log_cond_surv_01,
      log_cond_surv_02,
      log_cond_surv_12_s,
      log_cond_surv_12_t,
      log_cond_surv_03,
      log_cond_surv_32
    );
  }
}

generated quantities {
  #include "_ms_standalone_generated_quantities.stan"
}
