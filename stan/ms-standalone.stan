// ============================================================================
// SHARED STANDALONE MULTISTATE MODEL
// ============================================================================
// Illness-death model with GP baseline hazards and optional time-invariant
// covariates. No biomarker state-space dynamics (SLD, PSA, etc.).
//
// Used by both SCLC (tumor/SLD) and Pioneer (PSA) pipelines —
// the multistate structure is identical; only the data source differs.
//
// Transitions:
//   0→1: Progression (RECIST PD / PSA-PD)
//   0→2: On-trial death without progression
//   1→2: Post-progression death (semi-Markov)
//   0→3: Dropout
//   3→2: Off-trial death

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
}

data {
  #include "_hierarchy_data.stan"

  // Fit control — 0 for prior predictive, 1 for posterior
  int<lower=0, upper=1> fit_multistate_data;

  // Time grid length for hazard integration and KM/CIF computation
  int<lower=1> max_all_t;

  #include "_visit_data.stan"

  #include "modules/visits/data.stan"        // forecast_observation_interval
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/endpoints/data.stan"     // n_pfs_quantiles, pfs_quantiles, cond_group, etc.
  #include "modules/propensity/flags.stan"
  #include "modules/propensity/hyperparams.stan"
}

transformed data {
  // Standalone: all patients are forecast (no background split)
  int n_forecast_patients = n_patients;

  #include "_hierarchy_transformed_data.stan"
  #include "_forecast_routing_transformed_data.stan"

  array[n_forecast_patients] int forecast_patient_idx = linspaced_int_array(n_forecast_patients, 1, n_forecast_patients);
  array[n_trials + 1] int forecast_trial_patient_pos = trial_patient_pos;

  // Visit infrastructure: patient_visit_pos, n_total_visits, n_patient_screening_visits,
  // last_predict_visit, patient_last_obs_visit, forecast_visits_pos
  // max_all_t is a data input in standalone (not derived inline)
  #include "_visit_transformed_data.stan"
  #include "modules/visits/transformed_data.stan"

  // Time grid for multistate GP
  array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

  #include "_qr_decomposition.stan"
  #include "modules/multistate/transformed_data.stan"
  #include "modules/endpoints/transformed_data.stan"
  #include "modules/propensity/transformed_data.stan"
}

parameters {
  #include "modules/multistate/parameters.stan"
  #include "modules/propensity/parameters.stan"
}

transformed parameters {
  // No time-varying covariates in standalone — zero-sized placeholder to
  // satisfy the compiler's dimension check inside multistate/transformed_parameters.stan.
  array[enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 ? n_time_varying_covar : 0]
    matrix[n_patients, max_all_t] ms_time_varying_covar_01;

  #include "modules/multistate/transformed_parameters.stan"
  #include "modules/propensity/transformed_parameters.stan"
}

model {
  #include "modules/propensity/priors.stan"
  #include "modules/multistate/priors.stan"

  #include "modules/propensity/likelihood.stan"
  #include "modules/multistate/likelihood.stan"
}

generated quantities {
  #include "tumor/_ms_standalone_generated_quantities.stan"
}
