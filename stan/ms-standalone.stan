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
  #include "_base_hierarchy_data.stan"

  // Fit control — 0 for prior predictive, 1 for posterior
  int<lower=0, upper=1> fit_multistate_data;

  // Time grid length for hazard integration and KM/CIF computation
  int<lower=1> max_all_t;

  // Assessment visit schedule — needed for visit-gated 0→3 dropout hazard
  // and IC gap computation. Pass n_patient_visits (per-patient count);
  // patient_visit_pos is derived in transformed data.
  array[n_patients] int<lower=0> n_patient_visits;
  array[sum(n_patient_visits)] int t_patient_visits;  // includes pre-baseline weeks

  #include "modules/visits/data.stan"        // forecast_observation_interval
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/endpoints/data.stan"     // n_pfs_quantiles, pfs_quantiles, cond_group, etc.
  #include "modules/propensity/flags.stan"
  #include "modules/propensity/hyperparams.stan"
}

transformed data {
  // patient_visit_pos derived from n_patient_visits
  array[n_patients + 1] int<lower=1> patient_visit_pos = create_pos(n_patient_visits);

  // Standalone: all patients are forecast (no background split)
  int n_forecast_patients = n_patients;

  #include "_base_hierarchy_transformed_data.stan"
  #include "_forecast_routing_transformed_data.stan"

  array[n_forecast_patients] int forecast_patient_idx = linspaced_int_array(n_forecast_patients, 1, n_forecast_patients);
  array[n_trials + 1] int forecast_trial_patient_pos = trial_patient_pos;

  // Total visit count (needed for dummy biomarker arrays in GQ)
  int n_total_visits = sum(n_patient_visits);

  // Screening visit count per patient — visits with t <= 0 (mirrors _base_transformed_data.stan)
  array[n_patients] int<lower=0> n_patient_screening_visits = zeros_int_array(n_patients);
  for (i in 1:n_patients) {
    int v_start = patient_visit_pos[i];
    int v_end = patient_visit_pos[i + 1] - 1;
    for (v in v_start:v_end) {
      if (t_patient_visits[v] <= 0) n_patient_screening_visits[i] += 1;
    }
  }

  // Forecast visit infrastructure (mirrors _base_transformed_data.stan)
  int<lower=1> last_predict_visit = max_all_t;
  array[n_patients] int<lower=1> patient_last_obs_visit = get_max_pos(t_patient_visits, patient_visit_pos);
  array[n_patients] int<lower=0> n_patient_forecast_visits;
  for (i in 1:n_patients) {
    n_patient_forecast_visits[i] = last_predict_visit - patient_last_obs_visit[i];
  }
  array[n_patients + 1] int<lower=1> forecast_visits_pos = create_pos(n_patient_forecast_visits);
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
