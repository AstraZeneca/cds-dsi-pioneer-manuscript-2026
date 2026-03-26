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

  // Standalone: all patients are HMC (no Laplace split)
  int n_hmc_patients = n_patients;

  #include "_base_hierarchy_transformed_data.stan"
  #include "_hmc_routing_transformed_data.stan"

  array[n_hmc_patients] int hmc_patient_idx = linspaced_int_array(n_hmc_patients, 1, n_hmc_patients);
  array[n_trials + 1] int hmc_trial_patient_pos = trial_patient_pos;

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
  #include "modules/multistate/priors.stan"
  #include "modules/propensity/priors.stan"

  if (fit_multistate_data) {
    profile("multistate loglik") {
      target += reduce_sum(
          multistate_partial_sum,
          ms_final_state[hmc_patient_idx],  // auto-sliced
          1,                                 // grainsize=1: TBB auto-balances
          likelihood_weight[hmc_patient_idx],
          ms_time_01[hmc_patient_idx], ms_time_02[hmc_patient_idx], ms_time_12[hmc_patient_idx],
          ms_time_03[hmc_patient_idx], ms_time_32[hmc_patient_idx],
          ms_censored_01[hmc_patient_idx], ms_censored_02[hmc_patient_idx], ms_censored_12[hmc_patient_idx],
          ms_censored_32[hmc_patient_idx],
          ms_prog_deterministic[hmc_patient_idx],
          ms_ic_gap_01[hmc_patient_idx],
          t_patient_visits,
          patient_visit_pos,
          log_cond_surv_01,
          log_cond_surv_02,
          log_cond_surv_12_s,
          log_cond_surv_12_t,
          log_cond_surv_03,
          log_cond_surv_32,
          enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
          enable_ms_03, enable_ms_32,
          enable_ms_visit_gated_01
        );
    }
  }
}

generated quantities {
  #include "tumor/_ms_standalone_generated_quantities.stan"
}
