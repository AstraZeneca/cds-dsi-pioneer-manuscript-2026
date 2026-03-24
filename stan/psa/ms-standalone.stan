// ============================================================================
// PIONEER STANDALONE MULTISTATE MODEL
// ============================================================================
// Full model minus PSA state-space: multistate illness-death model only.
// Fits the survival structure (0→1, 0→2, 1→2, 0→3, 3→2) independently
// of PSA dynamics, using the same data and hierarchy as pioneer.stan.
//
// Uses _base_hierarchy_{data,transformed_data}.stan instead of the full
// _base_{data,transformed_data}.stan because the latter assumes all modules
// (tr, frac, init, state_space) are present. The few extra fields those
// files compute that are needed by shared modules are declared explicitly.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}

data {
  #include "_base_hierarchy_data.stan"
  #include "modules/visits/data.stan"
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/state_space/data.stan"

  // Visit schedule — needed by multistate likelihood (0→3 IC gap, visit-gated
  // dropout hazard) and GQ. Passed from the full PSA stan data.
  array[n_patients] int<lower=0> n_patient_visits;
  array[sum(n_patient_visits)] int t_patient_visits;  // includes pre-baseline (negative) weeks

  // max_all_t: passed as extend_max_all_t from the full PSA stan data
  int<lower=1> max_all_t;

  int<lower=0, upper=1> fit_multistate_data;
}

transformed data {
  // patient_visit_pos — computed from n_patient_visits (mirrors _base_transformed_data.stan)
  array[n_patients + 1] int<lower=1> patient_visit_pos = create_pos(n_patient_visits);

  // Standalone: all patients are HMC (no Laplace). Define before the includes
  // that compute HMC-dependent group sizes and position arrays.
  int n_hmc_patients = n_patients;

  #include "_base_hierarchy_transformed_data.stan"
  #include "_hmc_routing_transformed_data.stan"

  // hmc_patient_idx: trivial identity mapping (all patients are HMC)
  array[n_hmc_patients] int hmc_patient_idx = linspaced_int_array(n_hmc_patients, 1, n_hmc_patients);
  array[n_trials + 1] int hmc_trial_patient_pos = trial_patient_pos;

  // Time grid for multistate GP (in the full model from psa/transformed_data.stan)
  array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

  #include "_qr_decomposition.stan"
  #include "modules/multistate/transformed_data.stan"
  #include "modules/endpoints/transformed_data.stan"
}

parameters {
  #include "modules/multistate/parameters.stan"
}

transformed parameters {
  // No time-varying covariates: PSA-derived and only in the full joint model.
  // Declare zero-sized placeholder to satisfy the compiler's dimension check.
  array[enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 ? n_time_varying_covar : 0]
    matrix[n_patients, max_all_t] ms_time_varying_covar_01;

  #include "modules/multistate/transformed_parameters.stan"
}

model {
  #include "modules/multistate/priors.stan"

  if (fit_multistate_data) {
    profile("multistate loglik") {
      ms_final_state[hmc_patient_idx] ~ multistate(
        ones_vector(n_hmc_patients),   // no propensity weighting in standalone
        enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
        enable_ms_03, enable_ms_32,
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
        log_cond_surv_32
      );
    }
  }
}

generated quantities {
  #include "tumor/_ms_standalone_generated_quantities.stan"
}
