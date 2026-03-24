// ============================================================================
// PSA/PCWG3 CATEGORIZATION
// ============================================================================
// Converts generic state-space trajectories to PSA-specific quantities and
// PCWG3 categories. Included at GQ scope — all variables here are saved.
//
// Requires (from state_space/generated_quantities.stan):
//   rep_patient_log_obs, rep_mean_patient_log_obs
//   forecast_patient_log_obs, forecast_mean_patient_log_obs
//
// Requires (from modules/visits/):
//   forecast_obs_visits_pos, forecast_observation_interval
//
// Requires (from psa/transformed_data.stan):
//   UNDETECTABLE, PSA50, STABLE, PSA_PD, NE

// --- PSA aliases (for R targets backward compatibility) ---
vector[n_total_visits] rep_patient_log_psa = rep_patient_log_obs;
vector[n_total_visits] rep_mean_patient_log_psa = rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_psa = forecast_patient_log_obs;
vector[n_total_forecast_visits] forecast_mean_patient_log_psa = forecast_mean_patient_log_obs;

// --- PCWG3 classification ---
// rep_pcwg3: categories at observed visits (from mean, for visualization).
// forecast_obs_pcwg3: categories at assessment visits (subsampled from noisy
//   weekly forecast internally, for endpoint computation).
array[n_total_visits] int<lower=UNDETECTABLE, upper=NE> rep_pcwg3;
array[n_total_forecast_obs_visits] int forecast_obs_pcwg3;

(rep_pcwg3, forecast_obs_pcwg3) = calculate_all_patients_pcwg3(
  rep_mean_patient_log_obs,
  forecast_patient_log_obs,
  patient_visit_pos,
  forecast_visits_pos,
  forecast_obs_visits_pos,
  n_patient_screening_visits,
  t_patient_visits_day,
  patient_last_obs_visit,
  forecast_observation_interval,
  last_predict_visit,
  psa_undetectable_threshold
);
