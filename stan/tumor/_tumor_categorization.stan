// ============================================================================
// TUMOR/SLD CATEGORIZATION
// ============================================================================
// Converts generic state-space trajectories to SLD-specific quantities and
// RECIST categories. Included at GQ scope — all variables here are saved.
//
// Requires (from state_space/generated_quantities.stan):
//   rep_patient_log_obs, rep_mean_patient_log_obs
//   forecast_patient_log_obs, forecast_mean_patient_log_obs
//
// Requires (from modules/visits/):
//   forecast_obs_visits_pos, forecast_observation_interval

// --- SLD aliases (for R targets backward compatibility) ---
vector[n_total_visits] rep_patient_log_sld = rep_patient_log_obs;
vector[n_total_visits] rep_mean_patient_log_sld = rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_sld = forecast_patient_log_obs;
vector[n_total_forecast_visits] forecast_mean_patient_log_sld = forecast_mean_patient_log_obs;

// --- RECIST classification ---
// rep_recist: categories at observed visits (for R backward compat).
// forecast_obs_recist: categories at assessment visits (subsampled from noisy
//   weekly forecast internally, for endpoint computation and visualization).
array[n_total_visits] int<lower = CR, upper = PD + 1> rep_recist;
array[n_total_forecast_obs_visits] int forecast_obs_recist;

(rep_recist, forecast_obs_recist) = calculate_all_patients_recist(
  rep_patient_log_sld,
  forecast_patient_log_sld,
  patient_visit_pos,
  forecast_visits_pos,
  forecast_obs_visits_pos,
  forecast_observation_interval,
  n_patient_screening_visits
);
