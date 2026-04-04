// ============================================================================
// PSA TRAJECTORY ALIASES
// ============================================================================
// Renames generic state-space outputs to PSA-specific names for R targets.
// Included at GQ scope — all variables here are saved to CSV.
//
// PCWG3 categorization (rep_pcwg3, forecast_obs_pcwg3) is computed inside
// the endpoint local block in pioneer.stan to avoid saving large arrays.
//
// Requires (from state_space/generated_quantities.stan):
//   rep_patient_log_obs, rep_mean_patient_log_obs
//   forecast_patient_log_obs, forecast_mean_patient_log_obs

// --- PSA aliases (for R targets backward compatibility) ---
vector[n_total_visits] rep_patient_log_psa = rep_patient_log_obs;
vector[n_total_visits] rep_mean_patient_log_psa = rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_psa = forecast_patient_log_obs;
vector[n_total_forecast_visits] forecast_mean_patient_log_psa = forecast_mean_patient_log_obs;
