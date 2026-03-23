// ============================================================================
// VISITS MODULE TRANSFORMED DATA
// ============================================================================
// Two forecast grids:
//   1. Weekly dynamics grid — for multistate hazards, KM curves, visualization
//   2. Observation grid — for noisy burden measures and RECIST at assessment intervals

// --- Observation grid (assessment-visit subsampled from weekly grid) ---

array[n_patients] int<lower=0> n_patient_forecast_obs_visits;
for (i in 1:n_patients) {
  // ceil(n / interval) = (n - 1 + interval) %/% interval for n > 0
  n_patient_forecast_obs_visits[i] = n_patient_forecast_visits[i] > 0
    ? (n_patient_forecast_visits[i] - 1 + forecast_observation_interval) %/% forecast_observation_interval
    : 0;
}
int n_total_forecast_obs_visits = sum(n_patient_forecast_obs_visits);
array[n_patients + 1] int<lower=1> forecast_obs_visits_pos = create_pos(n_patient_forecast_obs_visits);

print("forecast_observation_interval = ", forecast_observation_interval);
print("n_total_forecast_obs_visits = ", n_total_forecast_obs_visits);
