// Biomarker-agnostic trajectory generation
// Contract: requires baseline_obs_per_patient[n_patients] (transformed data)
//           and measure_sd_obs (generated quantities scalar) to be defined

matrix[n_total_forecast_visits, 2] forecast_patient_states;
vector[n_total_visits] rep_patient_log_obs, rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_obs, forecast_mean_patient_log_obs;

profile("gen_quant_trajectories") {
  if (enable_patient_process_noise_tr) {
    // Process noise ON: Use states_full_grid (dense grid computed in transformed_parameters)
    (forecast_patient_states, rep_patient_log_obs, rep_mean_patient_log_obs,
     forecast_patient_log_obs, forecast_mean_patient_log_obs) =
      generate_all_patients_states_with_means_rng(
        states_full_grid,
        patient_visit_pos,
        patient_visit_m1_pos,
        forecast_visits_pos,
        patient_last_obs_visit,
        last_predict_visit,
        t_patient_visits,
        t_patient_visit_idx,
        baseline_obs_per_patient,
        measure_sd_obs,
        n_patient_screening_visits
      );
  } else {
    // Process noise OFF: Compute states on-the-fly using constant rates (original approach)
    // This avoids needing states_full_grid which is not computed when process noise is off
    for (i in 1:n_patients) {
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, i);
      int visit_size = get_pos_size(patient_visit_pos, i);

      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = get_pos_size(forecast_visits_pos, i);

      // Build forecast time WITH anchor (duplicate last observed week as element 1)
      array[forecast_size + 1] real forecast_time = linspaced_array(
        forecast_size + 1,
        patient_last_obs_visit[i],
        last_predict_visit);

      // Generate states using constant rates (computed on-the-fly, not from grid)
      matrix[forecast_size, 2] temp_forecast_patient_states;
      vector[visit_size] temp_rep_patient_log_obs;
      vector[visit_size] temp_rep_mean_patient_log_obs;
      vector[forecast_size] temp_forecast_patient_log_obs;
      vector[forecast_size] temp_forecast_mean_patient_log_obs;
      matrix[visit_size - 1, 2] temp_obs_process_noise;  // Unused but required by function

      (temp_forecast_patient_states, temp_rep_patient_log_obs, temp_rep_mean_patient_log_obs,
       temp_forecast_patient_log_obs, temp_forecast_mean_patient_log_obs, temp_obs_process_noise) =
        generate_patient_states_with_means_rng(
          states[visit_start:visit_end],  // Use states computed in transformed_parameters
          forecast_time,
          patient_log_decrease_rate[i, 1],  // Scalar rate (constant)
          patient_log_growth_rate[i, 1],    // Scalar rate (constant)
          baseline_obs_per_patient[i],
          negative_infinity(),  // growth lag (disabled)
          1.0,                  // growth transition
          rep_matrix(0.0, forecast_size, 2),  // No forecast process noise
          measure_sd_obs
        );

      // Store results
      forecast_patient_states[forecast_visit_start:forecast_visit_end] = temp_forecast_patient_states;
      rep_patient_log_obs[visit_start:visit_end] = temp_rep_patient_log_obs;
      rep_mean_patient_log_obs[visit_start:visit_end] = temp_rep_mean_patient_log_obs;
      forecast_patient_log_obs[forecast_visit_start:forecast_visit_end] = temp_forecast_patient_log_obs;
      forecast_mean_patient_log_obs[forecast_visit_start:forecast_visit_end] = temp_forecast_mean_patient_log_obs;
    }
  }
}
