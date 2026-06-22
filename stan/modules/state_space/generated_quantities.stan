// Back-transformed coefficients from QR space to original space
#include "modules/tr/generated_quantities.stan"
#include "modules/frac/generated_quantities.stan"
#include "modules/init/generated_quantities.stan"

// Biomarker-agnostic trajectory generation
// Contract: requires baseline_obs_per_patient[n_patients] (transformed data)
//           and measure_sd_obs (generated quantities scalar) to be defined

matrix[n_total_forecast_visits, 2] forecast_patient_states;
vector[n_total_visits] rep_patient_log_obs, rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_obs, forecast_mean_patient_log_obs;

// Per-patient static log-level for the combine function.
// Flag off OR no static patients -> -inf (static compartment absent).
// init_log_static_patient is indexed forecast-locally (size n_forecast_patients);
// map to the unified patient index p via forecast_patient_idx, matching states_full_grid.
vector[n_patients] static_log_level_per_patient = rep_vector(negative_infinity(), n_patients);
if (enable_static_init) {
  for (j in 1:n_forecast_patients) {
    static_log_level_per_patient[forecast_patient_idx[j]] = init_log_static_patient[j];
  }
}

profile("gen_quant_trajectories") {
  if (enable_patient_process_noise_tr) {
    // Process noise ON: Use states_full_grid (dense grid computed in transformed_parameters)
    (forecast_patient_states, rep_patient_log_obs, rep_mean_patient_log_obs,
     forecast_patient_log_obs, forecast_mean_patient_log_obs) =
      generate_all_patients_states_with_means_rng(
        states_full_grid,
        forecast_patient_idx,
        patient_visit_pos,
        patient_visit_m1_pos,
        forecast_visits_pos,
        patient_last_obs_visit,
        last_predict_visit,
        t_patient_visits,
        t_patient_visit_idx,
        baseline_obs_per_patient,
        static_log_level_per_patient,
        measure_sd_obs,
        n_patient_screening_visits
      );
  } else {
    // Process noise OFF: Compute states on-the-fly using constant rates (original approach)
    // This avoids needing states_full_grid which is not computed when process noise is off
    // Dual indexing: data arrays use unified patient_visit_pos[p]; states use forecast_visit_pos[j].
    for (j in 1:n_forecast_patients) {
      int p = forecast_patient_idx[j];  // unified patient index

      // Data positions (unified — index into visit-flat data arrays)
      int data_start, data_end;
      (data_start, data_end) = get_pos(patient_visit_pos, p);
      int visit_size = data_end - data_start + 1;

      // State positions (forecast-local — index into states[n_forecast_visits, 2])
      int state_start, state_end;
      (state_start, state_end) = get_pos(forecast_visit_pos, j);

      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, p);
      int forecast_size = forecast_visit_end - forecast_visit_start + 1;

      // Build forecast time WITH anchor (duplicate last observed week as element 1)
      array[forecast_size + 1] real forecast_time = linspaced_array(
        forecast_size + 1,
        patient_last_obs_visit[p],
        last_predict_visit);

      // Per-step Gompertz decay factor on the GROWTH rate only. The warp clock is
      // anchored at the patient's BASELINE week — the same origin the in-sample
      // states use in transformed_parameters — so the forecast CONTINUES the decay
      // the trajectory had already reached at the cutoff instead of re-accelerating.
      // kappa==0 (gr_decay off) => growth_warp(t,0)=t => factor==1 => identical to
      // the un-attenuated path, so the decay variant degrades gracefully.
      int baseline_visit_idx = data_start + n_patient_screening_visits[p] - 1;
      real baseline_week = t_patient_visits[baseline_visit_idx];
      real kappa_j = gr_decay_kappa[j];  // forecast-local index matches j
      vector[forecast_size + 1] forecast_tv_factor;
      for (t in 1:(forecast_size + 1)) {
        if (t == 1) {
          forecast_tv_factor[t] = 1.0;
        } else {
          real e_hi = forecast_time[t]     - baseline_week;
          real e_lo = forecast_time[t - 1] - baseline_week;
          real dphi = growth_warp(e_hi, kappa_j) - growth_warp(e_lo, kappa_j);
          real dt   = forecast_time[t] - forecast_time[t - 1];
          forecast_tv_factor[t] = dt > 0 ? dphi / dt : 1.0;
        }
      }

      // Generate states on-the-fly (not from grid). Uses the _decay RNG so the
      // Gompertz attenuation reaches the forecast — the non-decay variant grew the
      // growth arm at the FULL rate forever, producing implausibly high SLD.
      matrix[forecast_size, 2] temp_forecast_patient_states;
      vector[visit_size] temp_rep_patient_log_obs;
      vector[visit_size] temp_rep_mean_patient_log_obs;
      vector[forecast_size] temp_forecast_patient_log_obs;
      vector[forecast_size] temp_forecast_mean_patient_log_obs;
      matrix[visit_size - 1, 2] temp_obs_process_noise;  // Unused but required by function

      (temp_forecast_patient_states, temp_rep_patient_log_obs, temp_rep_mean_patient_log_obs,
       temp_forecast_patient_log_obs, temp_forecast_mean_patient_log_obs, temp_obs_process_noise) =
        generate_patient_states_with_means_decay_rng(
          states[state_start:state_end],    // forecast-local positions
          forecast_time,
          patient_log_decrease_rate[j, 1],  // forecast-local j
          patient_log_growth_rate[j, 1],    // forecast-local j
          baseline_obs_per_patient[p],      // unified
          static_log_level_per_patient[p],  // unified static log-level
          forecast_tv_factor,               // Gompertz decay weight on growth rate
          rep_matrix(0.0, forecast_size, 2),  // No forecast process noise
          measure_sd_obs
        );

      // Store results at unified data positions
      forecast_patient_states[forecast_visit_start:forecast_visit_end] = temp_forecast_patient_states;
      rep_patient_log_obs[data_start:data_end] = temp_rep_patient_log_obs;
      rep_mean_patient_log_obs[data_start:data_end] = temp_rep_mean_patient_log_obs;
      forecast_patient_log_obs[forecast_visit_start:forecast_visit_end] = temp_forecast_patient_log_obs;
      forecast_mean_patient_log_obs[forecast_visit_start:forecast_visit_end] = temp_forecast_mean_patient_log_obs;
    }
  }
}
