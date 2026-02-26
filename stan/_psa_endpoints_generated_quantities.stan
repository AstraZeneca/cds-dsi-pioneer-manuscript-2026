// PSA aliases for backward compatibility with R targets
vector[n_total_visits] rep_patient_log_psa = rep_patient_log_obs;
vector[n_total_visits] rep_mean_patient_log_psa = rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_psa = forecast_patient_log_obs;
vector[n_total_forecast_visits] forecast_mean_patient_log_psa = forecast_mean_patient_log_obs;

// PCWG3 classification for replicated and forecast PSA trajectories
// Categories: 1=Undetectable, 2=PSA50, 3=Stable, 4=PSA-PD, 5=Not Evaluable
array[sum(n_patient_visits)] int<lower=1, upper=5> rep_pcwg3 = rep_array(5, sum(n_patient_visits));
array[n_total_forecast_visits] int<lower=1, upper=4> forecast_pcwg3;

profile("psa_endpoints") {
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);
    int visit_size = get_pos_size(patient_visit_pos, i);
    int n_screening = n_patient_screening_visits[i];
    int n_treat = visit_size - n_screening;

    // Convert replicated log(PSA) means to absolute PSA for PCWG3 evaluation
    vector[visit_size] rep_psa_abs = exp(rep_mean_patient_log_obs[visit_start:visit_end]);

    if (n_treat > 0) {
      array[n_treat] int rep_cats = calculate_psa_category(
        rep_psa_abs,
        t_patient_visits[visit_start:visit_end],
        n_screening,
        psa_undetectable_threshold
      );
      rep_pcwg3[(visit_start + n_screening):visit_end] = rep_cats;
    }

    // Forecast PCWG3 categories
    int forecast_start, forecast_end;
    (forecast_start, forecast_end) = get_pos(forecast_visits_pos, i);
    int forecast_size = get_pos_size(forecast_visits_pos, i);

    if (forecast_size > 0) {
      // Concatenate last observed + forecast for nadir tracking continuity
      vector[1 + forecast_size] combined_psa;
      combined_psa[1] = rep_psa_abs[visit_size]; // last observed
      combined_psa[2:] = exp(forecast_mean_patient_log_obs[forecast_start:forecast_end]);

      array[1 + forecast_size] int combined_weeks;
      combined_weeks[1] = patient_last_obs_visit[i];
      combined_weeks[2:] = linspaced_int_array(
        forecast_size, patient_last_obs_visit[i] + 1, last_predict_visit);

      // Get nadir from observed post-treatment period for continuity
      real obs_nadir = n_treat > 0
        ? min(rep_psa_abs[(n_screening + 1):visit_size])
        : rep_psa_abs[visit_size];

      array[forecast_size] int forecast_cats = calculate_psa_category(
        combined_psa,
        combined_weeks,
        obs_nadir,
        1, // 1 "screening" visit = the anchor
        psa_undetectable_threshold
      );
      forecast_pcwg3[forecast_start:forecast_end] = forecast_cats;
    }
  }
}
