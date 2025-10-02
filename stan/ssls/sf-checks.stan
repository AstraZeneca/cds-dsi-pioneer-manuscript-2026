// sf-checks.stan (relocated to ssls)
for (i in 1:n_patients) {
    int visit_start, screening_visit_end, treat_visit_start, visit_end;
    (visit_start, screening_visit_end, treat_visit_start, visit_end) = 
      get_visit_pos(patient_visit_pos, i, n_patient_screening_visits[i]);
    int visit_size = get_pos_size(patient_visit_pos, i);
    int n_screen = n_patient_screening_visits[i];
    int n_treat = visit_size - n_screen;

    if (!(n_screen >= 0)) fatal_error("sf-checks: negative screening count", i, n_screen);
    if (!(treat_visit_start == visit_start + n_screen))
      fatal_error("sf-checks: treatment start offset mismatch", i, treat_visit_start, visit_start + n_screen);
    if (!(visit_end - visit_start + 1 == visit_size))
      fatal_error("sf-checks: visit size mismatch", i, visit_end - visit_start + 1, visit_size);

    array[visit_size] int curr_visits = get_int_sub_array(t_patient_visits, patient_visit_pos, i);

    for (k in 2:visit_size) {
      if (curr_visits[k] <= curr_visits[k-1])
        fatal_error("sf-checks: non-increasing observed visit weeks", i, curr_visits[k-1], curr_visits[k]);
    }

    if (!(patient_last_obs_visit[i] == curr_visits[visit_size]))
      fatal_error("sf-checks: last observed week mismatch", i, patient_last_obs_visit[i], curr_visits[visit_size]);

    if (!(n_patient_forecast_visits[i] >= 0))
      fatal_error("sf-checks: negative forecast visit count", i, n_patient_forecast_visits[i]);
    if (n_patient_forecast_visits[i] == 0) {
    } else {
      if (!(last_predict_visit > patient_last_obs_visit[i]))
        fatal_error("sf-checks: last_predict_visit not after last observed", i, last_predict_visit, patient_last_obs_visit[i]);
    }
}
