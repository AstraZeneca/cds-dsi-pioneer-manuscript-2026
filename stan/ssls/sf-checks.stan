// sf-checks.stan (relocated to ssls)
for (i in train_patients_pos:train_patients_end) {
    int train_idx = i - train_patients_pos + 1;
    int train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end;
    (train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end) = 
      get_visit_pos(train_patient_visit_pos, train_idx, n_patient_screening_visits[i]);
    int train_visit_size = get_pos_size(train_patient_visit_pos, train_idx);
    int n_screen = n_patient_screening_visits[i];
    int n_treat = train_visit_size - n_screen;

    if (!(n_screen >= 0)) fatal_error("sf-checks: negative screening count", i, n_screen);
    if (!(train_treat_visit_start == train_visit_start + n_screen))
      fatal_error("sf-checks: treatment start offset mismatch", i, train_treat_visit_start, train_visit_start + n_screen);
    if (!(train_visit_end - train_visit_start + 1 == train_visit_size))
      fatal_error("sf-checks: visit size mismatch", i, train_visit_end - train_visit_start + 1, train_visit_size);

    array[train_visit_size] int curr_visits = get_int_sub_array(train_patient_visits, train_patient_visit_pos, train_idx);

    for (k in 2:train_visit_size) {
      if (curr_visits[k] <= curr_visits[k-1])
        fatal_error("sf-checks: non-increasing observed visit weeks", i, curr_visits[k-1], curr_visits[k]);
    }

    if (!(patient_last_obs_visit[i] == curr_visits[train_visit_size]))
      fatal_error("sf-checks: last observed week mismatch", i, patient_last_obs_visit[i], curr_visits[train_visit_size]);

    if (!(n_patient_forecast_visits[i] >= 0))
      fatal_error("sf-checks: negative forecast visit count", i, n_patient_forecast_visits[i]);
    if (n_patient_forecast_visits[i] == 0) {
    } else {
      if (!(last_predict_visit > patient_last_obs_visit[i]))
        fatal_error("sf-checks: last_predict_visit not after last observed", i, last_predict_visit, patient_last_obs_visit[i]);
    }
}
