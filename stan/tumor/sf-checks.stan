// Data-only consistency checks executed once in transformed data.
// These avoid per-iteration overhead in model / generated quantities.
// Assumptions: helper functions (get_pos, get_visit_pos, get_int_sub_array) already declared.
// All referenced data objects are defined in previously included data/transformed data files.

// Loop over training patients only (indices align with arrays using global patient index i)
for (i in train_patients_pos:train_patients_end) {
    int train_idx = i - train_patients_pos + 1;

    // Visit partition and sizes
    int train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end;
    (train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end) = 
      get_visit_pos(train_patient_visit_pos, train_idx, n_patient_screening_visits[i]);
    int train_visit_size = get_pos_size(train_patient_visit_pos, train_idx);
    int n_screen = n_patient_screening_visits[i];
    int n_treat = train_visit_size - n_screen;

    if (!(n_screen >= 0)) fatal_error("sf-checks: negative screening count", i, n_screen);
  // Allow patients with zero post-screening (treatment) visits; just skip progression logic later.
  // Previously enforced n_treat > 0, which was too strict for patients dropping out after screening.
    if (!(train_treat_visit_start == train_visit_start + n_screen))
      fatal_error("sf-checks: treatment start offset mismatch", i, train_treat_visit_start, train_visit_start + n_screen);
    if (!(train_visit_end - train_visit_start + 1 == train_visit_size))
      fatal_error("sf-checks: visit size mismatch", i, train_visit_end - train_visit_start + 1, train_visit_size);

    // Extract all (screening + treatment) visit weeks for this patient
    array[train_visit_size] int curr_visits = get_int_sub_array(train_patient_visits, train_patient_visit_pos, train_idx);

    // Monotonicity of visit weeks
    for (k in 2:train_visit_size) {
      // Enforce strictly increasing week numbers (no duplicate weeks allowed)
      if (curr_visits[k] <= curr_visits[k-1])
        fatal_error("sf-checks: non-increasing observed visit weeks", i, curr_visits[k-1], curr_visits[k]);
    }

    // Last observed week consistency (should match patient_last_obs_visit)
    if (!(patient_last_obs_visit[i] == curr_visits[train_visit_size]))
      fatal_error("sf-checks: last observed week mismatch", i, patient_last_obs_visit[i], curr_visits[train_visit_size]);

    // Forecast sizing expectations
    if (!(n_patient_forecast_visits[i] >= 0))
      fatal_error("sf-checks: negative forecast visit count", i, n_patient_forecast_visits[i]);
    if (n_patient_forecast_visits[i] == 0) {
      // Nothing further to assert about future schedule for this patient
    } else {
      // Ensure at least one forecast step beyond last observed week
      if (!(last_predict_visit > patient_last_obs_visit[i]))
        fatal_error("sf-checks: last_predict_visit not after last observed", i, last_predict_visit, patient_last_obs_visit[i]);
    }

    // Cutoff bounds: each per-patient cutoff visit index must not exceed total potential (obs + future)
//   for (m in 1:n_mature_cutoffs_calendar_days) {
//     int cutoff_visit_idx = patient_relative_day_at_cutoff[i, m];
//     if (!(cutoff_visit_idx >= 0))
//       fatal_error("sf-checks: negative cutoff index", i, m, cutoff_visit_idx);
//   // Note: patient_relative_day_at_cutoff can legitimately exceed the number of
//   // (observed + planned forecast) visits; downstream code truncates via 'min'.
//   // We therefore do NOT assert an upper bound here.
//   }
}
