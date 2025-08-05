/**
 * Determine which patient visits occurred before a given calendar day cutoff
 *
 * IMPORTANT: The cutoff is defined as the first day/week that is out-of-sample and NOT trained on.
 * All data on or after the cutoff is excluded from training and considered for testing/forecasting only.
 *
 * This version works with patient visits rather than individual tumor measurements.
 * It finds the last visit where SLD was measured across all tumors before the cutoff.
 *
 * @param cutoff_calendar_day Calendar date to use as cutoff
 * @param patient_calendar_day Calendar date when each patient entered the study
 * @param t_patient_visits Week numbers for patient visits (ragged array)
 * @param t_patient_visits_day Day numbers for patient visits (ragged array)
 * @param patient_visit_pos Position array defining boundaries for each patient's visits
 * @return Tuple of (last_visit_day, last_visit_week, last_visit_calendar_day, last_visit_week_observed, cutoff_last_visit_idx)
 */
tuple(array[] int, array[] int, array[] int, array[] int, array[] int) cutoff_visits(
  int cutoff_calendar_day, 
  array[] int patient_calendar_day, 
  array[] int t_patient_visits,           // Week numbers for all patient visits
  array[] int t_patient_visits_day,       // Day numbers for all patient visits
  array[] int patient_visit_pos           // Position array for patient visits
) {
  int n_patients = size(patient_calendar_day);
  
  // 0 is the default sentinel value if last visit is before cutoff
  array[n_patients] int last_visit_day = zeros_int_array(n_patients);
  array[n_patients] int last_visit_week = zeros_int_array(n_patients), 
                        last_visit_week_observed = ones_int_array(n_patients); 
  array[n_patients] int last_visit_calendar_day;
  
  // Precompute for each patient the last visit index before or at cutoff_last_visit_week
  array[n_patients] int cutoff_last_visit_idx = zeros_int_array(n_patients);
  
  for (i in 1:n_patients) {
    // Get start and end positions for this patient's visits
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);
    int n_patient_visits = visit_end - visit_start + 1;
    if (n_patient_visits <= 0) {
      last_visit_day[i] = 0;
      last_visit_week[i] = 0;
      last_visit_calendar_day[i] = 0;
      cutoff_last_visit_idx[i] = 0;
      continue;
    }
    // Convert global cutoff to patient-specific study day
    int patient_cutoff_study_day = calendar_date_to_study_date(patient_calendar_day[i], cutoff_calendar_day);
    // Extract this patient's visit data
    array[n_patient_visits] int patient_t_visits_week = t_patient_visits[visit_start:visit_end];
    array[n_patient_visits] int patient_t_visits_day = t_patient_visits_day[visit_start:visit_end];
    // Sort visits by day to find chronological order
    array[n_patient_visits] int patient_visit_sort_idx = sort_indices_asc(patient_t_visits_day);
    // Find the last visit on or before the cutoff
    int t_idx = 0;
    while (t_idx < n_patient_visits && 
           patient_t_visits_day[patient_visit_sort_idx[t_idx + 1]] <= patient_cutoff_study_day) {
      t_idx += 1;
    }
    // If at least one visit occurred before cutoff, record it
    if (t_idx > 0) {  
      last_visit_day[i] = patient_t_visits_day[patient_visit_sort_idx[t_idx]]; 
      last_visit_week[i] = patient_t_visits_week[patient_visit_sort_idx[t_idx]]; 
    }

    last_visit_week_observed[i] = last_visit_week[i] > 0;

    // Calculate the calendar date of this patient's final visit (regardless of cutoff)
    last_visit_calendar_day[i] = patient_calendar_day[i] + 
                                 patient_t_visits_day[patient_visit_sort_idx[n_patient_visits]] - 1;

    // Calculate cutoff_last_visit_idx - the index of the last visit before or at cutoff_last_visit_week
    if (last_visit_week[i] > 0) {
      int j = visit_start;
      while (j <= visit_end && t_patient_visits[j] <= last_visit_week[i]) {
        cutoff_last_visit_idx[i] = j;
        j += 1;
      }
    } else {
      cutoff_last_visit_idx[i] = 0;
    }
  }
  
  return (last_visit_day, last_visit_week, last_visit_calendar_day, last_visit_week_observed, cutoff_last_visit_idx);
}

tuple(array[] int, array[] int, array[] int) fine_cutoff_visits(
  int cutoff_calendar_day, array[] int patient_calendar_day, array[] int t_measure, array[] int t_day_measure, array[] int patient_tumor_measure_pos
) {
  int n_patients = size(patient_calendar_day);
  
  // 0 is the default sentinel value if last visit is negative 
  array[n_patients] int last_visit_day = zeros_int_array(n_patients), last_visit_week = zeros_int_array(n_patients); 
  
  array[n_patients] int last_visit_calendar_day;
  
  for (i in 1:n_patients) {
    int t_measure_pos = patient_tumor_measure_pos[i];
    int t_measure_end = patient_tumor_measure_pos[i + 1] - 1;
    int n_patient_measures = t_measure_end - t_measure_pos + 1;
    if (n_patient_measures <= 0) {
      last_visit_day[i] = 0;
      last_visit_week[i] = 0;
      last_visit_calendar_day[i] = 0;
      continue;
    }
    // Defensive: only access arrays if n_patient_measures > 0
    int patient_cutoff_study_day = calendar_date_to_study_date(patient_calendar_day[i], cutoff_calendar_day);
    array[n_patient_measures] int patient_t_measure;
    array[n_patient_measures] int patient_t_day_measure;
    for (j in 1:n_patient_measures) {
      patient_t_measure[j] = t_measure[t_measure_pos + j - 1];
      patient_t_day_measure[j] = t_day_measure[t_measure_pos + j - 1];
    }
    array[n_patient_measures] int patient_measure_t_sort_idx = sort_indices_asc(patient_t_day_measure);
    int t_idx = 1;
    while (t_idx <= n_patient_measures && patient_t_day_measure[patient_measure_t_sort_idx[t_idx]] <= patient_cutoff_study_day) {
      t_idx += 1;
    }
    if (t_idx > 1) {
      last_visit_day[i] = patient_t_day_measure[patient_measure_t_sort_idx[t_idx - 1]];
      last_visit_week[i] = patient_t_measure[patient_measure_t_sort_idx[t_idx - 1]];
    }
    last_visit_calendar_day[i] = patient_calendar_day[i] + patient_t_day_measure[patient_measure_t_sort_idx[n_patient_measures]] - 1;
  }
  
  return (last_visit_day, last_visit_week, last_visit_calendar_day);
} 

/**
 * Get indices of out-of-sample patients for leave-future-out cross-validation
 * 
 * This function identifies which patients should be included in the testing set for each
 * cutoff date. It works with a pre-sorted array of patient last visit dates to efficiently
 * find the first patient with data after each cutoff. All subsequent patients in the 
 * sorted array will also have data after that cutoff and are included in testing.
 * 
 * The function enables temporal cross-validation where we train on data up to a cutoff
 * date and test on future observations. This mimics real-world scenarios where we want
 * to predict future outcomes based on historical data.
 * 
 * @param sorted_last_visit_calendar_day Array of calendar dates for each patient's last 
 *                                       visit, pre-sorted in ascending order
 * @param cutoff_calendar_day Array of calendar dates to use as temporal cutoffs for 
 *                            cross-validation splits
 * 
 * @return Array where element i contains the index (in sorted_last_visit_calendar_day) 
 *         of the first patient to include in the testing set for cutoff i. All patients
 *         from this index to the end of the sorted array have visits after cutoff i.
 * 
 * @example
 * sorted_last_visit_calendar_day = [100, 120, 150, 180, 200, 220]  // 6 patients
 * cutoff_calendar_day = [110, 160]  // 2 cutoffs
 * 
 * Returns: [2, 4]
 * - Cutoff 110: First patient with last visit > 110 is at index 2 (day 120)
 *   Testing set for cutoff 1: patients at indices [2, 3, 4, 5, 6]
 * - Cutoff 160: First patient with last visit > 160 is at index 4 (day 180)  
 *   Testing set for cutoff 2: patients at indices [4, 5, 6]
 * 
 * @note If no patients have visits after a cutoff, a warning is printed and the 
 *       corresponding index may not be set properly.
 */
array[] int get_oos_patients_idx(
  array[] int sorted_last_visit_calendar_day, 
  array[] int cutoff_calendar_day
) {
  int n_cutoffs = size(cutoff_calendar_day);
  int n_patients = size(sorted_last_visit_calendar_day);

  // Initialize output array with zeros
  array[n_cutoffs] int patient_idx = rep_array(0, n_cutoffs);
  
  int n_remaining_testing_patients = n_patients;
  int curr_patient_idx = 1;      // Current position in sorted patient array
  int patient_idx_pos = 1;       // Current position in output array

  // Iterate through sorted patients to find first patient after each cutoff
  while (curr_patient_idx <= n_patients && patient_idx_pos <= n_cutoffs) {
    // Skip patients whose last visit is on or before the current cutoff
    while (curr_patient_idx <= n_patients && 
           sorted_last_visit_calendar_day[curr_patient_idx] <= cutoff_calendar_day[patient_idx_pos]) {
      curr_patient_idx += 1;
    }

    // If we found a patient with visits after this cutoff, record the index
    if (curr_patient_idx <= n_patients) {
      patient_idx[patient_idx_pos] = curr_patient_idx;
      patient_idx_pos += 1;
    } else {
      // No more patients have visits after this cutoff
      print("Warning: no patients have any visits after cutoff ", patient_idx_pos);
    }
  }

  return patient_idx;
}

/**
 * Calculate testing visit week boundaries for leave-future-out cross-validation
 * 
 * This function determines the time windows (in actual week numbers) for evaluating
 * model predictions in temporal cross-validation. It supports nested evaluation windows
 * where you can train up to one cutoff and test between that cutoff and a later one.
 * 
 * The function handles the complexity of patient-specific study calendars by converting
 * global cutoff dates to patient-specific study days, then finding which visits fall
 * within the testing windows.
 * 
 * @param oos_patient_idx Array of starting indices (in sorted patient list) for 
 *                        out-of-sample patients at each cutoff
 * @param last_visit_calendar_day_sort_idx Indices that sort patients by their last 
 *                                         visit calendar date (ascending)
 * @param cutoff_calendar_day Array of calendar dates defining temporal cutoffs
 * @param patient_calendar_day Calendar date when each patient entered the study
 * @param t_measure Week numbers for all tumor measurements across all patients
 * @param t_day_measure Day numbers for all tumor measurements across all patients
 * @param patient_tumor_measure_pos Position array defining where each patient's 
 *                                  measurements start/end in the measurement arrays
 * 
 * @return Tuple of two arrays:
 *   1. first_testing_visit_week[n_futures, n_patients]: Week number of each patient's
 *      first visit after cutoff n. Zero if no visits after cutoff.
 *   2. last_testing_visit_week[n_futures, n_futures, n_patients]: Week number of each
 *      patient's last visit on or before cutoff m when training up to cutoff n.
 *      Initialized to min(t_measure) and updated only for valid (n,m) pairs where m > n.
 * 
 * @example
 * Patient has visits at weeks [2, 5, 8, 12, 15] on days [14, 35, 56, 84, 105]
 * Patient's study start: day 100
 * Cutoff 1: calendar day 140 (= study day 40)
 * Cutoff 2: calendar day 190 (= study day 90)
 * 
 * Results:
 * - first_testing_visit_week[cutoff1, patient] = 8 (first visit after day 40)
 * - last_testing_visit_week[cutoff1, cutoff2, patient] = 12 (last visit before day 90)
 * - Testing window for (cutoff1→cutoff2): weeks 8-12
 * 
 * @note 
 * - Returns actual week numbers, not indices
 * - Only processes patients identified as out-of-sample by oos_patient_idx
 * - For the last cutoff (m = n_futures), uses all remaining visits (up to max_all_t)
 * - Returns 0 as sentinel value when no testing visits found after cutoff
 */
tuple(array[,] int, array[,,] int, array[,] int, array[,,] int) get_testing_visit_week_bounds(
  array[] int oos_patient_idx, 
  array[] int last_visit_calendar_day_sort_idx,
  array[] int cutoff_calendar_day, 
  array[] int patient_calendar_day,
  array[] int t_patient_visits_week, 
  array[] int t_patient_visits_day, 
  array[] int patient_visit_pos
) {
  int n_patients = size(patient_calendar_day);
  int n_futures = size(oos_patient_idx);
  int min_all_t = min(t_patient_visits_week);
  
  // Initialize output arrays
  array[n_futures, n_patients] int first_testing_visit_week = rep_array(0, n_futures, n_patients);
  array[n_futures, n_futures, n_patients] int last_testing_visit_week = rep_array(min_all_t, n_futures, n_futures, n_patients);
  array[n_futures, n_patients] int testing_start_idx = rep_array(0, n_futures, n_patients);
  array[n_futures, n_futures, n_patients] int testing_end_idx = rep_array(0, n_futures, n_futures, n_patients);
  
  // Process each cutoff
  for (n in 1:n_futures) {
    // Get patients who have visits after cutoff n
    int n_curr_patients = n_patients - oos_patient_idx[n] + 1;
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[oos_patient_idx[n]:];
    
    // Convert cutoff n to patient-specific study days
    array[n_curr_patients] int lower_cutoff_visit_days = calendar_date_to_study_date(patient_calendar_day[curr_patients], cutoff_calendar_day[n]);

    // if (n == 1) {
    //   print("patient_calendar_day[curr_patients] = ", patient_calendar_day[curr_patients], ", cutoff_calendar_day[n] = ", cutoff_calendar_day[n]);
    // }
    
    // Find first testing visit for each patient after cutoff n
    for (i_idx in 1:n_curr_patients) {
      int i = curr_patients[i_idx];  // Actual patient ID
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, i);
      int n_patient_visits = visit_end - visit_start + 1;
      
      // Extract and sort patient's visits by day
      array[n_patient_visits] int patient_t_visits_week = t_patient_visits_week[visit_start:visit_end];
      array[n_patient_visits] int patient_t_visits_day = t_patient_visits_day[visit_start:visit_end];

      assert_strict_ascending(patient_t_visits_week);
      assert_strict_ascending(patient_t_visits_day);
      
      // Find first visit after cutoff
      int t_idx = 1;
      while (t_idx <= n_patient_visits && 
             patient_t_visits_day[t_idx] <= max(0, lower_cutoff_visit_days[i_idx])) {
        t_idx += 1;
      }
      
      if (t_idx <= n_patient_visits) {
        first_testing_visit_week[n, i] = patient_t_visits_week[t_idx];
        testing_start_idx[n, i] = visit_start + t_idx - 1;
      } 

      // if (i_idx <= 2 && n == 1) {
      //   print(i_idx, ": visit_start = ", visit_start, ", visit_end = ", visit_end, ", patient_t_visits_week = ", patient_t_visits_week,
      //   ", lower_cutoff_visit_days[i_idx] = ", lower_cutoff_visit_days[i_idx], 
      //         ", n_patient_visits = ", n_patient_visits, ", t_idx = ", t_idx);
      // }
    }
    
    // Find last testing visit for nested cross-validation windows
    for (m in (n + 1):n_futures) {
      // Convert cutoff m to patient-specific study days
      array[n_curr_patients] int upper_cutoff_visit_days = calendar_date_to_study_date(patient_calendar_day[curr_patients], cutoff_calendar_day[m]);
    
      for (i_idx in 1:n_curr_patients) {
        int i = curr_patients[i_idx];  // Actual patient ID
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);
        int n_patient_visits = visit_end - visit_start + 1;
        
        // Extract and sort patient's visits by day (descending for last visit)
        array[n_patient_visits] int rev_patient_t_visits_week = reverse(t_patient_visits_week[visit_start:visit_end]);
        array[n_patient_visits] int rev_patient_t_visits_day = reverse(t_patient_visits_day[visit_start:visit_end]);
        
        // Find last visit on or before cutoff m
        int t_idx = 1;
        while (t_idx <= n_patient_visits && 
               rev_patient_t_visits_day[t_idx] > max(0, upper_cutoff_visit_days[i_idx])) {
          t_idx += 1;
        }
        
        if (t_idx <= n_patient_visits) {
          last_testing_visit_week[n, m, i] = rev_patient_t_visits_week[t_idx];
          testing_end_idx[n, m, i] = visit_end - t_idx + 1;
        } 

        // if (i_idx <= 2 && n == 1 && m == 2) {
        //   print(i_idx, ": visit_start = ", visit_start, ", visit_end = ", visit_end,
        //   ", upper_cutoff_visit_days[i_idx] = ", upper_cutoff_visit_days[i_idx],
        //         ", n_patient_visits = ", n_patient_visits, ", t_idx = ", t_idx);
        // }
      }
    }
  }
  
  return (first_testing_visit_week, last_testing_visit_week, testing_start_idx, testing_end_idx);
}