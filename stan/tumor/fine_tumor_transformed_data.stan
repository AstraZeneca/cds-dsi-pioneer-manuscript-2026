
// Starting position of tumors for each patient in a flattened tumor array. This is used to traverse such data structures
// as n_measures.
// Diagram for patient_tumor_pos:
// [1, 3, 5, 6, ...]
//  ^  ^  ^  ^
//  |  |  |  |
//  |  |  |  Start of tumors for patient 4
//  |  |  Start of tumors for patient 3
//  |  Start of tumors for patient 2
//  Start of tumors for patient 1
array[n_patients + 1] int<lower = 1, upper = sum(n_patient_tumors) + 1> patient_tumor_pos = create_pos(n_patient_tumors);

// Starting position of measurements for each patient in a flattened measurement array
array[n_patients + 1] int<lower = 1, upper = sum(n_measures) + 1> patient_tumor_measure_pos = create_pos(n_measures, patient_tumor_pos);

// int min_all_t = min(t_measure); // Earliest measurement time across all patients
// int<lower = min_all_t> max_all_t = max(max(t_measure) + 1, extend_max_all_t); // Latest measurement time or extended time, whichever is greater
// int<lower = 0> max_t_width = max(t_measure) - min_all_t + 1;
array[sum(n_measures)] int<lower = 1> t_patient_measure_idx; // Index of each measurement time relative to the first measurement for each patient

print("max_all_t = ", max_all_t);

int<lower = 0> n_tumors = sum(n_patient_tumors); // Total number of tumors across all patients

// Information about missing measurements
array[n_tumors] int<lower = 0> n_missing_measures = calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors);  
array[sum(n_missing_measures)] int<lower = 1> t_patient_missing_measure_idx;  
array[n_tumors] int<lower = 0> n_full_measures; // Total number of measures (observed + missing) per tumor

// Information about screening measurements
array[n_tumors] int<lower = 0> n_screening_t = calc_n_screening_t(n_patient_tumors, n_measures, t_measure); // Number of pre-screening measures per tumor
array[n_patients] int<lower = 0> n_patient_screening_t = zeros_int_array(n_patients); // Number of pre-screening measures per patient

array[sum(n_measures)] int<lower = 1> t_measure_idx; // Measurement times indexed starting from 1

{
  int tumor_pos = 1;  // Current position in the tumor array
  int t_measure_pos = 1;  // Current position in the measurement array
  int t_missing_measure_pos = 1;  // Current position in the missing measurement array
  
  array[sum(n_missing_measures)] int t_missing_measure;  // Array to store missing measurement times
  array[sum(n_missing_measures)] int t_missing_measure_idx;  // Array to store indices of missing measurements
  
  // Calculate missing measurement times if there are any
  if (sum(n_missing_measures) > 0) {  
    t_missing_measure = calculate_t_missing_measure(n_measures, n_missing_measures, t_measure, n_patient_tumors);
  }
  
  // Iterate through all patients
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int save_t_measure_pos = t_measure_pos;
    int save_t_missing_measure_pos = t_missing_measure_pos;
    
    array[n_patient_tumors[i]] int min_t_idx;
    array[n_patient_tumors[i]] int max_t_idx;

    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1;  // End position of measurements for current tumor
      int t_missing_measure_end = t_missing_measure_pos + n_missing_measures[tumor_pos + j - 1] - 1;  // End position of missing measurements for current tumor
  
      // Calculate total number of measures (observed + missing) for current tumor
      n_full_measures[tumor_pos + j - 1] = n_measures[tumor_pos + j - 1] + n_missing_measures[tumor_pos + j - 1];  
  
      // Index measurement times relative to the earliest measurement time
      for (tp in t_measure_pos:t_measure_end) {
        t_measure_idx[tp] = t_measure[tp] - min_all_t + 1;
      }
  
      // Index missing measurement times relative to the earliest measurement time
      for (tp in t_missing_measure_pos:t_missing_measure_end) {
        t_missing_measure_idx[tp] = t_missing_measure[tp] - min_all_t + 1;
      }
  
      // Find minimum and maximum time indices for current tumor
      min_t_idx[j] = min(t_measure_idx[t_measure_pos:t_measure_end]);
      max_t_idx[j] = max(t_measure_idx[t_measure_pos:t_measure_end]);
  
      // Update positions for next tumor
      t_measure_pos = t_measure_end + 1;
      t_missing_measure_pos = t_missing_measure_end + 1;
    }
  
    // Calculate number of screening measurements for current patient
    n_patient_screening_t[i] = sum(n_screening_t[tumor_pos:tumor_end]);
  
    // Reset measurement positions
    t_measure_pos = save_t_measure_pos;
    t_missing_measure_pos = save_t_missing_measure_pos;
  
    // Iterate through all tumors for the current patient again
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1;  
      int t_missing_measure_end = t_missing_measure_pos + n_missing_measures[tumor_pos + j - 1] - 1;  
  
      // Index measurement times relative to the first measurement for each tumor
      for (tp in t_measure_pos:t_measure_end) {
        t_patient_measure_idx[tp] = t_measure_idx[tp] - min_t_idx[j] + 1;
      }
  
      // Index missing measurement times relative to the first measurement for each tumor
      for (tp in t_missing_measure_pos:t_missing_measure_end) {
        t_patient_missing_measure_idx[tp] = t_missing_measure_idx[tp] - min_t_idx[j] + 1;
      }
  
      // Update positions for next tumor
      t_measure_pos = t_measure_end + 1;
      t_missing_measure_pos = t_missing_measure_end + 1;
    }
  
    // Update tumor position for next patient
    tumor_pos = tumor_end + 1;
  }
}