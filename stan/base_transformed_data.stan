// Number of patients in each trial
array[n_trials] int<lower = 0, upper = n_patients> n_trial_patients = rep_array(0, n_trials);

for (i in 1:n_patients) {
  n_trial_patients[patient_trial[i]] += 1;
}

print("n_trial_patients = ", n_trial_patients);

// Starting position of patients for each trial in a flattened patient array
// Diagram for trial_patient_pos:
// [1, 4, 7, 12, ...]
//  ^  ^  ^  ^
//  |  |  |  |
//  |  |  |  Start of patients in trial 4
//  |  |  Start of patients in trial 3
//  |  Start of patients in trial 2
//  Start of patients in trial 1
array[n_trials + 1] int<lower = 1, upper = n_patients + 1> trial_patient_pos = create_pos(n_trial_patients);

// Starting position of tumors for each patient in a flattened tumor array
// Diagram for patient_tumor_pos:
// [1, 3, 5, 6, ...]
//  ^  ^  ^  ^
//  |  |  |  |
//  |  |  |  Start of tumors for patient 4
//  |  |  Start of tumors for patient 3
//  |  Start of tumors for patient 2
//  Start of tumors for patient 1
array[n_patients + 1] int<lower = 1, upper = sum(n_patient_tumors) + 1> patient_tumor_pos = create_pos(n_patient_tumors);

// Starting position of measurements for each trial in a flattened measurement array
array[n_trials + 1] int<lower = 1, upper = sum(n_patient_visits) + 1> trial_visit_pos = create_pos(n_patient_visits, trial_patient_pos);

// Starting position of measurements for each patient in a flattened measurement array
array[n_patients + 1] int<lower = 1, upper = sum(n_measures) + 1> patient_tumor_measure_pos = create_pos(n_measures, patient_tumor_pos);
array[n_patients + 1] int<lower = 1, upper = sum(n_patient_visits) + 1> patient_visit_pos = create_pos(n_patient_visits);

real delta = 1e-5; // Small value used for GP modeling to avoid numerical issues

int min_all_t = min(t_measure); // Earliest measurement time across all patients
int<lower = min_all_t> max_all_t = max(max(t_measure) + 1, extend_max_all_t); // Latest measurement time or extended time, whichever is greater
int<lower = 0> max_t_width = max(t_measure) - min_all_t + 1;
array[n_patients] int<lower = 0> patient_max_t_width; // Number of time intervals between first and last measurement for each patient
array[sum(n_measures)] int<lower = 1> t_patient_measure_idx; // Index of each measurement time relative to the first measurement for each patient
array[sum(n_patient_visits)] int<lower = 1> t_patient_visit_idx; // Index of each patient visit relative to the first visit for each patient
array[sum(n_patient_visits)] int<lower = 1, upper = max_t_width> t_visit_trial_idx = id2idx(t_patient_visits, trial_visit_pos); 

print("max_all_t = ", max_all_t);

int<lower = 0> n_tumors = sum(n_patient_tumors); // Total number of tumors across all patients

// Information about missing measurements
array[n_tumors] int<lower = 0> n_missing_measures = calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors);  
array[sum(n_missing_measures)] int<lower = 1> t_patient_missing_measure_idx;  
array[n_tumors] int<lower = 0> n_full_measures; // Total number of measures (observed + missing) per tumor

// Information about screening measurements
array[n_tumors] int<lower = 0> n_screening_t = calc_n_screening_t(n_patient_tumors, n_measures, t_measure); // Number of pre-screening measures per tumor
array[n_patients] int<lower = 0> n_patient_screening_t = zeros_int_array(n_patients); // Number of pre-screening measures per patient
array[n_patients] int<lower = 0> n_patient_screening_visits = zeros_int_array(n_patients);

array[sum(n_measures)] int<lower = 1> t_measure_idx; // Measurement times indexed starting from 1

{
  int tumor_pos = 1;
  int t_measure_pos = 1;
  int t_missing_measure_pos = 1;
  
  array[sum(n_missing_measures)] int t_missing_measure;
  array[sum(n_missing_measures)] int t_missing_measure_idx;
  
  if (sum(n_missing_measures) > 0) { 
    t_missing_measure = calculate_t_missing_measure(n_measures, n_missing_measures, t_measure, n_patient_tumors);
  }
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int save_t_measure_pos = t_measure_pos;
    int save_t_missing_measure_pos = t_missing_measure_pos;
    
    array[n_patient_tumors[i]] int min_t_idx;
    array[n_patient_tumors[i]] int max_t_idx;
   
    int first_visit = t_patient_visits[patient_visit_pos[i]];
    
    int curr_patient_visit_pos, curr_patient_visit_end;
    (curr_patient_visit_pos, curr_patient_visit_end) = get_pos(patient_visit_pos, i);
    
    for (v in curr_patient_visit_pos:curr_patient_visit_end) {
      t_patient_visit_idx[v] = t_patient_visits[v] - first_visit + 1;
      
      if (t_patient_visits[v] <= 0) {
        n_patient_screening_visits[i] += 1;
      }
    }
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1;
      int t_missing_measure_end = t_missing_measure_pos + n_missing_measures[tumor_pos + j - 1] - 1;
      
      n_full_measures[tumor_pos + j - 1] = n_measures[tumor_pos + j - 1] + n_missing_measures[tumor_pos + j - 1]; 
      
      for (tp in t_measure_pos:t_measure_end) {
        t_measure_idx[tp] = t_measure[tp] - min_all_t + 1;
      }
      
      for (tp in t_missing_measure_pos:t_missing_measure_end) {
        t_missing_measure_idx[tp] = t_missing_measure[tp] - min_all_t + 1;
      }
      
      min_t_idx[j] = min(t_measure_idx[t_measure_pos:t_measure_end]);
      max_t_idx[j] = max(t_measure_idx[t_measure_pos:t_measure_end]);
      
      t_measure_pos = t_measure_end + 1;
      t_missing_measure_pos = t_missing_measure_end + 1;
    }
    
    n_patient_screening_t[i] = sum(n_screening_t[tumor_pos:tumor_end]);
    
    patient_max_t_width[i] = max(max_t_idx) - min(min_t_idx) + 1;
    
    t_measure_pos = save_t_measure_pos;
    t_missing_measure_pos = save_t_missing_measure_pos;
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1; 
      int t_missing_measure_end = t_missing_measure_pos + n_missing_measures[tumor_pos + j - 1] - 1; 
      
      for (tp in t_measure_pos:t_measure_end) {
        t_patient_measure_idx[tp] = t_measure_idx[tp] - min_t_idx[j] + 1;
      }
      
      for (tp in t_missing_measure_pos:t_missing_measure_end) {
        t_patient_missing_measure_idx[tp] = t_missing_measure_idx[tp] - min_t_idx[j] + 1;
      }
      
      t_measure_pos = t_measure_end + 1;
      t_missing_measure_pos = t_missing_measure_end + 1;
    }
    
    tumor_pos = tumor_end + 1;
  }
}

// Array of measurement times used for GP modeling
array[max_t_width] real all_measure_t;

for (t in 1:max_t_width) {
  all_measure_t[t] = t / 12.0; // Scaling factor for time intervals. The 12 here is arbitrary (if it actually had any meaning at one point).
}

