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
array[n_trials + 1] int<lower = 1, upper = n_patients + 1> trial_patient_pos;
trial_patient_pos[1] = 1;

for (s in 1:n_trials) {
  trial_patient_pos[s + 1] = sum(n_trial_patients[:s]) + 1;  
}  

// Starting position of tumors for each patient in a flattened tumor array
// Diagram for patient_tumor_pos:
// [1, 3, 5, 6, ...]
//  ^  ^  ^  ^
//  |  |  |  |
//  |  |  |  Start of tumors for patient 4
//  |  |  Start of tumors for patient 3
//  |  Start of tumors for patient 2
//  Start of tumors for patient 1
array[n_patients + 1] int<lower = 1, upper = sum(n_patient_tumors) + 1> patient_tumor_pos;
patient_tumor_pos[1] = 1;

for (i in 1:n_patients) {
  patient_tumor_pos[i + 1] = sum(n_patient_tumors[:i]) + 1;  
}  

// Starting position of measurements for each patient in a flattened measurement array
array[n_patients + 1] int<lower = 1, upper = sum(n_measures) + 1> patient_tumor_measure_pos;
patient_tumor_measure_pos[1] = 1;

for (i in 2:(n_patients + 1)) {
  patient_tumor_measure_pos[i] = sum(n_measures[:(patient_tumor_pos[i] - 1)]) + 1;  
}  

real delta = 1e-9; // Small value used for GP modeling to avoid numerical issues

int min_all_t = min(t_measure); // Earliest measurement time across all patients
int<lower = min_all_t> max_all_t = max(max(t_measure) + 1, extend_max_all_t); // Latest measurement time or extended time, whichever is greater
array[sum(n_measures)] int<lower = 1> patient_t_measure_idx; // Index of each measurement time relative to the first measurement for each patient
array[n_patients] int<lower = 0> patient_max_t_width; // Number of time intervals between first and last measurement for each patient

print("max_all_t = ", max_all_t);

int<lower = 0> n_tumors = sum(n_patient_tumors); // Total number of tumors across all patients

// Information about missing measurements
array[n_tumors] int<lower = 0> n_missing_measures = calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors);  
array[sum(n_missing_measures)] int<lower = 1> patient_t_missing_measure_idx;  
array[n_tumors] int<lower = 0> n_full_measures; // Total number of measures (observed + missing) per tumor

// Information about screening measurements
array[n_tumors] int<lower = 0> n_screening_t = calc_n_screening_t(n_patient_tumors, n_measures, t_measure); // Number of pre-screening measures per tumor
array[n_patients] int<lower = 0> n_patient_screening_t = rep_array(0, n_patients); // Number of pre-screening measures per patient

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
        patient_t_measure_idx[tp] = t_measure_idx[tp] - min_t_idx[j] + 1;
      }
      
      for (tp in t_missing_measure_pos:t_missing_measure_end) {
        patient_t_missing_measure_idx[tp] = t_missing_measure_idx[tp] - min_t_idx[j] + 1;
      }
      
      t_measure_pos = t_measure_end + 1;
      t_missing_measure_pos = t_missing_measure_end + 1;
    }
    
    tumor_pos = tumor_end + 1;
  }
}

// Array of measurement times used for GP modeling
array[max(patient_max_t_width)] real all_measure_t;

for (t in 1:max(patient_max_t_width)) {
  all_measure_t[t] = t / 12.0; // Scaling factor for time intervals. The 12 here is arbitrary (if it actually had any meaning at one point).
}

// Variables for handling separate baseline and proportional hazards
int n_base_separate_trials = separate_baseline_hazard ? n_trials : 1;
int n_prop_separate_trials = (1 - no_prop_hazard) * (separate_prop_hazard ? n_trials : 1);

print("n_base_separate_trials = ", n_base_separate_trials, ", n_prop_separate_trials = ", n_prop_separate_trials);