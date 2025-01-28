array[n_trials] int<lower = 0, upper = n_patients> n_trial_patients = rep_array(0, n_trials); // How many patients per trial

for (i in 1:n_patients) {
  n_trial_patients[patient_trial[i]] += 1;
}

print("n_trial_patients = ", n_trial_patients);

array[n_trials + 1] int<lower = 1, upper = n_patients + 1> trial_patient_pos;
trial_patient_pos[1] = 1;

for (s in 1:n_trials) {
  trial_patient_pos[s + 1] = sum(n_trial_patients[:s]) + 1; 
} 

// Positions of the first n_measure per patient
array[n_patients + 1] int<lower = 1, upper = sum(n_patient_tumors) + 1> patient_tumor_pos;
patient_tumor_pos[1] = 1;

for (i in 1:n_patients) {
  patient_tumor_pos[i + 1] = sum(n_patient_tumors[:i]) + 1; 
} 

// Positions of the first t_measure or tumor size per patient
array[n_patients + 1] int<lower = 1, upper = sum(n_measures) + 1> patient_tumor_measure_pos;
patient_tumor_measure_pos[1] = 1;

for (i in 1:n_patients) {
  patient_tumor_measure_pos[i + 1] = sum(n_measures[:patient_tumor_pos[i]]) + 1; 
} 

real delta = 1e-9; // Used for GP modeling

int min_all_t = min(t_measure);
int<lower = min_all_t> max_all_t = max(max(t_measure) + 1, extend_max_all_t);
array[sum(n_measures)] int<lower = 1> patient_t_measure_idx; // The index of each t relative to the first t per patient
array[n_patients] int<lower = 0> patient_max_t_width; // The number of intervals from the first to the last measurement per patient

print("max_all_t = ", max_all_t);

int<lower = 0> n_tumors = sum(n_patient_tumors);

// Same as n_measures and t_measure but for the missing measurement intervals 
array[n_tumors] int<lower = 0> n_missing_measures = calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors); 
array[sum(n_missing_measures)] int<lower = 1> patient_t_missing_measure_idx; 
array[n_tumors] int<lower = 0> n_full_measures; // Missing and observed measures

array[n_tumors] int<lower = 0> n_screening_t = calc_n_screening_t(n_patient_tumors, n_measures, t_measure); // number of pre-screening measures per tumor
array[n_patients] int<lower = 0> n_patient_screening_t = rep_array(0, n_patients); // ...per patient

array[sum(n_measures)] int<lower = 1> t_measure_idx; // measure ts starting at 1.

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

array[max(patient_max_t_width)] real all_measure_t; // This is for the GP "proximity" between size measurement time intervals.

for (t in 1:max(patient_max_t_width)) {
  all_measure_t[t] = t / 12.0; // Why 12? Our intervals are weeks not months. Doesn't matter.  
}

int n_base_separate_trials = separate_baseline_hazard ? n_trials : 1;
int n_prop_separate_trials = (1 - no_prop_hazard) * (separate_prop_hazard ? n_trials : 1);

print("n_base_separate_trials = ", n_base_separate_trials, ", n_prop_separate_trials = ", n_prop_separate_trials);