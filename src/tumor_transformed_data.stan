real delta = 1e-9;

// int<lower = 0> n_all_tumor_measures; // The sum of all measures for all tumors if we observed all of them for each interval

int min_all_t = min(t_measure);
int<lower = min_all_t> max_all_t = max(t_measure);
array[sum(n_measures)] int<lower = 1> patient_t_measure_idx; // The index of each t relative to the first t per patient
array[n_patients] int<lower = 0> patient_max_t_width; // The number of intervals from the first to the last measurement per patient

int<lower = 0> n_tumors = sum(n_patient_tumors);

// Same as n_measures and t_measure but for the missing measurement intervals 
array[n_tumors] int<lower = 0> n_missing_measures = calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors); 
array[sum(n_missing_measures)] int<lower = 1> patient_t_missing_measure_idx; 
array[n_tumors] int<lower = 0> n_full_measures;

array[n_tumors] int<lower = 0> n_screening_t = calc_n_screening_t(n_patient_tumors, n_measures, t_measure); 
array[n_patients] int<lower = 0> n_patient_screening_t =  rep_array(0, n_patients);

array[sum(n_measures)] int<lower = 1> t_measure_idx;

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
      
      // print(patient_t_measure_idx[t_measure_pos:t_measure_end]);
      
      t_measure_pos = t_measure_end + 1;
      t_missing_measure_pos = t_missing_measure_end + 1;
    }
    
    tumor_pos = tumor_end + 1;
  }
}

array[max(patient_max_t_width)] real all_measure_t; // This is for the GP "proximity" between size measurement time intervals.

for (t in 1:max(patient_max_t_width)) {
  all_measure_t[t] = t / 12.0; 
} 

// // n_all_tumor_measures = to_int(to_row_vector(patient_max_t_width) * to_vector(n_patient_tumors));
// n_all_tumor_measures = sum(n_measures) + sum(n_missing_measures); 

// if (sum(n_missing_measures) > 0) { 
//   array[sum(n_missing_measures)] int t_missing_measure = calculate_t_missing_measure(n_measures, n_missing_measures, t_measure, n_patient_tumors);
//   patient_t_missing_measure_idx =
//   // patient_t_missing_measure_idx = calculate_t_missing_measure(n_measures, n_missing_measures, patient_t_measure_idx, n_patient_tumors);
// }