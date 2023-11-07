int<lower = 0> n_patients;

// How many times were tumors measured per patient.
// [measures_1, ..., measures_i, ..., measures_N]
array[n_patients] int<lower = 2> n_measures; // Must have at least one baseline measure and one follow-up measure. 

// The periods of each measurement per patient. t=0 is the baseline which is always present but not included in t_measures.
// So for each patient there will be (n_measures[i] - 1) t's in this array.
// [..., (t_{i,1}, ..., t_{i,measures_i - 1}), ... ]
array[sum(n_measures) - n_patients] int<lower = 1> t_measure;  

// How many tumors were measured. Tumors that weren't detected or disappeared = 0 size. 
// [num_tumors_1, ..., num_tumors_i, ...]
array[n_patients] int<lower = 1> n_patient_tumors; 