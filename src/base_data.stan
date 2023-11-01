int<lower = 0> n_patients;

// How many times were tumors measured per patient.
// [measures_1, ..., measures_i, ..., measures_N]
array[n_patients] int<lower = 1> n_measures; 

// The periods of each measurement per patient. t=0 is the baseline.
// [..., (t_{i,1}, ..., t_{i,measures_i}), ... ]
array[sum(n_measures)] int<lower = 0> t_measure; 

// How many tumors were measured. Tumors that weren't detected or disappeared = 0 size. 
// [num_tumors_1, ..., num_tumors_i, ...]
array[n_patients] int<lower = 0> n_patient_tumors; 