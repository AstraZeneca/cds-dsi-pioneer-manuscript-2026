int<lower = 1> n_trials;
int<lower = 0> n_patients;

int<lower = 1> n_tumor_locations;

array[n_patients] int<lower = 1, upper = n_trials> patient_trial;

// How many tumors were measured. Tumors that weren't detected or disappeared = 0 size. 
// [num_tumors_1, ..., num_tumors_i, ...]
array[n_patients] int<lower = 1> n_patient_tumors; 

// How many times were tumors measured per tumor.
// [..., (measures_{i,1} ..., measures_{i,n_patient_tumors[i]}), ...]
array[sum(n_patient_tumors)] int<lower = 1> n_measures; 

array[sum(n_patient_tumors)] int<lower = 1, upper = n_tumor_locations> tumor_location; 

// The periods of each measurement per tumor. 
// [..., (t_{i,j,1}, ..., t_{i,j,n_measures_{i,j}}}), ... ]
array[sum(n_measures)] int t_measure; // t <= 0 are screening (baseline) assessments (measurements)

// [..., ((tumor_size_{i,1,1}, ..., tumor_size_{i, 1, n_measures_i}), ..., (..., tumor_size_{i,j,t},...), ...), ...  ] 
vector<lower = 0>[sum(n_measures)] tumor_size; // cm 