int<lower = 0> n_patients;
int<lower = 1> n_measures;

array[n_patients] int<lower = 0> n_patient_tumors; // The overall number of tumors observed at any time per patient 
//array[sum(n_patient_tumors)] int<lower = 0, upper = n_measures> n_tumor_measures; // The number of times each tumor is measured