int n_tumor_covar_col = calc_n_tumor_covar_col(tumor_hazard_type); // number of columns in covariates design matrix, excluding the intercept. 

int max_measures = 2; 
array[n_tumors, max_measures] int<lower = min(t_measure), upper = max(t_measure)> tumor_covar_t; // ts (weeks) of the assessments used in the covar design matrix
array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> patient_max_2nd_tumor_t; // What t is the second assessment in the covar design matrix 
matrix[n_tumors, n_tumor_covar_col] tumor_covar; // This is the design matrix with the covar in the first two (or _n_) assessments.
matrix[n_tumors, n_tumor_covar_col] uncentered_tumor_covar; // Just scaled
vector[n_tumor_covar_col] tumor_covar_mean;
vector<lower = 0>[n_tumor_covar_col] tumor_covar_sd = rep_vector(0, n_tumor_covar_col);

(tumor_covar_t, patient_max_2nd_tumor_t, tumor_covar, uncentered_tumor_covar, tumor_covar_mean, tumor_covar_sd) = 
  prepare_early_tumors_covar(tumor_size, tumor_hazard_type, n_tumors, n_patient_tumors, n_measures, t_measure, n_screening_t, max_measures); 
