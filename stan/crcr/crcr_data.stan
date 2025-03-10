array[n_patients] int<lower = 0, upper = 1> confirmed_response; // Observed response if classified
array[n_patients] int<lower = 1> confirmed_response_day; 
array[n_patients] int<lower = 1> confirmed_response_week; 
array[n_patients] int<lower = 0, upper = 1> confirmed_response_censored; // Right censored
array[n_patients] int<lower = 0> confirmed_response_interval_censored;

int<lower = 1> extend_max_confresp_week; // Extrapolate survival inference between the last observed week in the data.

int n_tumor_covar; // Number of tumor size sum covariates
int<lower = 0> n_covar; // Number of covariates other than tumor size

matrix[n_patients, n_covar] covar_design_matrix; // Design matrix for covariates other than tumor size

array[n_patients] int<lower = 0, upper = 1> orr_pop;