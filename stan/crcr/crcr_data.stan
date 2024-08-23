array[n_patients] int<lower = 0, upper = 1> confirmed_response;
array[n_patients] int<lower = 1> confirmed_response_week;
array[n_patients] int<lower = 0, upper = 1> confirmed_response_censored;

int n_tumor_covar;
int<lower = 0> n_covar; 

matrix[n_patients, n_covar] covar_design_matrix;
  