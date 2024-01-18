int<lower = 0, upper = 1> fit_tumor_data; // Sample of prior prediction only
int<lower = 0, upper = 1> predict_missing_sizes; 
int<lower = 0, upper = 1> gen_tumor_sizes; // Generate data given prior/posterior of model parameters
int<lower = 0, upper = 1> multilevel_patient;
int<lower = 0, upper = 1> multilevel_tumor;

// Hyperparam
real<lower = 0> pop_tumor_gp_rho_alpha;
real<lower = 0> pop_tumor_gp_rho_beta;