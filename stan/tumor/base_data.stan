int<lower = 0, upper = 1> fit_tumor_data;
int<lower = 0> sf_rep_T;
int<lower = 0, upper = 1> debug;
int<lower = 1, upper = n_patients> n_shards;

// Note: add_trial_level_baseline_hazard moved to modules/other_events/flags.stan as oe_enable_trial_baseline_hazard

array[sum(n_patient_visits)] int<lower = 1, upper = 5> recist;

array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive. The last week observed with no progression.
array[n_patients] int<lower = 0, upper = 1> right_censored;
array[n_patients] int<lower = 0> interval_censored; // The number of weeks after `pfs` that actual progression could have happened. E.g., zero means progression happened the next week.

array[n_patients] int<lower = 0> target_pfs; // PFS based on target tumor SLD only 
array[n_patients] int<lower = 0, upper = 1> target_right_censored;

array[n_patients] int<lower = 0> death_week;

int<lower = 0> n_covar; 
matrix[n_patients, n_covar] covar_design_matrix;
