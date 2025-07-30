int<lower = 0, upper = 1> fit_tumor_data;
int<lower = 0> sf_rep_T;
int<lower = 0, upper = 1> debug;
int<lower = 0, upper = 1> pop_growth_lag_param_only;
int<lower = 0, upper = 1> pop_initial_states_param_only;
int<lower = 0, upper = 1> pop_rates_param_only;
int<lower = 0, upper = 1> pop_rho_param_only; 
int<lower = 0, upper = 1> pop_covar_coef_only;
int<lower = 0, upper = 1> independ_long_process_noise;
int<lower = 0, upper = 1> independ_cross_process_noise;
int<lower = 1, upper = n_patients> train_patients_pos, train_patients_end;
int<lower = 1, upper = n_patients> n_shards;

int<lower = 0, upper = 1> add_trial_level_net_rate; 
int<lower = 0, upper = 1> add_trial_level_baseline_hazard;
int<lower = 0, upper = 1> add_trial_level_prop;

array[sum(n_patient_visits)] int<lower = 1, upper = 5> recist;

array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive. The last week observed with no progression.
array[n_patients] int<lower = 0, upper = 1> right_censored;
array[n_patients] int<lower = 0> interval_censored; // The number of weeks after `pfs` that actual progression could have happened. E.g., zero means progression happened the next week.

array[n_patients] int<lower = 0> target_pfs; // PFS based on target tumor SLD only 
array[n_patients] int<lower = 0, upper = 1> target_right_censored;

array[n_patients] int<lower = 0> death_week;

int<lower = 0> n_covar; 
matrix[n_patients, n_covar] covar_design_matrix;
