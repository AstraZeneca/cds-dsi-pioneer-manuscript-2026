// Model settings
int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
int<lower = 0, upper = 1> crcr_ignore_interval_censoring; // Treat observed confirmed response week as true and ignore t_measure.
int<lower = 0, upper = 1> pfs_ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
int<lower = 0, upper = 1> gen_log_lik; // Calculate log likelihood for LOO-CV
int<lower = 0, upper = 1> prior_sense; // For prior sensitivity using {priorsense}
int<lower = 0, upper = 1> pfs_only; // Ignore the confirmed response model; don't include as predictor in proportional hazard.
int<lower = 0, upper = 1> no_prop_hazard;
int<lower = 0, upper = 1> no_tumor_effects; // Don't include tumor size as a predictor in proportional hazard.

// Hierarchical settings 
int<lower = 0, upper = 1> add_trial_level_baseline_hazard;
int<lower = 0, upper = 1> add_trial_level_prop_hazard;
int<lower = 0, upper = 1> separate_baseline_hazard;
int<lower = 0, upper = 1> separate_prop_hazard;

// This is the data that is shared with the tumor model 
#include "../base_data.stan"

#include "../bootstrap/leave_out_trial_bootstrap_data.stan"

// Calculating log likelihood for a single trial. Useful if you want to compare the preformance of a model using a single trial with one that is multilevel. 
int<lower = 0, upper = n_trials> log_lik_trial;

array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive. The last week observed with no progression.
array[n_patients] int<lower = 0, upper = 1> right_censored;
array[n_patients] int<lower = 0> interval_censored; // The number of weeks after `pfs` that actual progression could have happened. E.g., zero means progression happened the next week.
array[n_patients] int<lower = 1> admin_right_censored_week;

#include "../crcr/crcr_data.stan" 

// Hyperparam
#include "../baseline_hazard/baseline_hazard_hyperparam.stan"
#include "../crcr/crcr_hyperparam.stan"

array[separate_prop_hazard ? n_trials : 1] vector<lower = 0>[2] tumor_stim_pop_coef_sd;
array[separate_prop_hazard ? n_trials : 1] real conf_resp_effect_mean;
array[separate_prop_hazard ? n_trials : 1] real<lower = 0> conf_resp_effect_sd;
array[separate_prop_hazard ? n_trials : 1] vector[n_covar] covar_effect_mean;
array[separate_prop_hazard ? n_trials : 1] vector<lower = 0>[n_covar] covar_effect_sd;

// Multilevel hyperparameters for the proportional hazard parameters 
real<lower = 0> covar_trial_sd_sd;
real<lower = 0> covar_trial_corr_eta; // Correlation between parameters