// frac/flags.stan
// Fraction mix module gating flags.

int<lower=0,upper=1> enable_pop_cov_frac;
int<lower=0,upper=1> enable_trial_intercept_frac;
int<lower=0,upper=1> enable_trial_cov_frac;
int<lower=0,upper=1> enable_patient_intercept_frac;
int<lower=0,upper=1> enable_patient_cov_frac;
