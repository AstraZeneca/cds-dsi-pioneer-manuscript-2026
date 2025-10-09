// tr/flags.stan
// Module: Total Rate (tr)
// Defines gating flags for total rate module. Population intercept always on.
// See docs Section 2.3 / 2.4.

int<lower=0,upper=1> enable_pop_cov_tr;           // population covariate linear model
int<lower=0,upper=1> enable_trial_intercept_tr;    // trial random intercept hierarchy
int<lower=0,upper=1> enable_trial_cov_tr;          // trial slope deviations (QR)
int<lower=0,upper=1> enable_patient_intercept_tr;  // patient random intercept hierarchy
int<lower=0,upper=1> enable_patient_cov_tr;        // patient slope deviations (QR)
