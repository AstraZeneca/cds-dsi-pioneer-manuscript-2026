// tr/hyperparams.stan
// Data (hyperparameter) declarations for Total Rate (tr) module.
// NEW naming aligned with conventions:
//   pop_log_total_rate_mean            -> tr_loc_pop_mean
//   pop_log_total_rate_sd              -> tr_loc_pop_sd
//   trial_log_total_rate_sd_sd         -> tr_sd_trial_intercept_sd
//   patient_log_total_rate_sd_sd       -> tr_sd_patient_intercept_sd
// Future (covariate) hyperparams (not yet supplied by R side):
//   tr_coef_pop_mean, tr_coef_pop_sd
//   tr_sd_trial_slope_sd, tr_sd_patient_slope_sd
// These are declared now to stabilize interface; they will be populated when Step 9 updates R code.

real tr_loc_pop_mean;                 // mean prior for population log total rate (log scale)
real<lower=0> tr_loc_pop_sd;          // sd prior for population log total rate

// Hierarchical intercept prior scale hyperparameters
real<lower=0> tr_sd_trial_intercept_sd;    // prior SD for trial intercept SD
real<lower=0> tr_sd_patient_intercept_sd;  // prior SD for patient intercept SD

// QR-space coefficient hyperparameters (applied in model block)
vector[n_covar] tr_coef_qr_pop_mean;                   // mean for QR coefficients (typically 0)
vector<lower=0>[n_covar] tr_coef_qr_pop_sd;            // sd for QR coefficients (typically 1)

// Hierarchical slope SD hyperpriors (per covariate)
row_vector<lower=0>[n_covar] tr_sd_trial_slope_sd;    // prior SD for each trial-level slope SD
row_vector<lower=0>[n_covar] tr_sd_patient_slope_sd;  // prior SD for each patient-level slope SD

// Process noise hyperparameters (for AR(1) process on measurements)
// Patient-level process noise prior hyperparameters
real tr_log_sd_pop_process_noise_mean;         // Prior mean for log(σ) at population level (for patient process)
real tr_log_sd_pop_process_noise_sd;           // Prior SD for log(σ) at population level (for patient process)
real tr_logit_phi_pop_process_noise_mean;      // Prior mean for logit(φ) at population level (for patient process)
real tr_logit_phi_pop_process_noise_sd;        // Prior SD for logit(φ) at population level (for patient process)
real<lower=0> tr_log_sd_patient_process_noise_sd;  // Patient variation in log(σ)
real<lower=0> tr_phi_patient_process_noise_sd;     // Patient variation in logit(φ)

// Population-level time-varying process noise hyperparameters (shared across all patients)
real tr_log_sd_pop_process_noise_pop_mean;     // Prior mean for pop-level log(σ)
real tr_log_sd_pop_process_noise_pop_sd;       // Prior SD for pop-level log(σ)
real tr_logit_phi_pop_process_noise_pop_mean;  // Prior mean for pop-level logit(φ)
real tr_logit_phi_pop_process_noise_pop_sd;    // Prior SD for pop-level logit(φ)
