// tr/hyperparams.stan
// Data (hyperparameter) declarations for Total Rate (tr) module.
// Multi-level hierarchy: hyperparameters are arrays indexed by level (1..n_levels)

real tr_loc_pop_mean;         // mean prior for population log total rate (log scale)
real<lower=0> tr_loc_pop_sd;  // sd prior for population log total rate

// Hierarchical intercept prior scale hyperparameters - one per level
array[n_levels] real<lower=0> tr_sd_level_intercept_sd;

// Student-t hierarchy: nu prior hyperparameters (one per level)
array[n_levels] real<lower=0> tr_nu_level_prior_alpha;
array[n_levels] real<lower=0> tr_nu_level_prior_beta;

// QR-space coefficient hyperparameters (applied in model block)
vector[n_covar] tr_coef_qr_pop_mean;        // mean for QR coefficients (typically 0)
vector<lower=0>[n_covar] tr_coef_qr_pop_sd; // sd for QR coefficients (typically 1)

// Hierarchical slope SD hyperpriors - one vector per level
array[n_levels] row_vector<lower=0>[n_covar] tr_sd_level_slope_sd;

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
