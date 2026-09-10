// --- LFO CV specific ---
int<lower = 1> n_cutoffs;
array[n_cutoffs] int<lower = 1> cutoff_calendar_day;
int<lower = 1> max_n_rows;            // Number of n rows to compute (1 for exact LFO, n_cutoffs for PSIS)
int<lower = 1> max_forecast_horizon;  // Number of m columns to compute (2 for exact LFO, n_cutoffs for PSIS)
int<lower = 1> lfo_eval_trial;        // Trial index for held-out log-lik evaluation
