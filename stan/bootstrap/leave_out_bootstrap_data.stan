int<lower = 0, upper = n_trials> leave_out_trial;
int<lower = 0> n_bootstrap_sample;
int<lower = 0> n_bootstrap_cr_maturity_rates;
vector<lower = 0, upper = 1>[n_bootstrap_cr_maturity_rates] bootstrap_cr_maturity_rates;
int<lower = 0> n_bootstrap_pfs_maturity_rates;
vector<lower = 0, upper = 1>[n_bootstrap_pfs_maturity_rates] bootstrap_pfs_maturity_rates;
  