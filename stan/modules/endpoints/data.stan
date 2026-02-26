// ============================================================================
// Endpoint Computation Data
// Shared between full model (sf-ssm-log-space.stan) and standalone
// multistate model (ms-standalone.stan).
// ============================================================================

// Quantiles to estimate (e.g., 0.25, 0.5, 0.75 for Q1, median, Q3)
int<lower=0> n_pfs_quantiles;
vector<lower=0, upper=1>[n_pfs_quantiles] pfs_quantiles;

int<lower=0> n_pfs_timepoints;
array[n_pfs_timepoints] int<lower=0> pfs_timepoints;  // In months

// Conditioning groups or strata to estimate outcomes for a particular covariate
int<lower=0> n_cond_group;
array[n_cond_group] int<lower=0> cond_group_size;
array[sum(cond_group_size)] int<lower=1, upper=n_patients> cond_group;
