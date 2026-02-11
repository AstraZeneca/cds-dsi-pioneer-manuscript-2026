// ============================================================================
// STATE COMPUTATION FLAGS
// ============================================================================

// Whether to compute states on full time grid [n_patients x max_t_width]
// Set to 1 if: process noise is enabled, OR downstream modules need time-varying states
// (e.g., other_events module with time-varying biomarker covariates)
int<lower=0, upper=1> enable_states_full_grid;

// ============================================================================
// ENDPOINT CONFIGURATION
// ============================================================================

// Quantiles to estimate (e.g., 0.25, 0.5, 0.75 for Q1, median, Q3)
int<lower = 0> n_pfs_quantiles;
vector<lower = 0, upper = 1>[n_pfs_quantiles] pfs_quantiles;

int<lower = 0> n_pfs_timepoints;
array[n_pfs_timepoints] int<lower = 0> pfs_timepoints; // In months

// Conditioning groups or strata to estimate outcomes for a particular covar
int<lower = 0> n_cond_group;
array[n_cond_group] int<lower = 0> cond_group_size;
array[sum(cond_group_size)] int <lower = 1, upper = n_patients> cond_group;
