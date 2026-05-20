// ============================================================================
// STATE COMPUTATION FLAGS
// ============================================================================

// Whether to compute states on full time grid [n_patients x max_t_width]
// Set to 1 if: process noise is enabled, OR downstream modules need time-varying states
// (e.g., other_events module with time-varying biomarker covariates)
int<lower=0, upper=1> enable_states_full_grid;

#include "modules/endpoints/data.stan"
