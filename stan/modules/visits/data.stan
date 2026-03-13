// ============================================================================
// VISITS MODULE DATA
// ============================================================================
// Clinical assessment scheduling for forecast periods.
// This module defines the observation interval for forecast burden measures
// (SLD, PSA, ctDNA, etc.) — separate from the weekly dynamics grid used for
// hazard computation.

int<lower=1> forecast_observation_interval;  // weeks between forecast assessments (typically 6)
