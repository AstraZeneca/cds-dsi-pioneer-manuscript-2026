// Include fragment: transformed data block content for SSLS models
// This module is biomarker-agnostic - observation-specific normalization
// (e.g., normalized_sld) should be in the respective observation modules.
// Visit infrastructure (patient_visit_pos, n_total_visits, forecast_visits_pos, etc.)
// is in _visit_transformed_data.stan.

// ============================================================================
// Visit indexing
// ============================================================================

array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
int<lower = 1> n_total_forecast_visits = get_pos_total_size(forecast_visits_pos);

// RECIST response categories
// ============================================================================

int CR = 1;
int PR = 2;
int SD = 3;
int PD = 4;

int NT_CR = 1;
int NT_PD = 2;

#include "_qr_decomposition.stan"

// ============================================================================
// Create cumulative sum indicator matrix for batched state computation
// ============================================================================

int max_unique_visit = max(pop_unique_visits) - min(pop_unique_visits) + 1;
matrix[n_pop_unique_visits, max_unique_visit] visit_cumsum_mat;
{
  int min_visit = min(pop_unique_visits);
  visit_cumsum_mat = rep_matrix(0, n_pop_unique_visits, max_unique_visit);

  for (v in 1:n_pop_unique_visits) {
    int shifted_visit = pop_unique_visits[v] - min_visit + 1;
    visit_cumsum_mat[v, :shifted_visit] = ones_row_vector(shifted_visit);
  }
}

// ============================================================================
// Endpoints: PFS conditioning groups
// ============================================================================

#include "modules/endpoints/transformed_data.stan"
