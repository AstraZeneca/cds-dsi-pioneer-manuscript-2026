// Include fragment: transformed data block content for SSLS models

// ============================================================================
// Basic tumor measurement setup
// ============================================================================

// Normalize SLD by baseline for each patient
vector<lower = 0>[sum(n_patient_visits)] normalized_sld;
for (i in 1:n_patients) {
  int visit_pos, visit_end;
  (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
  normalized_sld[visit_pos:visit_end] = sum_tumor_size[visit_pos:visit_end] / sum_tumor_size[visit_pos]; 
}

// ============================================================================
// Visit indexing and counts
// ============================================================================

int<lower = 1> n_total_visits = sum(n_patient_visits);
array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
int<lower = 1> n_total_forecast_visits = get_pos_total_size(forecast_visits_pos);

// ============================================================================
// RECIST response categories
// ============================================================================

int CR = 1;
int PR = 2;
int SD = 3;
int PD = 4;

int NT_CR = 1;
int NT_PD = 2;

// ============================================================================
// Covariate design matrix (QR decomposition)
// ============================================================================

matrix[n_patients, n_covar] Q_covar_design_matrix;
matrix[n_covar, n_covar] R_covar_design_matrix;

if (n_covar > 0) {
  Q_covar_design_matrix = qr_thin_Q(covar_design_matrix) * sqrt(n_patients - 1);
  R_covar_design_matrix = qr_thin_R(covar_design_matrix) / sqrt(n_patients - 1);
}

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

array[n_cond_group + 1] int cond_group_pos = create_pos(cond_group_size);
