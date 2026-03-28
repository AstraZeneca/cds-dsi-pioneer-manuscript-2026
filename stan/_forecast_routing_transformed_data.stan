// ============================================================================
// Forecast Routing Transformed Data
// ============================================================================
// Computes parameter-sizing groups and position arrays for forecast patients.
//
// Requires in scope:
//   n_forecast_patients — either from data (full models) or as a transformed data
//                         constant (standalone: int n_forecast_patients = n_patients)
//   n_levels, n_groups_per_level — from _base_hierarchy_data.stan

// Parameter-sizing groups: mirrors n_groups_per_level but substitutes n_forecast_patients
// at the patient level (level n_levels). This ensures patient-level NCP parameters
// (tr_raw_level_intercept etc.) are sized by forecast patients only — background
// patients (when background) have no explicit parameters.
array[n_levels] int n_forecast_groups_per_level = n_groups_per_level;
n_forecast_groups_per_level[n_levels] = n_forecast_patients;

// Derived: total groups and position array for parameter-space indexing.
// Use n_forecast_level_pos (not level_pos) when indexing into parameter-sized arrays.
int n_forecast_total_groups = sum(n_forecast_groups_per_level);
array[n_levels + 1] int n_forecast_level_pos = create_pos(n_forecast_groups_per_level);
