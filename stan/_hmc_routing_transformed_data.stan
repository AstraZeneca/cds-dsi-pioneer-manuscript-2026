// ============================================================================
// HMC Routing Transformed Data
// ============================================================================
// Computes parameter-sizing groups and position arrays for HMC patients.
//
// Requires in scope:
//   n_hmc_patients — either from data (full models) or as a transformed data
//                    constant (standalone: int n_hmc_patients = n_patients)
//   n_levels, n_groups_per_level — from _base_hierarchy_data.stan

// Parameter-sizing groups: mirrors n_groups_per_level but substitutes n_hmc_patients
// at the patient level (level n_levels). This ensures patient-level NCP parameters
// (tr_raw_level_intercept etc.) are sized by HMC patients only — Laplace-marginalized
// patients have no explicit parameters and must not inflate the HMC parameter space.
array[n_levels] int n_hmc_groups_per_level = n_groups_per_level;
n_hmc_groups_per_level[n_levels] = n_hmc_patients;

// Derived: total groups and position array for parameter-space indexing.
// Use n_hmc_level_pos (not level_pos) when indexing into parameter-sized arrays.
int n_hmc_total_groups = sum(n_hmc_groups_per_level);
array[n_levels + 1] int n_hmc_level_pos = create_pos(n_hmc_groups_per_level);
