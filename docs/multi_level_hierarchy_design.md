# Multi-Level Hierarchical Design for Stan Modules

**Date**: October 2, 2025\
**Status**: Design Proposal\
**Author**: System Design

## Table of Contents

1.  [Overview](#overview)
2.  [Core Concepts](#core-concepts)
3.  [Data Structure Design](#data-structure-design)
4.  [Transformed Data Implementation](#transformed-data-implementation)
5.  [Module Parameter Structure](#module-parameter-structure)
6.  [Transformed Parameters](#transformed-parameters)
7.  [Prior Specification](#prior-specification)
8.  [R Data Preparation](#r-data-preparation)
9.  [Migration Strategy](#migration-strategy)
10. [Examples](#examples)
11. [Utility Functions](#utility-functions)

------------------------------------------------------------------------

## Overview {#overview}

### Motivation

The current Stan modules have a two-level hierarchy hardcoded: - **Trial level**: Random effects across trials - **Patient level**: Random effects across patients within trials

This design proposal generalizes this to support **arbitrary multi-level hierarchies** where: - Each patient belongs to exactly one group at each level - Levels can be cross-classified (not strictly nested) - Patients are treated as the finest-grained level in the hierarchy - All levels use the same unified parameter structure

### Key Design Principles

1.  **Patients as Level N**: Treat patients as just another level (the last one) rather than special-casing them
2.  **Unified Parameter Structure**: All levels use identical parameter declarations and priors
3.  **Ragged Arrays**: Use position arrays to efficiently handle varying group sizes
4.  **Enable Flags**: Each level can be independently enabled/disabled for each module
5.  **Backward Compatibility**: Current trial+patient design maps naturally to this framework

------------------------------------------------------------------------

## Core Concepts {#core-concepts}

### Level Hierarchy Examples

#### Example 1: Trial → Patient (Current System)

```         
Level 1: Trial (n=5)
Level 2: Patient (n=150)
```

#### Example 2: Trial → Region → Patient

```         
Level 1: Trial (n=5)
Level 2: Region (n=10)
Level 3: Patient (n=150)
```

#### Example 3: Trial → Region → Site → Patient

```         
Level 1: Trial (n=5)
Level 2: Region (n=10)
Level 3: Site (n=30)
Level 4: Patient (n=150)
```

### Cross-Classified vs. Nested

This design supports **both** patterns:

-   **Nested**: Sites nested within regions, regions nested within trials
-   **Cross-Classified**: Patients belong to one trial AND one region (not necessarily nested)

The `patient_level_groups` matrix handles both cases identically.

------------------------------------------------------------------------

## Data Structure Design {#data-structure-design}

### File: `stan/base_data.stan`

#### Core Level Metadata

``` stan
// Multi-level hierarchy configuration
// Example: Level 1 = Trial, Level 2 = Region, Level 3 = Site, Level 4 = Patient
int<lower = 1> n_levels;  // Total number of levels INCLUDING patients

// Number of groups at each level
// Example: [5 trials, 10 regions, 30 sites, 150 patients]
array[n_levels] int<lower = 1> n_groups_per_level;

// Total groups across all levels
int<lower = 0> n_total_groups = sum(n_groups_per_level);
```

#### Patient Group Membership

``` stan
// For each patient (row), their group membership at each level (columns)
// Dimensions: [n_patients, n_levels]
// 
// Example for patient 42:
//   patient_level_groups[42, 1] = 3  (trial 3)
//   patient_level_groups[42, 2] = 7  (region 7)  
//   patient_level_groups[42, 3] = 15 (site 15)
//   patient_level_groups[42, 4] = 42 (patient 42 - identity mapping at finest level)
//
// Note: At the patient level (level n_levels), this is always the identity:
//       patient_level_groups[i, n_levels] = i
array[n_patients, n_levels] int patient_level_groups;
```

**Important Constraint**: The last level is always patients, where:

``` stan
patient_level_groups[i, n_levels] == i  // Identity mapping
n_groups_per_level[n_levels] == n_patients
```

------------------------------------------------------------------------

## Transformed Data Implementation {#transformed-data-implementation}

### File: `stan/base_transformed_data.stan`

#### Position Arrays

``` stan
// Position array for slicing by level
// level_pos[lv] gives the start index for level lv in flattened group arrays
array[n_levels + 1] int level_pos = create_pos(n_groups_per_level);

// Example with n_groups_per_level = [5, 10, 150]:
// level_pos = [1, 6, 16, 166]
//             |   |   |    |
//             |   |   |    End of all groups
//             |   |   Start of patient groups (16-165)
//             |   Start of region groups (6-15)
//             Start of trial groups (1-5)
```

#### Validation

``` stan
// Validate patient group assignments
for (p in 1:n_patients) {
  for (lv in 1:n_levels) {
    if (patient_level_groups[p, lv] < 1 || 
        patient_level_groups[p, lv] > n_groups_per_level[lv]) {
      reject("Patient ", p, " has invalid group assignment at level ", lv);
    }
  }
  
  // Verify patient level is identity mapping
  if (patient_level_groups[p, n_levels] != p) {
    reject("Patient ", p, " must have identity mapping at patient level (level ", n_levels, ")");
  }
}

// Check that patient level has n_patients groups
if (n_groups_per_level[n_levels] != n_patients) {
  reject("Patient level must have n_groups = n_patients");
}

print("Multi-level hierarchy validated:");
print("  n_levels = ", n_levels);
print("  n_groups_per_level = ", n_groups_per_level);
print("  n_total_groups = ", n_total_groups);
```

#### Optional: Patient Reordering (Advanced)

For cache efficiency, optionally pre-sort patients by group membership:

``` stan
// Create a reordered patient index for each level
// This groups patients by their level membership
array[n_levels] array[n_patients] int patients_sorted_by_level;

for (lv in 1:n_levels) {
  array[n_patients] int sort_keys = patient_level_groups[, lv];
  patients_sorted_by_level[lv] = sort_indices_by_key(sort_keys);
}

// Position arrays for sorted patients within each group
array[n_levels] array[n_groups_per_level[lv] + 1] int sorted_patient_pos_by_level;

for (lv in 1:n_levels) {
  array[n_groups_per_level[lv]] int n_patients_per_group = 
    zeros_int_array(n_groups_per_level[lv]);
  
  for (p in 1:n_patients) {
    int grp = patient_level_groups[p, lv];
    n_patients_per_group[grp] += 1;
  }
  
  sorted_patient_pos_by_level[lv] = create_pos(n_patients_per_group);
}
```

------------------------------------------------------------------------

## Module Parameter Structure {#module-parameter-structure}

### Module-Specific Flags (Data Block)

Each module (tr, frac, init) has its own enable flags:

``` stan
// Example for Total Rate (tr) module
// File: stan/ssls/modules/tr/flags.stan

// Population-level (above hierarchy)
int<lower=0, upper=1> enable_pop_cov_tr;

// Per-level flags (one for each hierarchical level)
array[n_levels] int<lower=0, upper=1> enable_level_intercept_tr;
array[n_levels] int<lower=0, upper=1> enable_level_cov_tr;
```

**Example Configuration**:

``` r
# Enable trial and patient intercepts, but not region:
enable_level_intercept_tr <- c(1, 0, 1)  # trial=YES, region=NO, patient=YES
```

### Unified Parameter Declarations

### File: `stan/ssls/modules/tr/parameters.stan`

``` stan
// tr/parameters.stan — fully unified across all levels

// ===== POPULATION LEVEL (Level 0, above all hierarchy) =====
real tr_loc_pop;  // Population intercept

// Population covariate coefficients (QR space)
vector[enable_pop_cov_tr ? n_covar : 0] tr_coef_qr_pop;

// ===== HIERARCHICAL STRUCTURE (All Levels Including Patients) =====

// Random intercepts: one SD hyperparameter per level
array[n_levels] real<lower=0> tr_sd_level_intercept;

// Raw standard normal draws for intercepts
// This is a ragged array: concatenate all groups across all levels
// Size: n_total_groups if all enabled, or sum of enabled groups
vector[n_total_groups] tr_raw_level_intercept;

// Random slopes: one SD vector per level
array[n_levels] vector<lower=0>[n_covar] tr_sd_level_slope;

// Raw standard normal draws for slopes  
// Dimensions: [n_total_groups, n_covar]
matrix[n_total_groups, n_covar] tr_raw_level_slope;
```

**Note**: If using enable flags to conditionally size parameters:

``` stan
// Compute total enabled groups
int n_total_groups_enabled_tr = 0;
for (lv in 1:n_levels) {
  if (enable_level_intercept_tr[lv]) {
    n_total_groups_enabled_tr += n_groups_per_level[lv];
  }
}

// Parameters only for enabled groups
vector[n_total_groups_enabled_tr] tr_raw_level_intercept;
```

### File: `stan/ssls/modules/frac/parameters.stan`

``` stan
// frac/parameters.stan — same structure as tr module

real frac_logit_loc_pop;  // Population intercept (logit scale)
vector[enable_pop_cov_frac ? n_covar : 0] frac_coef_qr_pop;

array[n_levels] real<lower=0> frac_sd_level_intercept;
vector[n_total_groups] frac_raw_level_intercept;

array[n_levels] vector<lower=0>[n_covar] frac_sd_level_slope;
matrix[n_total_groups, n_covar] frac_raw_level_slope;
```

### File: `stan/ssls/modules/init/parameters.stan`

``` stan
// init/parameters.stan — same structure

real init_logit_loc_pop;
vector[enable_pop_cov_init ? n_covar : 0] init_coef_qr_pop;

array[n_levels] real<lower=0> init_sd_level_intercept;
vector[n_total_groups] init_raw_level_intercept;

array[n_levels] vector<lower=0>[n_covar] init_sd_level_slope;
matrix[n_total_groups, n_covar] init_raw_level_slope;
```

------------------------------------------------------------------------

## Transformed Parameters {#transformed-parameters}

### Building Linear Predictors

### File: `stan/ssls/modules/tr/transformed_parameters.stan`

``` stan
// tr/transformed_parameters.stan — unified linear predictor construction

// ===== INTERCEPT EFFECTS =====
// Build linear predictor by summing across ALL levels
vector[n_patients] tr_linpred_levels = rep_vector(0, n_patients);

// Scale and accumulate effects from each level
for (lv in 1:n_levels) {
  if (enable_level_intercept_tr[lv]) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(level_pos, lv);
    
    // Scale raw effects by SD
    vector[n_groups_per_level[lv]] level_effects = 
      tr_sd_level_intercept[lv] * tr_raw_level_intercept[lv_start:lv_end];
    
    // Add to each patient's predictor via their group membership
    tr_linpred_levels += level_effects[patient_level_groups[, lv]];
  }
}

// ===== COVARIATE SLOPE EFFECTS =====
vector[n_patients] tr_linpred_slopes = rep_vector(0, n_patients);

for (lv in 1:n_levels) {
  if (enable_level_cov_tr[lv] && n_covar > 0) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(level_pos, lv);
    
    // Scale raw slopes for this level
    matrix[n_groups_per_level[lv], n_covar] level_slopes_qr = 
      tr_raw_level_slope[lv_start:lv_end, :] .* 
      rep_matrix(tr_sd_level_slope[lv]', n_groups_per_level[lv]);
    
    // Accumulate slope contributions
    for (p in 1:n_patients) {
      int grp = patient_level_groups[p, lv];
      tr_linpred_slopes[p] += dot_product(
        Q_covar_design_matrix[p, :], 
        level_slopes_qr[grp, :]
      );
    }
  }
}

// ===== POPULATION COVARIATE EFFECTS =====
vector[n_patients] tr_linpred_pop = enable_pop_cov_tr ? 
  (Q_covar_design_matrix * tr_coef_qr_pop) : 
  rep_vector(0, n_patients);

// ===== FINAL LINEAR PREDICTOR =====
vector[n_patients] tr_loc_patient = 
  tr_loc_pop 
  + tr_linpred_pop
  + tr_linpred_levels 
  + tr_linpred_slopes;
```

**Key Points**: - Single loop handles all levels uniformly - Patient level (level n_levels) is handled identically to other levels - No special-casing of patient parameters

------------------------------------------------------------------------

## Prior Specification {#prior-specification}

### Hyperparameters

### File: `stan/ssls/modules/tr/hyperparams.stan`

``` stan
// tr/hyperparams.stan — prior hyperparameters

// Population priors
real tr_loc_pop_mean;
real<lower=0> tr_loc_pop_sd;

// Unified across all levels (including patients)
// One hyperparameter per level
array[n_levels] real<lower=0> tr_sd_level_intercept_sd;
array[n_levels] real<lower=0> tr_sd_level_slope_sd;

// Optional: population covariate slope SD
real<lower=0> tr_coef_pop_sd;  // If not using QR (otherwise implicitly 1)
```

**Example Values**:

``` r
# For 3 levels: trial, region, patient
tr_sd_level_intercept_sd <- c(
  0.5,  # Trial level: moderate variation
  0.3,  # Region level: smaller variation
  0.2   # Patient level: residual variation
)
```

### Prior Distributions

### File: `stan/ssls/modules/tr/priors.stan`

``` stan
// tr/priors.stan — unified priors across all levels

// Population intercept
tr_loc_pop ~ normal(tr_loc_pop_mean, tr_loc_pop_sd);

// Population covariate effects
if (enable_pop_cov_tr) {
  tr_coef_qr_pop ~ normal(0, 1);  // QR decomposition implies unit scale
}

// ===== UNIFIED LOOP OVER ALL LEVELS =====
for (lv in 1:n_levels) {
  // Intercept SD hyperprior (always specified)
  tr_sd_level_intercept[lv] ~ normal(0, tr_sd_level_intercept_sd[lv]);
  
  // Intercept raw effects (if enabled)
  if (enable_level_intercept_tr[lv]) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(level_pos, lv);
    tr_raw_level_intercept[lv_start:lv_end] ~ std_normal();
  }
  
  // Slope SD hyperprior (if covariates enabled)
  if (enable_level_cov_tr[lv] && n_covar > 0) {
    tr_sd_level_slope[lv] ~ normal(0, tr_sd_level_slope_sd[lv]);
    
    // Slope raw effects
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(level_pos, lv);
    to_vector(tr_raw_level_slope[lv_start:lv_end, :]) ~ std_normal();
  }
}
```

**Benefits**: - Single loop handles all levels (trial, region, patient, etc.) - No code duplication - Easy to add/remove levels

------------------------------------------------------------------------

## R Data Preparation {#r-data-preparation}

### Example 1: Current System (Trial + Patient)

``` r
# Backward-compatible with current code
n_trials <- 3
n_patients <- 150

# Current: patient_trial assignment
patient_trial <- sample(1:n_trials, n_patients, replace = TRUE)

# Convert to multi-level format
n_levels <- 2  # trial (level 1) + patient (level 2)

patient_level_groups <- cbind(
  patient_trial,      # Level 1: trial membership
  1:n_patients       # Level 2: patient identity
)

stan_data <- list(
  # Core hierarchy
  n_levels = n_levels,
  n_groups_per_level = c(n_trials, n_patients),
  patient_level_groups = patient_level_groups,
  
  # Enable flags for tr module
  enable_pop_cov_tr = 1,
  enable_level_intercept_tr = c(1, 1),  # Trial and patient random effects
  enable_level_cov_tr = c(1, 1),        # Both levels have random slopes
  
  # Hyperparameters
  tr_loc_pop_mean = 0,
  tr_loc_pop_sd = 1,
  tr_sd_level_intercept_sd = c(0.5, 0.2),  # Trial SD, Patient SD
  tr_sd_level_slope_sd = c(0.3, 0.1)       # Trial slope SD, Patient slope SD
)
```

### Example 2: Trial + Region + Patient

``` r
# Three-level hierarchy
n_trials <- 5
n_regions <- 10  
n_patients <- 150

# Assign patients to trials and regions (cross-classified)
patient_trial <- sample(1:n_trials, n_patients, replace = TRUE)
patient_region <- sample(1:n_regions, n_patients, replace = TRUE)

n_levels <- 3  # trial + region + patient

patient_level_groups <- cbind(
  patient_trial,      # Level 1: trial
  patient_region,     # Level 2: region
  1:n_patients       # Level 3: patient (identity)
)

stan_data <- list(
  n_levels = 3,
  n_groups_per_level = c(n_trials, n_regions, n_patients),
  patient_level_groups = patient_level_groups,
  
  # Enable trial and patient, but NOT region intercepts
  enable_level_intercept_tr = c(1, 0, 1),
  
  # Hyperparameters for 3 levels
  tr_sd_level_intercept_sd = c(0.5, 0.3, 0.2)
)
```

### Example 3: Nested Hierarchy (Trial → Region → Site → Patient)

``` r
# Four-level nested hierarchy
n_trials <- 5
n_regions_per_trial <- 2
n_sites_per_region <- 3
n_patients_per_site <- 10

# Generate nested structure
trial_data <- expand.grid(
  site_in_region = 1:n_sites_per_region,
  region_in_trial = 1:n_regions_per_trial,
  trial = 1:n_trials
)

trial_data$region <- (trial_data$trial - 1) * n_regions_per_trial + 
                     trial_data$region_in_trial
trial_data$site <- (trial_data$region - 1) * n_sites_per_region + 
                   trial_data$site_in_region

# Replicate for patients
patient_data <- trial_data[rep(1:nrow(trial_data), each = n_patients_per_site), ]
patient_data$patient <- 1:nrow(patient_data)

n_patients <- nrow(patient_data)
n_sites <- max(patient_data$site)
n_regions <- max(patient_data$region)

n_levels <- 4

patient_level_groups <- cbind(
  patient_data$trial,    # Level 1
  patient_data$region,   # Level 2
  patient_data$site,     # Level 3
  patient_data$patient   # Level 4 (identity)
)

stan_data <- list(
  n_levels = 4,
  n_groups_per_level = c(n_trials, n_regions, n_sites, n_patients),
  patient_level_groups = patient_level_groups,
  
  enable_level_intercept_tr = c(1, 1, 1, 1),  # All levels enabled
  tr_sd_level_intercept_sd = c(0.5, 0.3, 0.2, 0.15)  # Decreasing variation
)
```

### Helper Functions

Create helper functions in R for data preparation:

``` r
# File: r/multi_level_hierarchy.R

#' Create multi-level hierarchy structure for Stan
#' 
#' @param patient_assignments List of vectors, one per level (excluding patient level)
#' @param n_patients Number of patients
#' @return List with hierarchy structure for Stan
create_hierarchy_structure <- function(patient_assignments, n_patients) {
  n_levels <- length(patient_assignments) + 1  # +1 for patient level
  
  # Build patient_level_groups matrix
  patient_level_groups <- do.call(cbind, patient_assignments)
  patient_level_groups <- cbind(patient_level_groups, 1:n_patients)
  
  # Calculate n_groups_per_level
  n_groups_per_level <- c(
    sapply(patient_assignments, max),
    n_patients
  )
  
  list(
    n_levels = n_levels,
    n_groups_per_level = n_groups_per_level,
    patient_level_groups = patient_level_groups
  )
}

# Usage:
hierarchy <- create_hierarchy_structure(
  patient_assignments = list(
    trial = patient_trial,
    region = patient_region
  ),
  n_patients = n_patients
)

stan_data <- c(stan_data, hierarchy)
```

------------------------------------------------------------------------

## Migration Strategy {#migration-strategy}

### Phase 1: Infrastructure (Week 1)

1.  **Update `pos.stan`**: Add new utility functions
    -   `get_level_pos()`
    -   `get_global_group_idx()`
    -   `get_level_slice()`
2.  **Update `base_data.stan`**: Add multi-level structures
    -   `n_levels`
    -   `n_groups_per_level`
    -   `patient_level_groups`
    -   Keep `patient_trial` for backward compatibility (alias it in transformed data)
3.  **Update `base_transformed_data.stan`**: Add level validation and position arrays
    -   Create `level_pos`
    -   Validate hierarchy
    -   Optional: Create `trial_patient_pos` from `patient_level_groups[, 1]` for backward compatibility

### Phase 2: Single Module Migration (Week 2)

1.  **Migrate `tr` module** as proof of concept:
    -   Update `parameters.stan` to use level arrays
    -   Update `transformed_parameters.stan` to loop over levels
    -   Update `priors.stan` to loop over levels
    -   Update `hyperparams.stan` to use level arrays
2.  **Test with backward-compatible data**:
    -   Use `n_levels = 2` (trial + patient)
    -   Verify identical results to current implementation

### Phase 3: Remaining Modules (Week 3)

1.  Migrate `frac` module
2.  Migrate `init` module
3.  Update any legacy code that references `patient_trial` directly

### Phase 4: R Infrastructure (Week 4)

1.  Create helper functions in `r/multi_level_hierarchy.R`
2.  Update existing analysis scripts to use new format
3.  Create examples and documentation

### Phase 5: Testing & Validation (Week 5)

1.  Test with multi-level hierarchies (trial + region + patient)
2.  Validate against current production results
3.  Performance benchmarking
4.  Documentation updates

### Backward Compatibility Approach

In `base_transformed_data.stan`, maintain compatibility:

``` stan
// For backward compatibility with existing code
array[n_patients] int patient_trial;

if (n_levels >= 1) {
  patient_trial = patient_level_groups[, 1];
} else {
  patient_trial = rep_array(1, n_patients);  // Single trial
}

// Existing code can still use patient_trial
array[n_trials + 1] int trial_patient_pos = create_pos(n_trial_patients);
// ... etc
```

------------------------------------------------------------------------

## Examples {#examples}

### Example A: Simple Two-Level (Current System)

**Data**:

```         
n_levels = 2
n_groups_per_level = [3, 10]  # 3 trials, 10 patients
patient_level_groups = [
  [1, 1],  # Patient 1: trial 1
  [1, 2],  # Patient 2: trial 1  
  [1, 3],  # Patient 3: trial 1
  [2, 4],  # Patient 4: trial 2
  [2, 5],  # Patient 5: trial 2
  [2, 6],  # Patient 6: trial 2
  [3, 7],  # Patient 7: trial 3
  [3, 8],  # Patient 8: trial 3
  [3, 9],  # Patient 9: trial 3
  [3, 10]  # Patient 10: trial 3
]
```

**Parameters** (if all enabled):

```         
tr_raw_level_intercept = [
  trial1_effect,      # Index 1
  trial2_effect,      # Index 2
  trial3_effect,      # Index 3
  patient1_effect,    # Index 4
  patient2_effect,    # Index 5
  ...
  patient10_effect    # Index 13
]

level_pos = [1, 4, 14]
```

**Linear Predictor for Patient 4**:

``` stan
Level 1 (trial): tr_sd_level_intercept[1] * tr_raw_level_intercept[2]  // trial 2
Level 2 (patient): tr_sd_level_intercept[2] * tr_raw_level_intercept[7]  // patient 4 (index 4-1+4=7)
```

### Example B: Three-Level Cross-Classified

**Data**:

```         
n_levels = 3
n_groups_per_level = [2, 3, 10]  # 2 trials, 3 regions, 10 patients
patient_level_groups = [
  [1, 1, 1],   # Patient 1: trial 1, region 1
  [1, 1, 2],   # Patient 2: trial 1, region 1
  [1, 2, 3],   # Patient 3: trial 1, region 2
  [1, 3, 4],   # Patient 4: trial 1, region 3
  [2, 1, 5],   # Patient 5: trial 2, region 1 (cross-classified!)
  [2, 2, 6],   # Patient 6: trial 2, region 2
  [2, 2, 7],   # Patient 7: trial 2, region 2
  [2, 3, 8],   # Patient 8: trial 2, region 3
  [2, 3, 9],   # Patient 9: trial 2, region 3
  [2, 3, 10]   # Patient 10: trial 2, region 3
]
```

**Parameters**:

```         
tr_raw_level_intercept = [
  trial1, trial2,           # Indices 1-2
  region1, region2, region3, # Indices 3-5
  patient1, ..., patient10   # Indices 6-15
]

level_pos = [1, 3, 6, 16]
```

**Linear Predictor for Patient 5**:

``` stan
Level 1 (trial): effects[2]      # Trial 2
Level 2 (region): effects[3]     # Region 1  
Level 3 (patient): effects[10]   # Patient 5
```

------------------------------------------------------------------------

## Utility Functions {#utility-functions}

### New Functions for `pos.stan`

``` stan
/**
 * Get the global group index from level and local group ID
 * 
 * @param level_pos Position array for levels
 * @param level Level index (1-based)
 * @param group_id Local group ID within level (1-based)
 * @return Global group index in flattened group array
 */
int get_global_group_idx(array[] int level_pos, int level, int group_id) {
  return level_pos[level] + group_id - 1;
}

/**
 * Extract all groups for a given level from flattened array
 * 
 * @param values Flattened array of values across all levels/groups
 * @param level_pos Position array for levels
 * @param level Level index to extract
 * @return Sub-array for specified level
 */
vector get_level_slice(vector values, array[] int level_pos, int level) {
  int start, end;
  (start, end) = get_pos(level_pos, level);
  return values[start:end];
}

/**
 * Extract patient indices belonging to a specific group at a specific level
 * 
 * @param patient_level_groups Patient group membership matrix
 * @param level Level index
 * @param group Group ID at that level
 * @return Array of patient indices in this group
 */
array[] int get_patients_in_group(
  array[,] int patient_level_groups,
  int level,
  int group
) {
  int n_patients = dims(patient_level_groups)[1];
  array[n_patients] int mask;
  int count = 0;
  
  for (p in 1:n_patients) {
    if (patient_level_groups[p, level] == group) {
      count += 1;
      mask[count] = p;
    }
  }
  
  return mask[1:count];
}

/**
 * Create position array for multi-level hierarchy
 * Same as standard create_pos but documented for multi-level use
 */
array[] int create_level_pos(array[] int n_groups_per_level) {
  return create_pos(n_groups_per_level);
}
```

------------------------------------------------------------------------

## Benefits Summary

### Code Benefits

1.  **Unified Structure**: All levels use identical parameter and prior patterns
2.  **No Special Cases**: Patient level handled like any other level
3.  **Flexible**: Easy to add/remove levels
4.  **Maintainable**: Single code path for all hierarchical effects
5.  **Efficient**: Ragged arrays minimize memory usage

### Statistical Benefits

1.  **Arbitrary Hierarchies**: Support any number of levels
2.  **Cross-Classified**: Not limited to nested structures
3.  **Selective Modeling**: Enable/disable levels independently per module
4.  **Flexible Priors**: Different hyperparameters per level
5.  **Extensible**: Easy to add new grouping factors

### Practical Benefits

1.  **Backward Compatible**: Current trial+patient maps naturally
2.  **Clear Semantics**: Level structure explicit in data
3.  **Easy Testing**: Can validate with simple hierarchies first
4.  **Reusable**: Pattern works for all modules (tr, frac, init, etc.)
5.  **Documented**: Clear examples and migration path

------------------------------------------------------------------------

## Open Questions & Future Enhancements

### Questions to Resolve

1.  Should we support different n_levels per module? (e.g., tr uses 3 levels, frac uses 2)
2.  Do we need crossed random effects within a level? (e.g., trial × treatment interaction)
3.  Should hyperpriors be hierarchical? (e.g., SD of SDs across levels)

### Potential Enhancements

1.  **Adaptive Levels**: Automatically detect levels from data structure
2.  **Level Names**: Store level names for better output labeling (via metadata file)
3.  **Sparse Representation**: For levels with many zero-membership groups
4.  **Interaction Terms**: Cross-level interactions (trial × region effects)
5.  **Time-Varying Hierarchy**: Patients switching groups over time

------------------------------------------------------------------------

## References

### Related Documentation

-   [Current Parameter Naming Convention](./naming_convention_compliance_check.md)
-   [Module Migration Summary](./parameter_renaming_summary.md)
-   [Position Array Documentation](../stan/pos.stan)

### Stan Documentation

-   [Hierarchical Models](https://mc-stan.org/docs/stan-users-guide/hierarchical-models.html)
-   [Ragged Data Structures](https://mc-stan.org/docs/stan-users-guide/ragged-data-structures.html)
-   [QR Decomposition for Regression](https://mc-stan.org/docs/stan-users-guide/QR-reparameterization.html)

------------------------------------------------------------------------

**Document Version**: 1.0\
**Last Updated**: October 2, 2025\
**Status**: Proposed Design - Pending Implementation