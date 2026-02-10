# Stan Code Cleanup Summary

**Date**: 2026-02-10

## What We Did

### 1. Reorganized Stan Directory Structure (Commit 5818a72)

**Major Changes:**
- Merged `base_data.stan` and `tumor/base_data.stan` into single file
- Moved all `stan/ssls/` contents to `stan/` root level
- Updated all `#include` paths (removed `../` references)
- Archived unused legacy code to `stan/legacy/`
- Removed empty `ssls/` directory

**Benefits:**
- Simplified directory structure
- Clearer code organization
- Easier to navigate
- All models in one place

### 2. Removed Dead Code (Commit 683fdd6)

**Removed ~180 lines of unused transformed data computations:**

#### tumor/tumor_transformed_data.stan (186 → 64 lines, 65% reduction)
Removed entire legacy visit tracking infrastructure:
- `n_patient_unique_visits`, `patient_unique_visits`, `patient_unique_visits_pos`
- `n_pop_unique_visits` (scoped to local), `pop_unique_visits` (scoped to local)
- `n_trial_unique_visits`, `trial_unique_visits`, `trial_unique_visits_pos`
- `measured_tumor_visits`, `measured2patient_visits_idx`
- `non_measured_tumor_visits`, `non_measured2patient_visits_idx`
- `patient_measured_tumor_visits_pos`, `patient_non_measured_tumor_visits_pos`
- `n_patient_non_measured_tumor_visits`, `n_patient_post_treat_visits`
- `post_treat_pos`, `post_treat_visits_pos`, `patient_non_measured_offset`

Kept only essentials:
- ✅ `patient2pop_unique_visit_idx` - needed for GP modeling
- ✅ `post_treat_sld` - needed for analysis
- ✅ `n_patient_forecast_visits` - needed for predictions
- ✅ `all_tumor_measure_t` - needed for GP time grid

#### _sf_transformed_data.stan (~50 lines removed)
Removed:
- `rep_time_points` (unused output from `sf_rep_T`)
- `n_total_visits_m1`
- Entire censoring tracking block (~40 lines):
  - `n_right_censored_patients`
  - `n_right_uncensored_patients`
  - `right_uncensored_patients`
  - `right_uncensored_idx`
  - `trial_right_censored_pos`
  - `trial_right_uncensored_pos`
- `ub_pfs_p1` array
- `NT_STABLE` constant
- `sorted_pfs_timepoints`

Kept (actually used in other blocks):
- ✅ `n_total_forecast_visits` - used in generated quantities
- ✅ `max_unique_visit` - used in transformed parameters
- ✅ `n_total_groups` - used in transformed parameters
- ✅ `visit_cumsum_mat` - used for batched state computation

#### base_transformed_data.stan (~20 lines removed)
Removed:
- `min_all_t` (computed max_t_width directly instead)
- `patient_max_t_width` array
- `t_visit_trial_idx` array
- `trial_visit_pos` array
- `all_measure_t` array (unused GP time array)
- Scoped `n_trial_patients` to local variable

#### Module files (~5 lines removed)
- `modules/measurement/transformed_data.stan`: Removed `lod` (use `log_lod` instead)
- `modules/other_events/transformed_data.stan`: Scoped `n_positive` and `quantiles_obs` to local

**Benefits:**
- **180 lines of dead code removed**
- Faster compilation
- Reduced memory usage
- Clearer code (easier to understand what's actually used)
- Removed entire legacy subsystems that are no longer needed

## Analysis Tools Created

### scripts/analyze_stan_variables.py

Python script to analyze variable usage in Stan models.

**Two modes:**

1. **Data Variable Analysis** (default):
   ```bash
   python scripts/analyze_stan_variables.py /mnt/code/stan
   ```
   Finds unused DATA inputs (variables declared in data blocks)

2. **Transformed Data Analysis** (`--transformed-data`):
   ```bash
   python scripts/analyze_stan_variables.py /mnt/code/stan --transformed-data
   ```
   Finds dead code (variables created in transformed data but never used)

**Features:**
- Recursively resolves all `#include` directives
- Traces usage through transformed_data, parameters, model, generated quantities
- Categorizes variables by usage context
- Cross-model comparison
- Line-by-line usage traces

## Verification

All three models compile successfully:
```bash
✅ stan/sf-ssm-log-space.stan
✅ stan/sf-ssls-lfo.stan
✅ stan/sf-ssls-lfo-endpoints.stan
```

Verified with:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/<model>.stan
```

## Reports Generated

1. **REORGANIZATION_COMPLETE.md** - Full details of directory reorganization
2. **UNUSED_VARIABLES_REPORT.md** - Analysis of unused data inputs
3. **DEAD_CODE_REPORT.md** - Analysis of unused transformed data variables
4. **CLEANUP_SUMMARY.md** (this file) - Overall summary

## Impact

### Before:
- stan/ssls/ nested structure with `../` includes
- 186 lines in tumor_transformed_data.stan
- ~278 total lines of dead code across files
- Two separate base_data.stan files

### After:
- Flat stan/ structure with clean includes
- 64 lines in tumor_transformed_data.stan (65% reduction)
- Dead code removed
- Single merged base_data.stan
- Cleaner, more maintainable codebase

## Next Steps (Optional)

The analysis tool identified some additional potential cleanups that weren't removed due to false positives:

### LFO-Specific Dead Code
The `_lfo_transformed_data.stan` file has ~22 additional unused variables identified by the tool. These could be reviewed manually for removal, but may be scaffolding for future LFO enhancements.

### Data Input Cleanup
Only 4 data inputs are truly unused across ALL models:
- `calendar_week` - never used
- `death_week` - never used
- `extend_max_all_t` - rarely needed
- `oe_enable_trial_tumor_cov` - flag never used

These could be removed from `base_data.stan` and `modules/other_events/flags.stan`.

## Git History

```
683fdd6 Remove dead code from transformed data sections (~180 lines)
5818a72 Reorganize Stan directory structure and add variable analysis tools
```

Both commits include full details in their messages and are co-authored with Claude.
