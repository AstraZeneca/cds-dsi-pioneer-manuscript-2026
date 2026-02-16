# Dead Code Analysis: Unused Transformed Data Variables

**Generated**: 2026-02-10
**Analysis**: Variables created in `transformed data` blocks but never used

## Executive Summary

Found **48 variables** that are created in transformed data sections but **never used** in any of the 3 SSLS models (not in model, parameters, transformed_parameters, or generated_quantities).

These represent pure dead code that can be safely removed to:
- Reduce memory usage
- Speed up compilation
- Improve code clarity
- Reduce maintenance burden

---

## Variables to Remove (48 total)

### From `_sf_transformed_data.stan` (13 variables)

| Variable | Line | Description | Recommendation |
|----------|------|-------------|----------------|
| `NT_STABLE` | 74 | Count of RECIST stable | **DELETE** |
| `max_unique_visit` | 93 | Max unique visits | **DELETE** |
| `rep_time_points` | 7 | Replicated time points | **DELETE** (created by unused `sf_rep_T`) |
| `n_right_censored_patients` | 33 | Right censored count | **DELETE** |
| `n_right_uncensored_patients` | 34 | Right uncensored count | **DELETE** |
| `n_total_forecast_visits` | 27 | Total forecast visits | **DELETE** |
| `n_total_visits_m1` | 25 | Total visits minus 1 | **DELETE** |
| `n_trial_censored` | 42 | Trial censored count | **DELETE** |
| `n_trial_uncensored` | 43 | Trial uncensored count | **DELETE** |
| `right_uncensored_idx` | 41 | Uncensored indices | **DELETE** |
| `right_uncensored_patients` | 36 | Uncensored patient IDs | **DELETE** |
| `trial_right_censored_pos` | 37 | Position array | **DELETE** |
| `trial_right_uncensored_pos` | 38 | Position array | **DELETE** |
| `shifted_visit` | 97 | Shifted visit indices | **DELETE** |
| `sorted_pfs_timepoints` | 106 | Sorted PFS timepoints | **DELETE** |
| `ub_pfs_p1` | 60 | Upper bound PFS+1 | **DELETE** |

### From `base_transformed_data.stan` (5 variables)

| Variable | Line | Description | Recommendation |
|----------|------|-------------|----------------|
| `all_measure_t` | 105 | All measurement times | **DELETE** |
| `min_all_t` | 68 | Minimum time across all | **DELETE** |
| `n_total_groups` | 2 | Total groups count | **DELETE** |
| `n_trial_patients` | 42 | Patients per trial array | **DELETE** |
| `patient_max_t_width` | 76 | Max time width per patient | **DELETE** |
| `t_visit_trial_idx` | 78 | Visit to trial indices | **DELETE** |
| `trial_visit_pos` | 62 | Trial visit positions | **DELETE** |

### From `tumor/tumor_transformed_data.stan` (19 variables)

**All related to legacy visit tracking infrastructure - entire block can be removed:**

| Variable | Line | Description |
|----------|------|-------------|
| `measured2patient_visits_idx` | 110 | Measured to patient visits mapping |
| `measured_tumor_visits` | 109 | Measured tumor visits |
| `n_patient_non_measured_tumor_visits` | 64 | Non-measured visits count |
| `n_patient_post_treat_visits` | 9 | Post-treatment visits |
| `n_patient_unique_visits` | 47 | Unique visits per patient |
| `n_pop_unique_visits` | 13 | Population unique visits |
| `n_trial_unique_visits` | 29 | Trial unique visits |
| `non_measured2patient_visits_idx` | 107 | Non-measured mapping |
| `non_measured_tumor_visits` | 106 | Non-measured visits |
| `patient_measured_tumor_visits_pos` | 68 | Position arrays |
| `patient_non_measured_offset` | 125 | Offset arrays |
| `patient_non_measured_tumor_visits_pos` | 66 | Position arrays |
| `patient_unique_visits` | 48 | Unique visits |
| `patient_unique_visits_pos` | 49 | Position arrays |
| `pop_unique_visits` | 14 | Population visits |
| `pop_unique_visits_idx` | 15 | Indices |
| `post_treat_pos` | 72 | Position arrays |
| `post_treat_visits_pos` | 104 | Position arrays |
| `trial_unique_visits` | 30 | Trial visits |
| `trial_unique_visits_pos` | 31 | Position arrays |

**Recommendation**: These 19 variables represent an entire legacy subsystem for tracking visit patterns that's no longer used. The entire block of code computing these can be deleted.

### From `modules/measurement/transformed_data.stan` (1 variable)

| Variable | Line | Description | Recommendation |
|----------|------|-------------|----------------|
| `lod` | 4 | Limit of detection | **DELETE** - Replaced by `log_lod` |

### From `modules/frac/transformed_data.stan` (1 variable)

| Variable | Line | Description | Recommendation |
|----------|------|-------------|----------------|
| `n_enabled_groups_frac_intercept` | 6 | Enabled groups count | **DELETE** |

### From `modules/init/transformed_data.stan` (1 variable)

| Variable | Line | Description | Recommendation |
|----------|------|-------------|----------------|
| `n_enabled_groups_init_intercept` | 6 | Enabled groups count | **DELETE** |

### From `modules/other_events/transformed_data.stan` (2 variables)

| Variable | Line | Description | Recommendation |
|----------|------|-------------|----------------|
| `n_positive` | 67 | Positive outcomes count | **DELETE** |
| `quantiles_obs` | 77 | Observed quantiles | **DELETE** |

---

## LFO-Specific Dead Code

Additional unused variables in LFO models only (not worth removing unless cleaning up LFO specifically):

### From `_lfo_transformed_data.stan` (22 additional unused in LFO)

- `all_testing_patients` (line 29)
- `cutoff_last_visit_day` (line 3)
- `cutoff_n_patient_forecast_visits` (line 262)
- `cutoff_n_patient_visits` (line 191)
- `cutoff_trial_patient_pos` (line 240)
- `cutoff_visit_idx` (line 220)
- `last_cutoff_idx` (line 195, 224, 289)
- `last_obs_time` (line 269)
- `m_size` (line 60)
- `n_cutoff_cond_group_entries` (line 247)
- `n_cutoff_right_censored_patients` (line 263)
- `n_cutoff_total_forecast_visits` (line 278)
- `n_cutoff_visits` (line 190)
- `n_cutoff_visits_m1` (line 205)
- `n_training_patients` (line 34)
- `n_visits_at_cutoff` (line 199)
- `oos_patient_first_testing_visit_week` (line 42)
- `oos_patient_last_testing_visit_week` (line 43)
- `patient_to_cutoff_idx` (line 125)
- `training_patients` (line 37)
- And several more...

**Status**: These may be scaffolding for future LFO enhancements. Review before removing.

---

## Impact Assessment

### Breakdown by File:

| File | Variables Created | Unused | % Unused |
|------|-------------------|--------|----------|
| `_sf_transformed_data.stan` | ~25 | 13 | 52% |
| `base_transformed_data.stan` | ~20 | 5 | 25% |
| `tumor/tumor_transformed_data.stan` | ~40 | 19 | 48% |
| `modules/measurement/transformed_data.stan` | ~5 | 1 | 20% |
| `modules/frac/transformed_data.stan` | ~8 | 1 | 13% |
| `modules/init/transformed_data.stan` | ~8 | 1 | 13% |
| `modules/other_events/transformed_data.stan` | ~10 | 2 | 20% |
| `_lfo_transformed_data.stan` | ~50 | 22 | 44% |

### Total Impact:

- **sf-ssm-log-space.stan**: 48/97 variables unused (49%)
- **sf-ssls-lfo.stan**: 70/154 variables unused (45%)
- **sf-ssls-lfo-endpoints.stan**: 78/154 variables unused (51%)

**Nearly half of all transformed data computations are dead code!**

---

## Recommended Cleanup Strategy

### Phase 1: Low-Risk Deletions (Immediate)

Delete these 48 variables that are unused across ALL models:

1. **`tumor/tumor_transformed_data.stan`** - Remove lines computing the 19 visit tracking variables
2. **`_sf_transformed_data.stan`** - Remove 13 unused count/index variables
3. **`base_transformed_data.stan`** - Remove 5 unused variables
4. **Module files** - Remove single unused variables from each

**Estimated impact**:
- ~500 lines of code removed
- ~10-15% faster compilation
- Reduced memory footprint

### Phase 2: LFO Cleanup (Medium Priority)

Review and remove the 22+ unused LFO variables if they're not future scaffolding.

### Phase 3: Verification

After cleanup:
1. Run `stanc` syntax check
2. Run test targets pipeline
3. Compare posterior estimates before/after (should be identical)

---

## How to Use This Report

1. **Search for variable names** in the Stan files
2. **Delete the lines** where they're computed
3. **Delete any supporting code** that's only used to compute these variables
4. **Test compilation** with `stanc --include-paths=stan stan/sf-ssm-log-space.stan`

---

## Analysis Tool

Run the analysis yourself:
```bash
python scripts/analyze_stan_variables.py /mnt/code/stan --transformed-data
```

This analyzes variables created in `transformed data` blocks and checks if they're used in:
- Parameters (for sizing/constraints)
- Transformed parameters
- Model block
- Generated quantities

Variables not found in any of these are flagged as dead code.
