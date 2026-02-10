# Dead Functions Removal Report

**Date**: 2026-02-10

## Summary

Identified and removed 15 dead functions from Stan function libraries (16.1% of all defined functions). These functions were never called in the three active SSLS models.

## Analysis Method

Created automated analysis tool: `scripts/analyze_stan_functions.py`

**Key features:**
- Recursively resolves `#include` directives to build complete dependency graph
- Extracts function definitions from function libraries
- Traces function calls across all Stan files
- Excludes legacy directory from analysis
- Detects function references (e.g., passed to `map_rect`) in addition to direct calls

## Functions Removed

### stan/util.stan (2 functions, 13.3% of file)

1. `num_leq()` - Get number of values in array ≤ threshold
2. `assert_greater()` - Strict greater-than assertion (Note: `assert_greater_than_or_equal()` is still used)

### stan/pos.stan (3 functions, 27.3% of file)

1. `get_sub_vector()` - Extract vector sub-array for group
2. `get_sub_row_vector()` - Extract row vector sub-array for group
3. `get_last_int()` - Get last element from group

### stan/pfs_functions.stan (2 functions, 22.2% of file)

1. `calc_pfs_n()` - Calculate proportion surviving beyond time n
2. `calc_admin_brier_score()` - Administrative Brier Score for competing risks (~70 lines)

### stan/gp.stan (6 function overloads, 40.0% of file)

1-3. `ncp_gp_matern32()` - 3 overloaded variants (vector intercept, real intercept, row_vector)
4-5. `ncp_gp_matern52()` - 2 overloaded variants (vector, row_vector)
6. `multi_normal_lcdf()` - Conditional multivariate normal CDF

### stan/_sf_functions.stan (2 functions, 18.2% of file)

1. `sf_log_space_trajectory_lpdf()` - Log likelihood for complete trajectory (~30 lines)
2. `assert_matching_states()` - Debug assertion for state matching (~30 lines)

## Lines Removed

- **util.stan**: ~12 lines
- **pos.stan**: ~33 lines
- **pfs_functions.stan**: ~90 lines
- **gp.stan**: ~35 lines
- **_sf_functions.stan**: ~60 lines

**Total: ~230 lines of dead code removed**

## Functions Initially Flagged but Actually Used

During analysis, one function was initially misidentified as dead:

- `calc_patient_states_rect()` - Used as function reference in `map_rect()` call

This was caught by improving the pattern matcher to detect function references in addition to direct calls.

## Benefits

1. **Faster Compilation**: ~230 fewer lines to parse and check
2. **Reduced Complexity**: Clearer which functions are actually used
3. **Easier Maintenance**: Less code to understand and maintain
4. **Historical Preservation**: All removed functions archived in `stan/legacy/dead_functions.stan`

## Verification

All three models compile successfully after removal:

```bash
✅ stan/sf-ssm-log-space.stan
✅ stan/sf-ssls-lfo.stan
✅ stan/sf-ssls-lfo-endpoints.stan
```

Verified with:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/<model>.stan
```

## Files Modified

- `stan/util.stan`
- `stan/pos.stan`
- `stan/pfs_functions.stan`
- `stan/gp.stan`
- `stan/_sf_functions.stan`

## Files Created

- `stan/legacy/dead_functions.stan` - Archive of all removed functions
- `scripts/analyze_stan_functions.py` - Dead function analysis tool
- `stan/DEAD_FUNCTIONS_REPORT.md` (this file)

## Impact by File

| File | Functions Before | Functions After | Dead Functions | % Dead |
|------|------------------|-----------------|----------------|--------|
| util.stan | 15 | 13 | 2 | 13.3% |
| pos.stan | 11 | 8 | 3 | 27.3% |
| pfs_functions.stan | 9 | 7 | 2 | 22.2% |
| gp.stan | 15 | 9 | 6 | 40.0% |
| _sf_functions.stan | 11 | 9 | 2 | 18.2% |
| **Total** | **61** | **46** | **15** | **24.6%** |

(Note: Analysis covers function files only, not all 93 functions across entire codebase)

## Related Work

This cleanup follows two previous commits:

1. **5818a72**: Reorganized Stan directory structure
2. **683fdd6**: Removed ~180 lines of dead code from transformed data sections

Combined with this work, the Stan codebase cleanup has removed:
- ~180 lines of dead transformed data computations
- ~230 lines of dead function definitions
- **Total: ~410 lines of dead code removed**

## Future Work

The analysis tool can be extended to:
- Detect unused parameters within used functions
- Find functions only called from legacy code
- Identify redundant function overloads
