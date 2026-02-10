# Stan Directory Reorganization - COMPLETED

**Date**: 2026-02-10

## Summary

Successfully reorganized the Stan directory structure by:
1. Merging duplicate `base_data.stan` files
2. Moving all SSLS code from `stan/ssls/` to `stan/`
3. Archiving unused legacy code
4. Updating all include paths

## Changes Made

### 1. Merged base_data.stan ✅
- **Before**: Two separate files
  - `/stan/base_data.stan` - General structure
  - `/stan/tumor/base_data.stan` - Tumor-specific data
- **After**: Single merged file at `/stan/base_data.stan`
- **Archived**: Originals moved to `/stan/legacy/`

### 2. Moved SSLS Models Up ✅
All files moved from `/stan/ssls/` → `/stan/`:
- **Main Models**:
  - `sf-ssm-log-space.stan` - Primary state-space model
  - `sf-ssls-lfo.stan` - Leave-future-out CV model
  - `sf-ssls-lfo-endpoints.stan` - LFO with endpoints

- **Supporting Files**:
  - `_endpoints_generated_quantities.stan`
  - `_lfo_endpoints_generated_quantities.stan`
  - `_lfo_transformed_data.stan`
  - `_mature_cutoffs_transformed_data.stan`
  - `_sf_accuracy_generated_quantities.stan`
  - `_sf-checks.stan`
  - `_sf_functions.stan`
  - `_sf_outcomes_info.stan`
  - `_sf-ssls-lfo-data.stan`
  - `_sf_transformed_data.stan`
  - `_sf_transformed_parameters.stan`

- **Directories**:
  - `modules/` - Modular parameter system (frac, init, measurement, other_events, tr)

- **Documentation**:
  - `MODULE_DESIGN.md`
  - `NAMING_CONVENTION.md`
  - `README.md`

### 3. Archived Legacy Code ✅
Moved to `/stan/legacy/`:

**Unused Directories**:
- `baseline_hazard/` - Baseline hazard models
- `crcr/` - Competing risks models
- `breast/` - Breast trials
- `pfs-confirmed-response/` - PFS-CR models
- `recruit/` - Recruitment models
- `ssls_legacy/` - Old SSLS code
- `tumor/` - Unused tumor model files

**Unused Files**:
- `base_data_original.stan` - Original pre-merge
- `state_space.stan`
- `pfs_transformed_data.stan`
- `test_pos_all.stan`

### 4. Updated Include Paths ✅
Changed all `#include` statements from:
```stan
#include "../util.stan"
#include "../tumor/base_data.stan"
```

To:
```stan
#include "util.stan"
#include "base_data.stan"
```

### 5. Cleaned Up Directories ✅
- Removed empty `/stan/ssls/` directory
- Kept only `tumor_transformed_data.stan` in `/stan/tumor/`

## Final Structure

```
/mnt/code/stan/
├── base_data.stan                         # MERGED
├── base_transformed_data.stan
├── util.stan
├── pos.stan
├── gp.stan
├── pfs_functions.stan
├── lfo.stan
├── recist.stanfunctions
│
├── sf-ssm-log-space.stan                  # Main SSLS model
├── sf-ssls-lfo.stan                       # LFO model
├── sf-ssls-lfo-endpoints.stan             # LFO endpoints
│
├── _endpoints_generated_quantities.stan
├── _lfo_endpoints_generated_quantities.stan
├── _lfo_transformed_data.stan
├── _mature_cutoffs_transformed_data.stan
├── _sf_accuracy_generated_quantities.stan
├── _sf-checks.stan
├── _sf_functions.stan
├── _sf_outcomes_info.stan
├── _sf-ssls-lfo-data.stan
├── _sf_transformed_data.stan
├── _sf_transformed_parameters.stan
│
├── MODULE_DESIGN.md
├── NAMING_CONVENTION.md
├── README.md
│
├── modules/
│   ├── frac/
│   ├── init/
│   ├── measurement/
│   ├── other_events/
│   └── tr/
│
├── tumor/
│   └── tumor_transformed_data.stan
│
└── legacy/                                # ARCHIVED
    ├── base_data_original.stan
    ├── baseline_hazard/
    ├── crcr/
    ├── breast/
    ├── pfs-confirmed-response/
    ├── recruit/
    ├── ssls_legacy/
    ├── tumor/
    ├── state_space.stan
    ├── pfs_transformed_data.stan
    └── test_pos_all.stan
```

## Verification

All three main Stan models compile successfully:
```bash
cd /mnt/code
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/sf-ssls-lfo.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/sf-ssls-lfo-endpoints.stan
```

All pass without errors ✅

## Benefits

1. **Simplified structure**: SSLS code no longer nested in subdirectory
2. **Cleaner includes**: No more `../` in include paths
3. **Reduced duplication**: Single `base_data.stan` instead of two
4. **Clear legacy separation**: Unused code archived, not deleted
5. **Maintained functionality**: All models compile and work identically

## Next Steps

1. Update any R scripts that reference Stan file paths
2. Update CLAUDE.md if it references old paths
3. Consider cleaning up utility files to remove unused functions (future optimization)
4. Git commit these changes with descriptive message
