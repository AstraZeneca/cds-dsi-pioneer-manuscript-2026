# .stanfunctions Naming Convention

**Date**: 2026-02-10

## Summary

Standardized naming convention for Stan files that contain only function definitions by using the `.stanfunctions` extension and removing redundant "functions" from filenames.

## Motivation

1. **Clarity**: `.stanfunctions` extension makes it immediately clear the file contains only functions (no data, parameters, model blocks)
2. **Consistency**: Aligns with existing `recist.stanfunctions` convention
3. **Cleaner names**: Removes redundant "functions" from filenames (e.g., `pfs_functions.stan` → `pfs.stanfunctions`)

## Changes

| Before | After | Notes |
|--------|-------|-------|
| `util.stan` | `util.stanfunctions` | General utility functions |
| `pos.stan` | `pos.stanfunctions` | Position array functions |
| `gp.stan` | `gp.stanfunctions` | Gaussian process functions |
| `pfs_functions.stan` | `pfs.stanfunctions` | PFS/survival functions (removed "functions") |
| `lfo.stan` | `lfo.stanfunctions` | Leave-future-out functions |
| `sf_state_space.stan` | `sf_state_space.stanfunctions` | State-space model functions |
| `modules/tumor/functions.stan` | `modules/tumor/tumor.stanfunctions` | Tumor-specific functions |
| `recist.stanfunctions` | (no change) | Already followed convention |

## Convention Rules

### Use `.stanfunctions` when:
- File contains ONLY function definitions
- No `data`, `parameters`, `transformed parameters`, `model`, or `generated quantities` blocks

### Use `.stan` when:
- File contains data/parameter/model blocks
- File is a complete model
- File contains mixed content (e.g., `base_data.stan`, `_sf_transformed_data.stan`)

### Naming pattern:
- `<descriptive_name>.stanfunctions` (e.g., `util.stanfunctions`, `gp.stanfunctions`)
- For module function files: `modules/<module>/<module>.stanfunctions` (e.g., `modules/tumor/tumor.stanfunctions`)
- Remove redundant "functions" from name (e.g., `pfs_functions` → `pfs`)

## Files Updated

### Model files (includes updated):
- `stan/sf-ssm-log-space.stan`
- `stan/sf-ssls-lfo.stan`
- `stan/sf-ssls-lfo-endpoints.stan`

### Example include section after changes:
```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "sf_state_space.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
  #include "recist.stanfunctions"
}
```

## Benefits

1. **Immediate clarity**: File extension tells you content type
2. **Consistent convention**: All function-only files follow same pattern
3. **Cleaner organization**: Easy to identify function libraries vs other file types
4. **Future-proof**: When adding PSA module, follow same pattern: `modules/psa/psa.stanfunctions`

## Verification

All three models compile successfully:
```bash
✅ stan/sf-ssm-log-space.stan
✅ stan/sf-ssls-lfo.stan
✅ stan/sf-ssls-lfo-endpoints.stan
```

## Future Pattern

When adding new function-only files:
- Use `.stanfunctions` extension
- Name descriptively without "functions" suffix
- For modules: `modules/<module>/<module>.stanfunctions`

Example for PSA module:
```
modules/psa/
├── data.stan                    # Data inputs (not functions)
├── transformed_data.stan        # Data processing (not functions)
└── psa.stanfunctions           # PSA-specific functions
```
