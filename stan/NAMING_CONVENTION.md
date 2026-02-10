# Naming Convention Compliance Check
**Date:** September 30, 2025  
**Branch:** karim/module-param-update  

## Status Summary

The modular refactoring has implemented the new naming convention **partially**. The module structure (tr, frac, init) is in place with proper gating flags, but several **legacy names are preserved for backward compatibility** with downstream R code.

## Convention Overview

From `tumor_total_fraction_init_inventory.md` Section 2.1:

**Pattern:** `<module>_<scale?>_<component>_<level>`

- **module**: `tr` (total rate), `frac` (fraction mix), `init` (initial proportion)
- **scale**: optional `logit`, `log`, etc.
- **component**: `loc`, `coef_qr`, `coef`, `linpred`, `sd`, `raw`, `effect`
- **level**: `pop`, `trial`, `patient`

## Compliance Status by Module

### ✅ Fully Compliant (New Infrastructure)

| Category | Examples | Status |
|----------|----------|--------|
| Flags | `enable_pop_cov_tr`, `enable_trial_intercept_frac`, etc. | ✅ All follow convention |
| Module params | `tr_coef_qr_pop`, `frac_sd_trial_intercept`, `init_raw_patient_slope` | ✅ All new params follow convention |
| Linpreds | `tr_linpred_pop`, `frac_linpred_patient_dev`, `init_linpred_pop` | ✅ Correctly named |
| Effects (new) | `frac_effect_trial_intercept`, `init_effect_patient_intercept` | ✅ Follow convention |

### ⚠️ Partially Compliant (Legacy Names Preserved)

| Current Name | Target Name (per inventory) | Location | Module | Reason |
|--------------|----------------------------|----------|--------|---------|
| `pop_log_total_rate` | `tr_loc_pop` | tr/parameters.stan | tr | Backward compat (documented) |
| `pop_decrease_frac_logit` | `frac_logit_loc_pop` | frac/parameters.stan | frac | Backward compat (commented) |
| `pop_decrease_prop_logis` | `init_logit_loc_pop` | init/parameters.stan | init | Backward compat (commented) |
| `patient_decrease_frac_logit` | `frac_logit_loc_patient` | frac/transformed_parameters.stan | frac | Used in legacy code |
| `patient_decrease_prop_logis` | `init_logit_loc_patient` | init/transformed_parameters.stan | init | Used in legacy code |
| `patient_log_decrease_frac` | `frac_log_decrease_patient` | frac/transformed_parameters.stan | frac | Used in calc_states |
| `patient_log_growth_frac` | `frac_log_growth_patient` | frac/transformed_parameters.stan | frac | Used in calc_states |
| `patient_log_decrease_prop` | `init_log_decrease_patient` | init/transformed_parameters.stan | init | Used in calc_states |
| `patient_log_growth_prop` | `init_log_growth_patient` | init/transformed_parameters.stan | init | Used in calc_states |
| `trial_log_total_rate_effect` | `tr_effect_trial_intercept` | tr/transformed_parameters.stan | tr | Used in legacy aggregator |
| `patient_log_total_rate_effect` | `tr_effect_patient_intercept` | tr/transformed_parameters.stan | tr | Used in legacy aggregator |
| `tr_loc_patient` | `patient_log_total_rate` (or vice versa?) | tr/transformed_parameters.stan | tr | Mixed usage |

### 📍 Key Dependencies

The legacy names are **actively used** in:
1. **`stan/ssls/legacy/sf-ssls-transformed_parameters.stan`** - References `patient_log_decrease_frac`, `patient_log_growth_frac`, `patient_log_decrease_prop`, `patient_log_growth_prop`
2. **`stan/tumor/` files** - Various references to old names
3. **R extraction code** - Likely expects legacy parameter names from Stan output
4. **`calc_states()` function** - Takes `patient_log_*` parameters as arguments

## Recommendations

### Option A: Add Alias Layer (Recommended for now)
In each module's `transformed_parameters.stan`, add aliases at the end:

```stan
// Backward compatibility aliases
vector[n_train_patients] patient_log_decrease_frac = frac_log_decrease_patient;
vector[n_train_patients] patient_log_growth_frac = frac_log_growth_patient;
```

**Pros:** No downstream changes needed yet  
**Cons:** Duplicates names in output, increases memory slightly

### Option B: Global Find-Replace
Update all downstream code (R scripts, Stan files, initializers) simultaneously.

**Pros:** Clean, follows convention fully  
**Cons:** Large coordinated PR, potential for missed references

### Option C: Phased Migration
1. Keep current state with documented legacy names
2. Update R extraction code to use new names
3. Remove aliases in subsequent PR

**Pros:** Incremental, testable  
**Cons:** Takes longer, interim state has two naming schemes

## Current State Assessment

The refactoring is **architecturally sound**:
- ✅ Module structure is clean and well-separated
- ✅ Flags follow convention precisely  
- ✅ New parameters follow naming convention
- ✅ No orphan parameters (everything gated properly)
- ⚠️ Legacy names preserved for continuity

## Next Steps

1. **Decision needed:** Choose migration strategy (A, B, or C above)
2. **Document dependencies:** Grep for all uses of legacy names in R code
3. **Create migration map:** List all downstream files that need updates
4. **Add tests:** Ensure parameter extraction works with either naming scheme

## Notes

- The **numerical stability improvements** from the rebase (log_inv_logit, log1m_inv_logit) **were successfully applied** to the module files
- The flags system is **fully implemented** and follows the convention exactly
- The internal module code is **clean and consistent**
- Only the **interface to legacy code** uses old names

---
**Conclusion:** The implementation follows the convention for all new infrastructure. Legacy names are preserved intentionally for compatibility. A migration plan is needed to complete the naming transition.
