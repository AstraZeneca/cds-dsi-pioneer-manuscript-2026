# GitHub Copilot Instructions for Sclc Project

## Project Overview

This is the **sclc** project within the broader Pioneer CDS-DSI codebase. Sclc is a Bayesian hierarchical model for analyzing tumor dynamics and progression-free survival in oncology trials, specifically focused on the HISTORICAL and HISTORICAL-2 trials.

## Project-Specific Context

### What Sclc Uses

**Stan Models:**
- Primary model: `stan/ssls/sf-ssls-lfo.stan` and `stan/ssls/sf-ssm-log-space.stan`
- Modular architecture with components in `stan/ssls/modules/`
- Uses **n_causes = 1** (single combined competing risk for all non-target events)

**R Code Paths:**
- Data preparation: `r/sclc/prepare_analysis_data.R` (uses `prepare_tumor_stan_data()`)
- Initialization: `create_tumor_ssls_initializer()` in `r/initializers.R`
- Targets workflow: `targets/sclc_targets.R`
- Shared utilities: `r/priors.R`, `r/state_space.R`, `r/plot_functions.R`

**Key Design Decisions:**
- No patient-level covariate effects in other events model (only population and trial levels)
- Tumor covariates (`tumor_sum_covar`) are scaffolded but currently set to 0 (pending SSM integration)
- Uses modular Stan code with feature flags (e.g., `oe_enable_*` flags for other events)
- QR decomposition for numerical stability in covariate effects

### What Sclc Does NOT Use

**R Functions Not Used by Sclc:**
- `base_prepare_pfs_stan_data()` in `r/prepare_analysis_data.R` (used by endometrial-to-lung)
- `prepare_confirmed_resp_stan_data()` in `r/prepare_analysis_data.R` (used by other projects)
- `create_crcr_initializer()` and `create_crcr_pfs_initializer()` in `r/initializers.R`

**Stan Models Not Used:**
- Competing risks models in `stan/crcr/` (sclc uses SSLS, not CRCR)
- PFS-confirmed-response models in `stan/pfs-confirmed-response/`
- Breast-breast specific models

**Important:** When making changes to shared files like `r/prepare_analysis_data.R` or `r/priors.R`, ensure changes work for ALL projects (sclc, endometrial-to-lung, etc.), or make sclc-specific changes only in `r/sclc/` directory.

## Code Architecture

### Stan Module System

The Stan code uses a modular architecture with `#include` directives. Modules are organized by feature:

- **State Space Module** (`stan/ssls/modules/`):
  - `tr/` - Tumor regression (decrease) dynamics
  - `frac/` - Growth fraction dynamics
  - `init/` - Initial state modeling
  - `other_events/` - Competing risks for non-target events

- **Module Structure Pattern:**
  Each module has 7 standard files:
  1. `flags.stan` - Feature switches (e.g., `enable_trial_intercept_tr`)
  2. `data.stan` - Data declarations specific to module
  3. `hyperparams.stan` - Prior hyperparameters
  4. `transformed_data.stan` - Data preprocessing
  5. `parameters.stan` - Parameter declarations
  6. `transformed_parameters.stan` - Derived quantities
  7. `priors.stan` - Prior distributions

### Naming Conventions

**Module Prefixes:**
- Other events: `oe_*` for parameters/hyperparameters
- Tumor regression: `tr_*` for parameters/hyperparameters
- Growth fraction: `frac_*` for parameters/hyperparameters
- Initial state: `init_*` for parameters/hyperparameters

**Feature Flags:**
- Use verb "enable": `<module>_enable_<feature>`
- Examples: `oe_enable_trial_baseline_hazard`, `enable_trial_intercept_tr`

**Hierarchical Parameters:**
- Population level: `*_pop` (e.g., `tr_coef_qr_pop`)
- Trial level SD: `*_sd_trial_*` (e.g., `tr_sd_trial_intercept`)
- Trial level raw effects: `*_raw_trial_*` (e.g., `tr_raw_trial_intercept`)
- Patient level SD: `*_sd_patient_*` (e.g., `tr_sd_patient_intercept`)
- Patient level raw effects: `*_raw_patient_*` (e.g., `tr_raw_patient_intercept`)

**For detailed naming conventions, see:** `docs/multi_level_hierarchy_design.md`

### R Code Organization

- **Base functions** (`r/*.R`): Shared across all projects
- **Sclc-specific** (`r/sclc/*.R`): Only for sclc
- **Project workflows** (`targets/*.R`): Project-specific pipelines

## Documentation Maintenance Requirements

### When to Update Documentation

1. **Always check and update relevant documentation** when making code changes:
   - Module changes → Update module-specific docs
   - Naming changes → Update `docs/multi_level_hierarchy_design.md`
   - Architecture changes → Update this file and `README.md`
   - New features → Update relevant markdown files in `docs/`

2. **README files to maintain:**
   - `/mnt/code/README.md` - Main project overview
   - `/mnt/code/targets/README.md` - Targets workflow documentation
   - Module-specific README files (if they exist)

3. **Technical documentation files:**
   - `docs/ARCHITECTURE.md` - System architecture, naming conventions, and optimization techniques
   - `docs/OTHER_EVENTS_MODEL.md` - Other events model design and implementation
   - `docs/CHANGELOG.md` - Major changes and design decisions
   - `docs/CODEOWNERS` - Code ownership and review requirements
   - `sld_state_space_model.md` - State space model documentation

4. **Quarto documentation:**
   - Files in `quarto/` directory should be updated when analysis methods change
   - HTML outputs may need regeneration after changes

### Documentation Update Guidelines

- **Before implementing changes:** Read relevant markdown files to understand current design
- **After implementing changes:** Update affected documentation to reflect new behavior
- **When adding features:** Document design decisions in appropriate markdown files
- **When deprecating code:** Note deprecation in documentation and explain migration path
- **Keep examples current:** Update code examples in documentation when APIs change

### Markdown Formatting Guidelines

- **Lists after text:** Always include a blank line before a list when it follows text or a heading. This ensures proper rendering.
  ```markdown
  Some introductory text:
  
  - First item
  - Second item
  ```
  **Incorrect:**
  ```markdown
  Some introductory text:
  - First item
  - Second item
  ```

## Testing and Validation

### Before Committing Changes

1. **Compilation check:** Ensure Stan models compile using `stanc` directly
2. **Baseline regression:** Changes with all flags=FALSE should match previous behavior
3. **Incremental testing:** Test new features with flags enabled one at a time
4. **Documentation review:** Verify all relevant docs are updated

### Stan Model Validation

**Use `stanc` for syntax checking:**
```bash
~/.cmdstan/cmdstan-2.37.0/bin/stanc --include-paths=stan,stan/ssls stan/ssls/sf-ssm-log-space.stan
```

This is faster than full compilation and sufficient for checking syntax correctness. Only do full compilation with `cmdstanr::cmdstan_model()` when you need the executable.

### Key Files for Testing

- Test definitions: `tests/testthat/`
- Test runner: `tests/testthat.R`
- Example workflows: Files in `targets/` directory

## Common Operations

### Adding a New Module Parameter

1. Add feature flag in `modules/<module>/flags.stan`
2. Add hyperparameters in `modules/<module>/hyperparams.stan` with `<module>_` prefix
3. Add parameter declarations in `modules/<module>/parameters.stan`
4. Implement transformed parameter logic in `modules/<module>/transformed_parameters.stan`
5. Add priors in `modules/<module>/priors.stan`
6. Update `r/priors.R` with default hyperparameter values
7. Update relevant initializer in `r/initializers.R`
8. Update `targets/sclc_targets.R` with default flag value
9. **Update documentation** in `docs/ARCHITECTURE.md` and `docs/CHANGELOG.md`

### Adding Data to Stan Models

1. If module-specific: Add to `modules/<module>/data.stan`
2. If shared: Add to appropriate base data file (`stan/base_data.stan`, `stan/tumor/base_data.stan`, etc.)
3. Update R data preparation in `r/sclc/prepare_analysis_data.R`
4. Ensure `stan_data` list in targets file includes new data
5. **Update documentation** describing the new data and its purpose

### Modifying Shared Code

When editing files used by multiple projects (`r/priors.R`, `r/prepare_analysis_data.R`, etc.):

1. Check which projects use the function (search in `targets/` directory)
2. Ensure changes are backwards compatible OR make project-specific versions
3. For sclc-specific changes, use files in `r/sclc/` directory
4. Test with other projects if making shared changes
5. **Document breaking changes** in commit messages and relevant README files

## Version Control

- **Current branch:** Work is typically on `karim/non-target` or feature branches
- **Repository:** azu-oncology-rd/cds-dsi-pioneer-sclc-01-2025
- **Commit messages:** Should reference updated documentation files when applicable

## Stan Coding Guidelines

### Stan Best Practices

**Use Built-in Zero Constructors:**
- Use `zeros_vector()`, `zeros_row_vector()`, `zeros_int_array()`, `zeros_real_array()`, etc.
- Never manually create arrays/vectors of zeros

**Function Arguments:**
- Don't pass array/vector sizes as arguments - use `size()` or `num_elements()` inside the function
- For pos arrays, don't pass pos section sizes - infer them using appropriate pos functions

**Code Organization:**
- Stan code uses modular `#include` architecture
- Code fragments may be included in different sections than where they're defined
- Ignore linter warnings about code needing to be in specific sections (functions, transformed data, etc.)

### Stan Module Structure

Each module in `stan/ssls/modules/` follows this 7-file pattern:
1. `flags.stan` - Feature switches
2. `data.stan` - Module-specific data declarations
3. `hyperparams.stan` - Prior hyperparameters
4. `transformed_data.stan` - Data preprocessing
5. `parameters.stan` - Parameter declarations
6. `transformed_parameters.stan` - Derived quantities
7. `priors.stan` - Prior distributions

## R Coding Guidelines

### Code Style

- Follow tidyverse style guide
- Use modern pipe operator `|>` (not `%>%`)
- Prefer `purrr` and `dplyr` functions over base R loops
- Use `rlang::list_assign()` for list modifications

### Testing

- Use `testthat` framework for unit tests
- Tests located in `tests/testthat/`
- Run with `tests/testthat.R`

## References

### Configuration Files
- **Quarto config:** `_quarto.yml` for documentation generation settings
- **Targets config:** `_targets.yaml` and project-specific yaml files for workflow configuration
- **R environment:** `renv.lock` for package versions

### Key Documentation
- **Architecture & Design:** `docs/ARCHITECTURE.md`
- **Other Events Model:** `docs/OTHER_EVENTS_MODEL.md`
- **Change History:** `docs/CHANGELOG.md`
- **State space model:** `sld_state_space_model.md`
- **Code ownership:** `docs/CODEOWNERS`

---

**Last Updated:** November 14, 2025
**Maintainer:** This file should be updated whenever significant architectural or process changes occur.
