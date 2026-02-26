# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Sclc is a Bayesian hierarchical modeling system for analyzing tumor dynamics and progression-free survival (PFS) in oncology trials. The codebase supports three projects:

1. **SCLC-01**: Longitudinal tumor analysis with state-space modeling (primary focus)
2. **Endometrial to LUNG**: Trial outcome predictions
3. **Breast-01 to -04**: Cross-validation predictions

### GitHub Repositories (azu-oncology-rd org)
- **This repo**: `cds-dsi-pioneer-core` — shared modeling core
- `cds-dsi-pioneer-sclc-01-2025` — SCLC-01 analysis
- `cds-dsi-pioneer-lung-2024` — LUNG analysis
- `cds-dsi-pioneer-pioneer-2026` — Pioneer analysis

## Build and Development Commands

### Initial Setup
```bash
# Restore R packages
Rscript -e 'renv::restore()'

# Install CmdStan
Rscript -e 'cmdstanr::install_cmdstan()'
```

### Running Analysis Pipelines
```bash
# Set project via TAR_PROJECT env var (see _targets.yaml for options)
export TAR_PROJECT=sclc
Rscript -e 'targets::tar_make()'
```

### Stan Model Syntax Check (fast)
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/sf-ssm-log-space.stan
```

### Running Tests
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Quarto Website

### Location and Purpose
The SCLC-01 analysis results are published as a Quarto website located in `quarto/website/`. The site is titled **"PIONEER: SCLC-01 Forecasting"** and includes:
- Documentation (onboarding tutorial, model specification)
- Analysis results (descriptive stats, model validation, outcome predictions)
- Interactive visualizations using R plots

### Building the Website Locally

**IMPORTANT**: Always render from the project root directory (`/mnt/code`), not from inside `quarto/website`. This is because the `_quarto.yml` sets `execute-dir: project`, which means R code needs access to the project-level renv and data.

```bash
# From project root
quarto render quarto/website
```

The rendered site is output to `quarto/website/_site/`. To preview:
```bash
quarto preview quarto/website
```

### Publishing to RStudio Connect

See `docs/PUBLISHING.md` for full publishing instructions (API keys, rsconnect setup, troubleshooting).

**Quick reference:**
```bash
# Render first (always from project root)
quarto render quarto/website

# Deploy via R
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/website", server = "az-connect", account = "kmjq089")'
```

### Key Configuration Files
- `_quarto.yml` - Site configuration, navigation, theme settings
- `az-theme.scss` - AstraZeneca color scheme (navy, gold, turquoise, etc.)
- `styles.css` - Custom CSS for hero section and layout
- `_freeze/` - Cache directory for executed R code (speeds up rebuilds)
- `images/` - PIONEER helmet logo and favicon

### Theme and Styling
The site uses **AZ corporate colors** defined in `az-theme.scss` (navy `#003865`, gold `#F0AB00`, turquoise `#68D2DF`, plum `#830051`, pink `#D0006F`, platinum `#9DB0AC`).

### Rebuilding After Changes
If you update plot functions in `r/plot_functions.R`, you may need to clear the freeze cache:
```bash
rm -rf quarto/website/_freeze/analysis/
```

Then re-render the affected pages.

### TikZ Diagrams
- Compile diagrams: `/home/ubuntu/.TinyTeX/bin/x86_64-linux/pdflatex -output-directory=quarto/website/images quarto/website/images/<name>.tex`
- Convert to SVG: `pdf2svg quarto/website/images/<name>.pdf quarto/website/images/<name>-tikz.svg`
- Clean artifacts: `rm -f quarto/website/images/<name>.{aux,log,pdf}`
- Use `/tikz-diagram` skill for guided workflow with prerequisite checks

### Documentation Sync
- Run `/sync-docs` after code changes to identify outdated documentation
- Updates doc-mappings.yaml when new source files are added to major features
- Architecture diagrams in `quarto/website/images/` may need regeneration when model structure changes

## Architecture

### Stan Module System

Stan code uses modular `#include` architecture in `stan/modules/`:
- **tr/** - Tumor regression (decrease) dynamics
- **frac/** - Growth fraction dynamics
- **init/** - Initial state modeling
- **multistate/** - Illness-death multistate model (replaced "other_events" module; 3 transitions: 0→1, 0→2, 1→2)
- **measurement/** - Observation model (measurement error)

Each module follows a 7-file pattern:
1. `flags.stan` - Feature switches
2. `data.stan` - Module-specific data
3. `hyperparams.stan` - Prior hyperparameters
4. `transformed_data.stan` - Data preprocessing
5. `parameters.stan` - Parameter declarations
6. `transformed_parameters.stan` - Derived quantities
7. `priors.stan` - Prior distributions

### Three-Level Hierarchical Parameters

Parameters follow a population → trial → patient hierarchy:
- Population level: `*_pop` (e.g., `tr_intercept_pop`)
- Trial level SD: `*_sd_trial_*` (e.g., `tr_sd_trial_intercept`)
- Trial level raw effects: `*_raw_trial_*` (non-centered parameterization)
- Patient level SD: `*_sd_patient_*`
- Patient level raw effects: `*_raw_patient_*`

### Module Prefixes
- `tr_*` - Tumor regression
- `frac_*` - Growth fraction
- `init_*` - Initial state
- `ms_*` - Multistate hazard

### Key Stan Models
- `stan/sf-ssm-log-space.stan` - Main state-space longitudinal survival model
- `stan/sf-ssls-lfo.stan` - Leave-future-out cross-validation variant

## Key Files

### R Code
- `r/priors.R` - Prior specifications (shared)
- `r/initializers.R` - Stan model initializers (shared)
- `r/util.R` - Utility functions including `select_draws(fit, matches(...))` for efficient parameter extraction from CmdStanR fits
- `r/sclc/prepare_analysis_data.R` - SCLC-specific data prep
- `targets/sclc_targets.R` - SCLC pipeline definition

### Configuration
- `_targets.yaml` - Workflow configurations for all projects
- `renv.lock` - Package versions
- `_quarto.yml` - Documentation settings

## Analysis Results Storage

Analysis results are stored in `/mnt/data/analysis-results/karim_naguib/sclc/<run_name>/`:
- **Targets store**: `/mnt/data/analysis-results/karim_naguib/sclc/<run_name>/_targets`
- **Fit CSVs**: `/mnt/data/analysis-results/karim_naguib/sclc/<run_name>/fit`

### Store Selection via TAR_BRANCH

Analysis results are organized by data cut-off (DCO) using the `TAR_BRANCH` environment variable:
- `export TAR_BRANCH=dco3` - January 26, 2026 DCO (current, includes pdl1_central)
- `export TAR_BRANCH=dco2` - August 2025 DCO
- `export TAR_BRANCH=dco1` - April 2025 DCO

Store path: `/mnt/data/analysis-results/<user>/sclc/<TAR_BRANCH>/_targets`

**Example:**
```bash
export TAR_BRANCH=dco3
Rscript -e 'targets::tar_make(sclc_patient_data_jan26)'
```

When working with stored targets directly (e.g., in standalone scripts), always specify the store path explicitly:
```r
tar_read(sclc_patient_data_jan26, store = "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets")
```

**Verifying fit versions**: `tar_outdated()` can be unreliable for `*_res_*` targets. Check timestamps directly:
```bash
ls -lt /mnt/data/analysis-results/.../fit/tumor_ssls_ctdna_jan26/*.csv | head -5
targets::tar_meta(tumor_ssls_km_*_jan26, store = "...")$time
```

### SCLC Data Files

Patient and visit data files use date-based suffixes indicating when they were processed:
- `cooked_patient_data_200226.csv` - Feb 26, 2026 (includes pdl1_central column)
- `cooked_patient_data_220126.csv` - Jan 22, 2026 (older, missing pdl1_central)
- `assessment_visit_data_200226.csv` - Feb 26, 2026

**Important columns:**
- `pdl1` - Site-reported PDL1 values (original baseline)
- `pdl1_central` - Centrally assessed PDL1 values (added Feb 2026, from PDL1CBL)
- `pdl1_hi` / `pdl1_high` - Binary flag (≥50% threshold, derived from pdl1)

**Note**: Central vs site PDL1 can show significant discordance. Always check which column is being used for stratification.

### PFS Endpoint Definitions

Model outputs three PFS variants (in `tumor_ssls_km_rvar_*` targets):
- `target_km_est` - RECIST PD only (death is censored) - **NOT comparable to clinical PFS**
- `ms_km_est` - Multistate hazard (0→1 transition)
- `km_est` - Combined (progression OR death, includes 0→2 events) - **use for comparison against observed PFS**

**Critical**: Observed PFS includes death without progression (0→2 events). HISTORICAL has 89 such events (14%), SCLC-01 has 4 (3%). Always use `km_est` when comparing model predictions to observed PFS KM.

### Observed KM Data Structure

`km_trial_pfs` and `km_trial_os` targets have `btype` column:
- `btype == "lb"` - Lower bound: events placed at `time + 1` (earliest possible)
- `btype == "ub"` - Upper bound: events placed at `time + interval_censored + 1` (latest possible)

Plots always use `btype == "ub"`. For PFS, `interval_censored` captures visit-gap uncertainty (up to 6 weeks). **For OS targets, always add `interval_censored = 0L` in the `mutate()` before calling `get_km_res()`** — deaths are observed exactly and the PFS interval_censored column must not be carried over, or it shifts death times forward and inflates the observed OS KM.

`dco-comparisons.qmd` requires `km_trial_pfs_*_ctdna_{apr25,aug25}` targets (5 variants × 2 DCOs). These are fast to build but are not auto-built — run `tar_make()` for them explicitly if missing from the store.

## Coding Guidelines

### General
- **No backward-compatibility aliases**: Do not create variable or function aliases for backward compatibility unless explicitly requested. When renaming, update all references directly instead of adding shims or aliases.

### Stan
- Use built-in zero constructors: `zeros_vector()`, `zeros_int_array()`
- Don't pass array/vector sizes as arguments; use `size()` internally
- Ignore linter warnings about code sections (modular `#include` architecture places code fragments across sections)
- Use non-centered parameterization (NCP) for hierarchical parameters
- Use `fatal_error()` instead of `reject()` for data validation errors in transformed data

### R
- Use modern pipe operator `|>` (not `%>%`)
- Follow tidyverse style guide
- Prefer `purrr` and `dplyr` over base R loops
- Use `testthat` for unit tests
- **NEVER hardcode subject IDs** (usubjid, patient_id, etc.) - always use dynamic selection or filtering

### Targets
- **NEVER use `tar_config_set(store = ...)`** - it changes global state and causes conflicts
- Always use explicit `store` argument: `tar_read(name, store = "path/_targets")`
- Same applies to all targets functions: `tar_meta()`, `tar_load()`, etc.
- **IMPORTANT**: `TAR_BRANCH` environment variable does NOT work with `tar_make()` - always use explicit `store="/path/_targets"` argument
- `_targets.yaml` sclc store uses `!expr` with `DOMINO_STARTING_USERNAME` and `TAR_BRANCH` — never hardcode username or branch in this file
- `tumor_ssls_draws_pop` selection: `time_invariant_coef_qr_*` and `time_varying_coef_*` params don't follow the `_pop` suffix — they need `matches("^(time_invariant|time_varying)_coef")` added to the `select_draws` call
- **NEVER inline complex code in targets** - extract to helper functions in `r/` directory
  - Target commands should be simple function calls, not multi-line code blocks
  - Example: Use `tar_target(name, my_function(arg))` not `tar_target(name, { ... complex code ... })`
  - Helper functions belong in appropriate `r/` subdirectories (e.g., `r/sclc/plot_functions.R`)

### Bash and Command Execution
- **NEVER pipe long-running commands to `head`, `tail`, or similar** when running in background - it prevents real-time output monitoring
- If you need to capture output while preserving streaming, use `tee` instead: `command | tee output.log`
- Background tasks automatically capture output to a file - don't truncate it with pipes
- Example:
  ```bash
  # ❌ DON'T: User can't see real-time progress
  Rscript -e 'targets::tar_make()' | head -100

  # ✓ DO: Full streaming output visible
  Rscript -e 'targets::tar_make()'

  # ✓ ALTERNATIVE: If you need to save output too
  Rscript -e 'targets::tar_make()' | tee build.log
  ```

### Data Pipeline Scripts
- **NEVER run data preparation scripts without explicit confirmation from the user**
- A hookify rule (`.claude/hookify.data-pipeline-confirmation.local.md`) blocks accidental execution of:
  - `r/data_preparation_pipeline/*/main_pipeline_data_*.R`
  - Scripts containing `wrangle_data`, `prepare.*data`, or `pipeline.*` patterns
- **Why this matters**: Data pipeline scripts regenerate CSV files and can take hours to run
- **Before running**: Check if processed CSV files already exist at the expected location
- **Better alternative**: Use `targets::tar_make()` with specific target names to rebuild only what changed
- **Pattern**: Data pipelines should only run when raw data is updated, not for routine analysis

### Quarto and Documentation
- **Always use "SCLC-01"** when referring to the trial in user-facing text (documentation, plots, presentations)
- Use lowercase "sclc" only for code identifiers (variable names, trial codes, file paths)
- Example: Write "SCLC-01 trial" in figure captions, but `filter(trial == "sclc")` in R code
- **Use automatic section numbering**: Set `number-sections: true` in frontmatter, don't use manual numbers (1.1, 2.3) in headings
- **Cross-references**: Use section IDs `{#sec-name}` and reference with `@sec-name`, never hardcode "Section X.Y.Z"
- **Model specification is the blueprint**: `quarto/website/documentation/model-specification.qmd` is the authoritative specification for everything in the Stan model. Code and documentation must always match:
  - When changing Stan code, update the model specification to reflect the change
  - When the specification defines behavior (e.g., index conventions, endpoint formulas, routing logic), the code must not violate those definitions without updating the spec first
  - If a proposed code change contradicts the specification, flag the discrepancy before implementing
  - Treat the specification as a contract: it documents what the model *should* do, not just what it *happens* to do

### Git Worktree Workflow
When working with multiple git worktrees, follow this pattern to avoid duplicate commits:

1. **Make changes in your current worktree** - Don't `cd` to other worktrees to make the same changes
2. **Commit locally** - Commit your work in the worktree where you made the changes
3. **Merge from the target worktree** - Switch to the other worktree and merge or rebase

**Example:**
```bash
# In worktree A (/mnt/code/worktrees/code-feature-x): make changes and commit
git add .claude/settings.json
git commit -m "Update plugin configuration"

# In worktree B (/mnt/code on main branch): merge the changes
cd /mnt/code
git merge feature-x
git push origin main
```

**Why this pattern?**
- Avoids duplicate commits (same change, different SHAs)
- Keeps cleaner git history
- More efficient than manually replicating changes across worktrees
- Leverages git's merge/rebase capabilities

**Common mistake:**
```bash
# ❌ DON'T DO THIS:
# Making the same change in multiple worktrees separately
cd /mnt/code && edit file && git commit
cd /mnt/code/worktrees/code-feature-x && edit file && git commit  # Duplicate!

# ✓ DO THIS INSTEAD:
# Make change once, then merge
cd /mnt/code/worktrees/code-feature-x && edit file && git commit
cd /mnt/code && git merge feature-x
```

### Adding Module Parameters
1. Add feature flag in `modules/<module>/flags.stan`
2. Add hyperparameters in `modules/<module>/hyperparams.stan`
3. Add parameters in `modules/<module>/parameters.stan`
4. Implement transforms in `modules/<module>/transformed_parameters.stan`
5. Add priors in `modules/<module>/priors.stan`
6. Update `r/priors.R` with defaults
7. Update initializer in `r/initializers.R`
8. Update `targets/sclc_targets.R` with flag value

## Pull Request Checklist

Before creating or submitting a PR, verify all items below. These are **mandatory requirements**:

### Critical Rules (Will Block Merge)

- [ ] **NO hardcoded subject IDs anywhere** - No patient IDs, usubjid, subject identifiers, or any specific subject codes in:
  - Production code
  - Debug statements
  - Example code
  - Comments or documentation
  - Test fixtures (use synthetic IDs instead)
  - Use dynamic selection: `slice_sample(n=1)`, `filter()`, or parameterize the ID

- [ ] **NO install.packages() calls** - Dependencies must be managed via `renv`:
  - Remove all `install.packages()` and `if (!require(pkg)) install.packages(pkg)`
  - Add new packages to `renv.lock` using `renv::snapshot()`

- [ ] **Use native pipe `|>` not `%>%`** - Modern R pipe operator required throughout

- [ ] **NO syntax errors** - Code must parse without errors:
  - Run basic syntax check: `Rscript -e 'source("your_file.R")'`
  - For Stan: `stanc --include-paths=stan,stan/ssls your_model.stan`

### Code Quality Requirements

- [ ] **Follow tidyverse style guide** for R code

- [ ] **NO global assignments (`<<-`)** in functions - Return values instead

- [ ] **Proper function defaults** - Don't reference undefined variables in parameter defaults

- [ ] **Documentation matches implementation** - Roxygen docs must match actual return values and parameters

- [ ] **NO tar_config_set()** - Use explicit `store=` argument in all targets functions

### Recommended Pre-PR Workflow

```bash
# 1. Before committing - Quick validation
Rscript -e 'styler::style_file("r/your_modified_file.R")'  # Auto-format
grep -r "install.packages" r/  # Check for install calls
grep -r "%>%" r/  # Check for old pipe operator

# 2. Search for hardcoded subject IDs (common patterns)
grep -rE "(E[0-9]{7,}|[A-Z][0-9]{4}[0-9]{3}[0-9]{3})" r/ --include="*.R"

# 3. Test your changes
Rscript -e 'source("r/your_file.R")'  # Basic syntax check
Rscript -e 'testthat::test_dir("tests/testthat")'  # Run tests

# 4. Before creating PR - Automated review
# In Claude Code CLI:
/pr-review-toolkit:review-pr code
```

### Using Automated PR Review

The repository includes automated code review agents. Before creating your PR:

```bash
# In Claude Code
/pr-review-toolkit:review-pr code  # Check code quality and compliance
/pr-review-toolkit:review-pr tests  # Verify test coverage
/pr-review-toolkit:review-pr all    # Comprehensive review (recommended)
```

This will catch most issues automatically before reviewers see your PR.

## GitHub Project Management

Issues across all PIONEER repos are tracked in the **PIONEER 2026** GitHub Project (project number 56, owner `azu-oncology-rd`). See `docs/GITHUB_PROJECT.md` for full reference (project/field IDs, `gh` commands, GraphQL queries).

## Pioneer Claude Marketplace

- Plugin registry: `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/.claude-plugin/marketplace.json`
- **Adding a plugin**: create files in `plugins/<name>/` AND register it in `marketplace.json` — without the registry entry it won't appear in `/plugin`
- **Version bumps are required for deployment**: bump version in BOTH `plugins/<name>/.claude-plugin/plugin.json` AND `marketplace.json` whenever adding/changing tools — the plugin manager won't reinstall otherwise
- MCP tool naming convention: `mcp__plugin_<plugin-name>_<server-name>__<tool-name>` (e.g. `mcp__plugin_domino-toolkit_domino__start_job`)
- Marketplace repo: `https://github.com/azu-oncology-rd/cds-dsi-pioneer-claude-marketplace`

## Domino Environment

- Key env vars auto-set by Domino: `DOMINO_USER_API_KEY`, `DOMINO_USER_HOST`, `DOMINO_PROJECT_ID`, `DOMINO_PROJECT_NAME`
- Jobs API: list/get via `GET /api/jobs/beta/jobs`, logs via `GET /api/jobs/beta/jobs/{id}/logs`, start via `POST /v4/jobs/start`, stop via `POST /v4/jobs/stop`
- `stop_job` requires both `projectId` AND `jobId` in the request body

## Documentation

Key technical docs in `docs/`:
- `ARCHITECTURE.md` - System architecture and optimization details
- `ar1_process_noise.md` - AR(1) time-varying process noise implementation
- `OTHER_EVENTS_MODEL.md` - Other events model design
- `CHANGELOG.md` - Major changes and design decisions
