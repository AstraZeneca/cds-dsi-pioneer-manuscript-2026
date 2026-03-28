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
# Tumor models:
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan

# PSA models:
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan
```

### Running Tests
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Architecture

### Stan Module System

Stan code uses modular `#include` architecture in `stan/modules/`:
- **tr/** - Tumor regression (decrease) dynamics
- **frac/** - Growth fraction dynamics
- **init/** - Initial state modeling
- **multistate/** - Illness-death multistate model (replaced "other_events" module; 3 transitions: 0→1, 0→2, 1→2)
- **measurement/** - Observation model (measurement error)
- **propensity/** - Propensity-weighted borrowing from RWD (pioneer only; 6-file pattern — no `data.stan`)

Each module follows a 7-file pattern (propensity uses 6 — no `data.stan`):
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
- `stan/tumor/sf-ssm-log-space.stan` - Main state-space longitudinal survival model
- `stan/tumor/sf-ssls-lfo.stan` - Leave-future-out cross-validation variant
- `stan/ms-standalone.stan` - Standalone multistate model (requires initializer; always pass `likelihood_weight` to `multistate_lpmf`)
- `stan/psa/pioneer.stan` - PSA/RWD borrowing model with propensity weighting

## Key Files

### R Code
- `r/priors.R` - Prior specifications (shared)
- `r/initializers.R` - Stan model initializers (shared)
- `r/util.R` - Utility functions including `select_draws(fit, matches(...))` for efficient parameter extraction from CmdStanR fits
- `r/sclc/prepare_analysis_data.R` - SCLC-specific data prep
- `r/sclc/multistate.R` - Multistate classification and field derivation
- `r/pioneer/prepare_analysis_data.R` - Pioneer data prep (propensity, visit-gated covariate)
- `targets/sclc_targets.R` - SCLC pipeline definition
- `targets/pioneer_targets.R` - Pioneer pipeline definition

### Configuration
- `_targets.yaml` - Workflow configurations for all projects
- `renv.lock` - Package versions
- `sclc_targets.sh` / `pioneer_targets.sh` - Shell scripts for launching Domino jobs

## Pioneer Model Variants

Defined in `tar_map()` inside `targets/pioneer_targets.R`:

| Model | `hist` | `propensity` | `visit_gated_01` | Description |
|-------|--------|--------------|-----------------|-------------|
| `combined` | TRUE | FALSE | FALSE | Trial + RWD, trial-level hierarchy |
| `no_hist` | FALSE | FALSE | FALSE | Trial only |
| `propensity` | TRUE | TRUE | TRUE | Density-ratio weighted RWD borrowing |
| `propensity_ungated` | TRUE | TRUE | FALSE | Propensity without visit-gated PSA covariate |
| `visit_gated_no_trial` | TRUE | FALSE | TRUE | Visit-gated covariate, no propensity |

Target naming follows nested `tar_map`: outer = model, inner = type (`posterior`/`prior`). E.g. `pioneer_fit_res_posterior_propensity`.

## Analysis Results Storage

### Project and Store Context (ALWAYS CHECK FIRST)

**NEVER assume a project or store path without confirming.** This codebase supports multiple projects — do not default to sclc.

The active project is set via `TAR_PROJECT` (see `_targets.yaml`). Store paths by project:
- **sclc**: `/mnt/data/analysis-results/<user>/sclc/<TAR_BRANCH>/_targets`
- **pioneer**: `/mnt/data/analysis-results/<user>/pioneer/<TAR_BRANCH>/_targets`

When the project or store name is ambiguous, list available stores first:
```bash
ls /mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/<project>/
```

### Store Selection via TAR_BRANCH

`TAR_BRANCH` selects the analysis run (named by data cut-off or feature branch). For sclc:
- `export TAR_BRANCH=dco3` - January 26, 2026 DCO (current, includes pdl1_central)
- `export TAR_BRANCH=dco2` - August 2025 DCO
- `export TAR_BRANCH=dco1` - April 2025 DCO

For pioneer, common branches include: `main`, `rwd-filtered-1`, `laplace`, `multistate`, etc.

Full store path pattern: `/mnt/data/analysis-results/<user>/<TAR_PROJECT>/<TAR_BRANCH>/_targets`

### Sclc Fit Output

Sclc results are stored in `/mnt/data/analysis-results/<username>/sclc/<run_name>/`:
- **Targets store**: `/mnt/data/analysis-results/<username>/sclc/<run_name>/_targets`
- **Fit CSVs**: `/mnt/data/analysis-results/<username>/sclc/<run_name>/fit`

**Example:**
```bash
export TAR_BRANCH=dco3
Rscript -e 'targets::tar_make(sclc_patient_data_jan26)'
```

When working with stored targets directly (e.g., in standalone scripts), always specify the store path explicitly:
```r
tar_read(sclc_patient_data_jan26, store = file.path("/mnt/data/analysis-results", Sys.getenv("DOMINO_STARTING_USERNAME"), "sclc/dco3/_targets"))
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

## GitHub Project Management

Issues across all PIONEER repos are tracked in the **PIONEER 2026** GitHub Project (project number 56, owner `azu-oncology-rd`). See `docs/GITHUB_PROJECT.md` for full reference (project/field IDs, `gh` commands, GraphQL queries).

### Team Members (azu-oncology-rd)

| Username | Name |
|---|---|
| `kmjq089_azu` | Naguib, Karim |
| `kfvz858_azu` | Berché, Roger |
| `kjmr060_azu` | Metcalfe, Paul |
| `kzht939_azu` | Bevan, Antonia |
| `kqdr852_azu` | Li, Lu |

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
- **`get_job_logs` MUST always use `tail=N`** — never call without it; full logs are 700+ lines and will flood the context window. Use `tail=30` for status checks, `tail=50` for error diagnosis. See `pioneer-toolkit:read-job-logs` skill for the full pattern.

## Documentation

Key technical docs in `docs/`:
- `ARCHITECTURE.md` - System architecture and optimization details
- `ar1_process_noise.md` - AR(1) time-varying process noise implementation
- `OTHER_EVENTS_MODEL.md` - Other events model design
- `CHANGELOG.md` - Major changes and design decisions
- `propensity-borrowing-design.md` - Propensity-weighted borrowing design rationale
