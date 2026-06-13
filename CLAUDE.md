# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Bayesian hierarchical modeling system for analyzing cancer burden dynamics, progression-free survival (PFS), and overall survival (OS) in oncology trials. The codebase supports projects:

1. **SCLC-01**: Longitudinal tumor analysis with state-space modeling (primary focus)
2. **Endometrial to LUNG**: Trial outcome predictions
3. **Breast-01 to -04**: Cross-validation predictions
4. **Pioneer** Using PSA as the disease burden

### GitHub Repositories (azu-oncology-rd org)
- `cds-dsi-pioneer-core` — shared modeling core
- `cds-dsi-pioneer-sclc-01-2025` — SCLC-01 analysis
- `cds-dsi-pioneer-lung-2024` — LUNG analysis
- `cds-dsi-pioneer-pioneer-2026` — Pioneer analysis

Look at the TAR_PROJECT env var in .Renviron to see which project we are in.

## Build and Development Commands

### Initial Setup
```bash
# Restore R packages
Rscript --no-init-file -e 'renv::restore()'

# Install CmdStan
Rscript -e 'cmdstanr::install_cmdstan()'
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

Stan code uses modular `#include` architecture in `stan/modules/<name>/`. Each module contains a subset of standard files (`flags.stan`, `data.stan`, `hyperparams.stan`, `parameters.stan`, `priors.stan`, `transformed_data.stan`, `transformed_parameters.stan`, `generated_quantities.stan`, `likelihood.stan`) — which files exist varies by module. Browse `stan/modules/` to see what's available. Shared top-level includes live in `stan/_*.stan` and `stan/*.stanfunctions`. See `stan/MODULE_DESIGN.md` and `stan/NAMING_CONVENTION.md` for the design rationale.

### Hierarchical Parameters

The hierarchy is level-agnostic: `n_levels` intermediate levels between population and patient, configured via `n_groups_per_level` and `patient_level_groups` (see `stan/_hierarchy_data.stan`). Each level's intercept mode (`none`/`fe`/`re`/`re_gp`) is set per module via `enable_level_intercept_*` flags.

Naming convention:
- Population: `*_pop` (e.g., `tr_loc_pop`)
- Per-level SD: `*_sd_level_*` (e.g., `tr_sd_level_intercept`)
- Per-level raw NCP effects: `*_raw_level_*`
- Patient: `*_raw_patient_*`, `*_sd_patient_*`

### Module Prefixes
- `tr_*` - Tumor regression
- `frac_*` - Growth fraction
- `init_*` - Initial state
- `ms_*` - Multistate hazard

## Coding Style

- **Never use `!!!` (rlang splice) inside a `tar_target()` command.** It forces
  evaluation at manifest/parse time — *before* targets resolves dependencies — so
  any other target referenced inside the splice fails with `object '<name>' not
  found` during job setup (before `tar_make` runs). Use a function that takes the
  list as a plain argument instead, so the dependency stays a normal unquoted
  target reference: e.g. `x |> modifyList(other_target)` to override-merge, NOT
  `x |> list_modify(!!!other_target)`. (`c(x, other_target)` also works but
  *appends* — duplicate keys, not overrides.)
- **Always use tidyverse** — `dplyr`, `purrr`, `tidyr`, `ggplot2`, `stringr`, etc.
- **Use `|>` (native pipe)**, never `%>%`
- **Use MCP targets tools** (`mcp__plugin_targets-toolkit_targets__*`) for all targets operations — never `Rscript -e 'targets::...'` via Bash
- **Use `/start-job`** to launch Domino jobs — never raw curl

## Key Files

### R Code
- `r/util.R` - Utility functions including `select_draws(fit, matches(...))` for efficient parameter extraction from CmdStanR fits

Per-project code lives in `r/<project>/`; pipelines in `targets/<project>_targets.R`.

### Configuration
- `renv.lock` - Package versions
- `sclc_targets.sh` / `pioneer_targets.sh` / `publication_targets.sh` - Shell scripts for launching Domino jobs

**Always use the project-matching script**: `publication_targets.sh` for `TAR_PROJECT=publication`, `sclc_targets.sh` for sclc, `pioneer_targets.sh` for pioneer. Never use `sclc_targets.sh` for the publication pipeline.

## Analysis Results Storage

### Project and Store Context (ALWAYS CHECK FIRST)

**NEVER assume a project or store path without confirming.** This codebase supports multiple projects — do not default to sclc.

The active project is set via `TAR_PROJECT`. Store paths follow: `/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/<TAR_PROJECT>/<TAR_RUN>/_targets`

`TAR_RUN` selects the analysis run (named by data cut-off or feature branch). List available runs:
```bash
ls /mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/<project>/
```

Project-specific details (TAR_RUN values, data files, endpoint definitions) are in `.claude/rules/sclc.md` and `.claude/rules/pioneer.md`.

## GitHub Project Management

Issues across all PIONEER repos are tracked in the **PIONEER** GitHub Project (project number 56, owner `azu-oncology-rd`). See `docs/GITHUB_PROJECT.md` for full reference (project/field IDs, `gh` commands, GraphQL queries).

## Pioneer Claude Marketplace

- **All edits go in the marketplace repo clone at `/home/ubuntu/claude-marketplace/`**, never in this analysis repo. Commit and push from there.
- **Adding a plugin**: create files in `plugins/<name>/` AND register it in `.claude-plugin/marketplace.json` — without the registry entry it won't appear
- **Version bumps are required for deployment**: bump version in BOTH `plugins/<name>/.claude-plugin/plugin.json` AND `.claude-plugin/marketplace.json` (all paths relative to `/home/ubuntu/claude-marketplace/`) — the plugin manager won't reinstall without a bump in both
- MCP tool naming convention: `mcp__plugin_<plugin-name>_<server-name>__<tool-name>` (e.g. `mcp__plugin_domino-toolkit_domino__start_job`)

## Domino Environment

- Key env vars auto-set by Domino: `DOMINO_USER_API_KEY`, `DOMINO_USER_HOST`, `DOMINO_PROJECT_ID`, `DOMINO_PROJECT_NAME`
- Jobs API: list/get via `GET /api/jobs/beta/jobs`, logs via `GET /api/jobs/beta/jobs/{id}/logs`, start via `POST /v4/jobs/start`, stop via `POST /v4/jobs/stop`
- `stop_job` requires both `projectId` AND `jobId` in the request body
- **`get_job_logs` MUST always use `tail=N`** — never call without it; full logs are 700+ lines and will flood the context window. Use `tail=30` for status checks, `tail=50` for error diagnosis. See `pioneer-toolkit:read-job-logs` skill for the full pattern.

## Documentation

Key technical docs in `docs/`:
- `CHANGELOG.md` - Major changes and design decisions
- `propensity-borrowing-design.md` - Propensity-weighted borrowing design rationale
- `multi_level_hierarchy_design.md` - Three-level hierarchy design notes
- `pfs-os-relationship-spec.md` / `pfs-os-surrogacy-discussion.md` - PFS-OS linkage
- `stan_state_space_optimization.md` - Stan SSM performance notes
