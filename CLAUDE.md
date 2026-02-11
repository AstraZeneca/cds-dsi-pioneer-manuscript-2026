# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Sclc is a Bayesian hierarchical modeling system for analyzing tumor dynamics and progression-free survival (PFS) in oncology trials. The codebase supports three projects:

1. **SCLC-01**: Longitudinal tumor analysis with state-space modeling (primary focus)
2. **Endometrial to LUNG**: Trial outcome predictions
3. **Breast-01 to -04**: Cross-validation predictions

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
~/.cmdstan/cmdstan-2.37.0/bin/stanc --include-paths=stan stan/sf-ssm-log-space.stan
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

**Server URL**: `https://rstudio-connect.seml.scp.astrazeneca.net/connect/`
**Account**: `kmjq089`

#### Step-by-Step Publishing Instructions

**IMPORTANT**: Always render the site locally first to ensure everything works. Run from the project root:
```bash
# From /mnt/code (project root)
quarto render quarto/website
```

#### Method 1: Using R rsconnect (Recommended)

This is the most reliable method. Follow these steps:

**Step 1: Get your API key**
1. Go to `https://rstudio-connect.seml.scp.astrazeneca.net/connect/` in your browser
2. Sign in with SSO
3. Click your name (top right) → "API Keys"
4. Click "New API Key" and copy it

**Step 2: Configure credentials (one-time setup)**
```r
library(rsconnect)

# Add the server
rsconnect::addConnectServer(
  url = "https://rstudio-connect.seml.scp.astrazeneca.net",
  name = "az-connect"
)

# Add your API key (paste your actual key here)
rsconnect::connectApiUser(
  account = "kmjq089",
  server = "az-connect",
  apiKey = "YOUR_API_KEY_HERE"
)
```

This stores credentials in `~/.rsconnect/` (outside the git repo).

**Step 3: Deploy the website**

**IMPORTANT**: Always run from the project root directory (where renv is configured), not from inside the website directory.

```r
# From /mnt/code directory
library(rsconnect)
rsconnect::deploySite(
  siteDir = "quarto/website",
  server = "az-connect",
  account = "kmjq089"
)
```

**Step 4: Find your published site**
After deployment completes, look for the URL in the output or go to:
`https://rstudio-connect.seml.scp.astrazeneca.net/connect/#/content/listing?q=owner:kmjq089`

#### Method 2: Using Quarto CLI

This requires browser authentication (SSO) and may not work in remote CLI environments.

```bash
# From /mnt/code (project root)
quarto publish connect quarto/website --server https://rstudio-connect.seml.scp.astrazeneca.net/connect/
```

This will open your browser for SSO authentication.

#### Updating an Existing Deployment

Once you've published the first time, subsequent deployments are simple:

```r
# From /mnt/code directory
library(rsconnect)
rsconnect::deploySite(
  siteDir = "quarto/website",
  server = "az-connect",
  account = "kmjq089"
)
```

Or from bash:
```bash
# From /mnt/code (project root)
quarto publish connect quarto/website
```

#### Troubleshooting

- **API key expired**: Generate a new one and re-run `rsconnect::connectApiUser()`
- **Publishing fails**: Make sure `quarto render quarto/website` (from project root) completes successfully first
- **Missing plots**: Delete `_freeze/` cache and re-render from project root
- **IMPORTANT**: Always run quarto commands from `/mnt/code` (project root), not from inside `quarto/website/`
- **Git tracking warnings**: Confirm `rsconnect/` and `_publish.yml` are in `.gitignore`
- **Can't find account/server**: Make sure you're using the correct server name (`az-connect`). Run `rsconnect::accounts()` to check configured accounts.
- **Package not found errors**: Always run deployment from `/mnt/code` (project root) where renv is configured

**SECURITY**: API keys should NEVER be committed to git. The `rsconnect/` directory and `_publish.yml` are already in `.gitignore`.

### Key Configuration Files
- `_quarto.yml` - Site configuration, navigation, theme settings
- `az-theme.scss` - AstraZeneca color scheme (navy, gold, turquoise, etc.)
- `styles.css` - Custom CSS for hero section and layout
- `_freeze/` - Cache directory for executed R code (speeds up rebuilds)
- `images/` - PIONEER helmet logo and favicon

### Theme and Styling
The site uses **AZ corporate colors** defined in `az-theme.scss`:
- Primary (navy): `#003865`
- Gold: `#F0AB00`
- Turquoise: `#68D2DF`
- Plum: `#830051`
- Pink: `#D0006F`
- Platinum: `#9DB0AC`

Home page features a 2/3 text + 1/3 image layout with the PIONEER helmet logo. Favicon is a white "P" on navy circle background.

### Site Structure
```
├── index.qmd                           # Home page
├── documentation/
│   ├── onboarding-tutorial.qmd         # Getting started guide
│   └── model-specification.qmd         # Full model specification
└── analysis/
    ├── descriptive-statistics.qmd      # Data exploration
    ├── model-validation.qmd            # Posterior checks, LFO results
    └── outcome-predictions.qmd         # PFS predictions
```

### Rebuilding After Changes
If you update plot functions in `r/plot_functions.R`, you may need to clear the freeze cache:
```bash
rm -rf quarto/website/_freeze/analysis/
```

Then re-render the affected pages.

## Architecture

### Stan Module System

Stan code uses modular `#include` architecture in `stan/modules/`:
- **tr/** - Tumor regression (decrease) dynamics
- **frac/** - Growth fraction dynamics
- **init/** - Initial state modeling
- **other_events/** - Non-target progression and death events
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
- `oe_*` - Other events

### Key Stan Models
- `stan/sf-ssm-log-space.stan` - Main state-space longitudinal survival model
- `stan/sf-ssls-lfo.stan` - Leave-future-out cross-validation variant

## Key Files

### R Code
- `r/priors.R` - Prior specifications (shared)
- `r/initializers.R` - Stan model initializers (shared)
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

## Coding Guidelines

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

### Targets
- **NEVER use `tar_config_set(store = ...)`** - it changes global state and causes conflicts
- Always use explicit `store` argument: `tar_read(name, store = "path/_targets")`
- Same applies to all targets functions: `tar_meta()`, `tar_load()`, etc.

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

## Documentation

Key technical docs in `docs/`:
- `ARCHITECTURE.md` - System architecture and optimization details
- `ar1_process_noise.md` - AR(1) time-varying process noise implementation
- `OTHER_EVENTS_MODEL.md` - Other events model design
- `CHANGELOG.md` - Major changes and design decisions
