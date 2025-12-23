This repository is used by PIONEER projects for predicting PFS and ORR outcomes in small and/or ongoing trials using data from older trials. Currently, there are three active projects:

1. **Breast-01 to -04**: Predicting outcomes from one trial using the other three.
2. **Endometrial to LUNG**: Predicting the LUNG trial outcomes from the Endometrial trial.
3. **SCLC-01**: Longitudinal tumor analysis with PFS and ORR predictions.

We predict outcomes using a Bayesian multilevel piece-wise constant proportional hazard survival model for PFS and, similarly, a competing risks model for confirmed response status (which we use as predictor in the PFS model). The tumor dynamics are modeled using a state-space longitudinal survival (SSLS) model with modular parameter organization.

# Repo Structure

``` bash
.
├── targets/                                                          # {targets} pipeline definitions
│   ├── sclc_targets.R                                            # SCLC project pipeline
│   ├── endometrial_to_lung_targets.R                                    # Endometrial to LUNG pipeline
│   └── test_targets.R                                                # Testing pipeline
├── quarto/                                                           # Analysis notebooks
│   ├── sclc-tumor-analysis.qmd                                   # SCLC tumor analysis
│   ├── confirmed-response-pfs.qmd                                    # Breast notebook
│   └── endometrial-to-lung.qmd                                         # Endometrial to LUNG notebook
├── r/                                                                # R functions and utilities
│   ├── sclc/                                                     # SCLC-specific functions
│   ├── breast/                                               # Breast-specific functions
│   ├── endometrial-to-lung/                                            # Endometrial to LUNG-specific functions
│   ├── entimice/                                                     # Scripts to download analysis datasets
│   ├── data_preparation_pipeline/                                    # Data preparation utilities
│   ├── initializers.R                                                # Stan model initializers
│   ├── priors.R                                                      # Prior specifications
│   ├── posterior.R                                                   # Posterior extraction utilities
│   ├── plot_functions.R                                              # Visualization functions
│   ├── table_functions.R                                             # Table generation ({gt})
│   └── util.R                                                        # Common utilities
└── stan/                                                             # Stan statistical models
    ├── ssls/                                                         # State-space longitudinal survival (SSLS) models
    │   ├── modules/                                                  # Modular parameter organization
    │   │   ├── tr/                                                   # Total rate module
    │   │   ├── frac/                                                 # Fraction mix module
    │   │   ├── init/                                                 # Initial proportions module
    │   │   └── other_events/                                         # Other events (non-target PD, death)
    │   ├── sf-ssm-log-space.stan                                     # Main SSLS model
    │   ├── sf-ssls-lfo.stan                                          # Leave-future-out cross-validation
    │   └── legacy/                                                   # Legacy monolithic implementations
    ├── tumor/                                                        # Tumor-specific model components
    ├── crcr/                                                         # Confirmed response competing risks
    ├── baseline_hazard/                                              # GP PFS baseline hazard
    ├── pfs-confirmed-response/                                       # PFS with confirmed response
    ├── recruit/                                                      # Recruitment modeling
    ├── breast/                                               # Breast-specific files
    └── util.stan                                                     # General utility functions

# Documentation

For detailed technical documentation:

- **[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)** - System architecture, hierarchical parameter design, and Stan optimization
- **[`docs/OTHER_EVENTS_MODEL.md`](docs/OTHER_EVENTS_MODEL.md)** - Other events model (non-target PD, death) design and implementation
- **[`docs/CHANGELOG.md`](docs/CHANGELOG.md)** - Major changes and design decisions
- **[`docs/CODEOWNERS`](docs/CODEOWNERS)** - Code ownership and review requirements
- **`sld_state_space_model.md`** - State space model documentation

# Workflow

1. Packages
    - Make sure the package {renv} is installed on your system.
    - Run `renv::restore()`.
    - If you encounter installation problems for a particular package, run `renv::install("<package name>@<version>", type = "source")` where the version is same as located in the "renv.lock" file.
3. Install CmdStan by running `cmdstanr::install_cmdstan()`.
4. Run the project-specific script (below) to download your data.
5. Make analysis pipeline.
    - Make sure the correct TAR_PROJECT value is used for your project (see projects in "_targets.yaml").
    - Execute `targets::tar_make()`.
6. Build the project-specific notebook.

# Projects

## SCLC-01 Longitudinal Tumor Analysis

* {targets} project: `sclc`
* Targets file: `targets/sclc_targets.R`
* Notebook: `quarto/sclc-tumor-analysis.qmd`
* Description: Longitudinal tumor analysis with state-space modeling, including PFS and ORR predictions

## Breast-01 to -04

* Data script: `r/entimice/download-breast-entimice-data.R`
* {targets} project: (in development)
* Notebook: `quarto/confirmed-response-pfs.qmd`
* Description: Cross-validation predictions across Breast trials

## Predicting LUNG from Endometrial

* Data scripts: `r/entimice/download-endometrial-entimice-data.R` and `r/entimice/download-lung-entimice-data.R`
* {targets} project: `endometrial_to_lung`
* Targets file: `targets/endometrial_to_lung_targets.R`
* Notebook: `quarto/endometrial-to-lung.qmd`
* Published URL: https://rstudio-connect.seml.scp.astrazeneca.net/endometrial-to-lung/endometrial-to-lung.html
* Description: Predicting LUNG trial outcomes using Endometrial trial data
