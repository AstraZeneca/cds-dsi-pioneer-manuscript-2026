This project uses Breast01 to -Breast04 to investigate whether early assessments for tumor sizes can be used for a meaningful prognostic prediction of progression-free survival (PFS). We model PFS using a Bayesian hierarchical piece-wise constant proportional hazard model.

This the directory structure for this project is:

``` bash
quarto/
├── deck.qmd                              # revealjs slides deck. Not longer using this, but it still has some useful plots.
└── early-predict-bc-survival.qmd         # This is the primary notebook in which all analysis is done.
stan/
├── base_data.stan                        # This is the data section to be shared by all models.
├── pfs_functions.stan                    # Stan functions to be used in the PFS model.
├── pfs.stan                              # The main PFS model file.
├── tumor_data.stan                       # This is the tumor model settings data, useful for jointly fitting PFS and tumor models.
├── tumor_model.stan                      # The code in the model section of the tumor model.
├── tumor_parameters.stan                 # Parameters used in the tumor model.
├── tumor.stan                            # Main tumor model file (standalone).
├── tumor_transformed_data.stan           # Tranformed data section part for the tumor model.
├── tumor_transformed_parameters.stan     # Transformed parameters section for the tumor model.
└── util.stan                             # Utility functions to be shared by all models.
r/
├── download-entimice-data.R              # Script to download data.
├── sbc.R                                 # Standalone script to run simulation-based calibration.
├── sbc.sh                                # SLURM script to run sbc.R.
└── util.R                                # Shared functions.
```

These are the steps needed to conduct this analysis:

1. Use `r/download-entimice-data.R` to download the studies data.
2. Render or execute the code in `quarto/early-predict-bc-survival.qmd`. This file is divided into two main sections:
    i. A tumor-size dynamics model using Gaussian processes.
    ii. The main survival model conditional on early tumor-size assessments.
3. There are a few supplementary plots also generated in `quarto/deck.qmd`.

## Code Review

To review the code here are the most critical places to look:

* In `quarto/early-predict-bc-survival.qmd`:  
    - The `load-trial-data`, `prepare-analysis-data`, and `pfs-stan-data` labeled code blocks. There will be a number of plots generated after those blocks that summarize the data.
    - To see how our primary model is fit, see the code block `two-tumor-pfs-fit`. Analysis of the results are in the following code blocks.
* The main survival model is in `stan/pfs.stan`.
    - Various supporting Stan files are included for shared functions, data structure, and the jointly fit tumor model. The latter of which is not currently used. 
    - A bit of work is done in the `transformed data` section to prepare the data for analysis.
    - While the `parameters` section has all the core random variables, its the `transformed parameters` that does all the heavy lifting with calculating the interval specific conditional probability of disease progression.
    - You will notice various variables ending with "_raw": these are used in non-centered hierarchical modeling (non-centering is a strategy to avoid divergent transitions). The actual random variables are calculated in the `transformed parameters` section (dropping the "_raw" suffix).
    - The `generated quantities` section is where the main outcomes of the analysis are generated.
    - I would say that the heart of the code is in the patient-tumor level loop in the `transformed parameters` block and the `calc_pch_loglik()` function located in `stan/pfs_functions.stan`.