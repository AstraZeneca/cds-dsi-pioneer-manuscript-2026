# PIONEER — Bayesian Joint Modelling of Mechanistic Tumour Growth and Time-to-Event Endpoints

Publication repository for the PIONEER framework manuscript.

## Overview

PIONEER is a Bayesian joint modelling framework that couples a mechanistic two-component state-space submodel of longitudinal tumour size dynamics to a multistate proportional-hazard submodel for competing clinical events. The framework enables dynamic prediction of ongoing oncology trials from immature data.

**Authors:** Karim Naguib, Roger Berché, Lu Li, Antonia Bevan, Sajan Khosla, Jessica Davies, Paul Metcalfe

## Repository Structure

```
writing/              Manuscript source (Quarto)
r/                    R functions (data prep, model helpers, plotting)
  publication/        Publication-specific data preparation
stan/                 Stan statistical models
  modules/            Modular parameter components
targets/              {targets} pipeline definitions
  publication_targets.R   Main publication pipeline
quarto/               Analysis notebooks and publication websites
  publication/        Publication output (SCLC case study)
bib/                  Bibliography
data/                 Pre-computed inverse metrics for LFO warmstart
```

## Running the Pipeline

1. Install packages: `renv::restore()`
2. Install CmdStan: `cmdstanr::install_cmdstan()`
3. Set `TAR_PROJECT=publication` in environment
4. Run: `targets::tar_make()`

## Case Study

Applied to extensive-stage small-cell lung cancer (ES-SCLC) using patient-level data from Project Data Sphere (two trials, N = 497). Leave-future-out cross-validation demonstrates calibrated PFS forecasts from month 4 of enrolment and OS convergence by month 11.
