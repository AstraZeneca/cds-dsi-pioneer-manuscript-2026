# Publication Website Design

**Date:** 2026-05-27
**Scope:** Add a Quarto website for the publication pipeline, mirroring the sclc site structure with minimal duplication.

---

## Context

The codebase hosts two existing Quarto websites:
- `quarto/sclc/website/` — full SCLC-01 forecasting site (analysis + documentation)
- `quarto/pioneer/website/` — Pioneer analysis site

Both source shared infrastructure from `quarto/_shared/` (colors, base `_setup.R`) and follow the same pattern: own `_quarto.yml`, own `_setup.R` sourcing `_shared/_setup.R`, own `az-theme.scss` importing `_shared/az-colors.scss`.

The publication pipeline (`targets/publication_targets.R`, `TAR_PROJECT=publication`) runs the same SF-SSM Stan model on a combined 497-patient dataset with 6 covariates but without PDL1 stratification, ctDNA covariates, or conditional groups. Its target namespace uses plain names (`tumor_ssls_km_rvar_posterior`) rather than the sclc `_ctdna_jan26` suffix convention.

---

## Goals

- Results-only internal website (no documentation section — sclc site is canonical docs)
- Minimize duplication: same pattern as pioneer, reuse existing R plot/table functions
- Add missing pipeline targets so all website pages have data

---

## 1. Website Structure

New directory: `quarto/publication/website/`

```
quarto/publication/website/
├── _quarto.yml           # site config, navbar, theme refs
├── _metadata.yml         # date-modified: last-modified
├── _setup.R              # sources _shared/_setup.R + pub store + sclc plot/table fns
├── az-theme.scss         # @import "../../_shared/az-colors"; (identical to pioneer)
├── styles.css            # copied from sclc (hero, cards, figure constraints)
├── index.qmd             # landing page
└── analysis/
    ├── _metadata.yml
    ├── data-exploration.qmd
    ├── survival-analysis.qmd
    ├── overall-survival.qmd
    ├── clinical-endpoints.qmd
    ├── parameter-inference.qmd
    ├── patient-dynamics.qmd
    └── diagnostics.qmd
```

No `documentation/` folder. Navbar right side links to the sclc site for docs.

### `_quarto.yml` nav structure

```yaml
website:
  title: "PIONEER: Publication"
  navbar:
    left:
      - icon: house-fill / Data / Analysis menu
    right:
      - icon: book  # links to sclc site
      - icon: github  # pub repo
```

Analysis menu mirrors sclc: Data | Results (Clinical Endpoints, Survival, OS) | Model Details (Parameter Inference, Patient Dynamics) | Validation (Diagnostics).

---

## 2. R Infrastructure

### `r/publication/target_stores.R` (new file)

Defines a single store path for the publication pipeline:

```r
pub_store <- file.path(output_path, "publication", Sys.getenv("TAR_RUN", "main"), "_targets")
```

### `_setup.R`

```r
here::i_am("quarto/publication/website/_setup.R")
source(here::here("quarto", "_shared", "_setup.R"))
source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "publication", "target_stores.R"))
```

Reuses sclc's plot/table functions directly — publication targets produce the same data shapes (same `fit_type`, `trial`, `sample_*`/`spop_*` column conventions). No publication-specific plot functions needed at this point.

---

## 3. Missing Pipeline Targets

The following targets must be added to `targets/publication_targets.R`. All belong in the **inside the `tar_map`** (per-type: prior/posterior) section, then combined in the **Combined prior + posterior** section at the bottom.

### 3a. Inside `tar_map` (add alongside existing km_rvar, km_os_rvar, trial_pfs_quant)

**ORR:**
```r
tar_target(
  tumor_ssls_orr_rvar,
  tumor_ssls_draws_endpoints |>
    recover_types(select(all_analysis_data, trial)) |>
    spread_rvars(sample_target_orr[trial], spop_target_orr[trial]) |>
    mutate(fit_type = type)
)
```

**PFS-n at fixed timepoints:**
```r
tar_target(
  tumor_ssls_forecast_target_pfs_n_rvar,
  tumor_ssls_draws_endpoints |>
    recover_types(select(all_analysis_data, trial)) |>
    spread_rvars(
      sample_target_pfs_n[trial, n], spop_target_pfs_n[trial, n],
      sample_ms_pfs_n[trial, n],    spop_ms_pfs_n[trial, n],
      sample_pfs_n[trial, n],       spop_pfs_n[trial, n]
    ) |>
    left_join(pfs_timepoints_pub, by = "n") |>
    mutate(fit_type = type)
)
```

**Patient subsamples + SLD/RECIST trajectories**: In sclc these live in a `tar_map` over `event_type` (right_censored / uncensored) *nested inside* the outer `tar_map` over `type` (prior/posterior). The publication `tar_map` must gain the same inner `tar_map`. The inner targets (`state_patient_subsample`, `tumor_ssls_rep_sld_rvar`, etc.) get names like `tumor_ssls_rep_sld_rvar_right_censored_prior`.

```r
# Inside the existing tar_map(type = c("prior","posterior"), ...):
tar_map(
  tibble(
    event_type   = c("right_censored", "uncensored"),
    event_cond   = c(expr(right_censored), expr(!right_censored)),
    event_slicer = c(\(d, n) d, \(d, n) slice_sample(d, n = n)),
  ),
  names = "event_type",
  tar_target(state_patient_subsample,
    get_state_patients(all_analysis_data, by = trial,
                       cond = event_cond, slicer = event_slicer, sample_size = 40)),
  # SLD and RECIST targets follow here (see 3a below)
)
```

**SLD and RECIST trajectories** (inside the same inner `tar_map`):
```r
tar_target(tumor_ssls_rep_sld_rvar,
  get_sld(tumor_ssls_draws_sld_recist, state_patient_subsample) |> mutate(fit_type = type))
tar_target(tumor_ssls_forecast_sld_rvar,
  get_forecast_sld(tumor_ssls_draws_sld_recist, all_analysis_data,
                   state_patient_subsample, forecast_extent = extend_max_all_t) |>
  mutate(fit_type = type))
tar_target(tumor_ssls_staged_sld_rvar,
  bind_rows(obs = tumor_ssls_rep_sld_rvar |> rename(patient_sld = rep_patient_sld),
            forecast = tumor_ssls_forecast_sld_rvar |> rename(patient_sld = forecast_patient_sld),
            .id = "stage"))
tar_target(tumor_ssls_recist_rvar,
  get_recist(tumor_ssls_draws_sld_recist, state_patient_subsample) |> mutate(fit_type = type))
tar_target(tumor_ssls_forecast_recist_rvar,
  get_forecast_recist(tumor_ssls_draws_sld_recist, all_analysis_data,
                      state_patient_subsample, forecast_extent = extend_max_all_t) |>
  mutate(fit_type = type))
tar_target(tumor_ssls_staged_recist_rvar,
  bind_rows(obs = tumor_ssls_recist_rvar |> filter(!is.na(response)) |>
                    rename(recist = rep_recist),
            forecast = tumor_ssls_forecast_recist_rvar |> rename(recist = forecast_obs_recist),
            .id = "stage"))
```

### 3b. Combined prior + posterior (add alongside existing `all_tumor_ssls_*`)

```r
tar_target(all_tumor_ssls_orr_rvar,
  bind_rows(tumor_ssls_orr_rvar_prior, tumor_ssls_orr_rvar_posterior))
tar_target(all_tumor_ssls_forecast_target_pfs_n_rvar,
  bind_rows(tumor_ssls_forecast_target_pfs_n_rvar_prior,
            tumor_ssls_forecast_target_pfs_n_rvar_posterior))
tar_target(all_tumor_ssls_rep_sld_rvar,
  bind_rows(tumor_ssls_rep_sld_rvar_right_censored_prior,
            tumor_ssls_rep_sld_rvar_right_censored_posterior,
            tumor_ssls_rep_sld_rvar_uncensored_prior,
            tumor_ssls_rep_sld_rvar_uncensored_posterior))
tar_target(all_tumor_ssls_forecast_sld_rvar,
  bind_rows(tumor_ssls_forecast_sld_rvar_right_censored_prior,
            tumor_ssls_forecast_sld_rvar_right_censored_posterior,
            tumor_ssls_forecast_sld_rvar_uncensored_prior,
            tumor_ssls_forecast_sld_rvar_uncensored_posterior))
tar_target(all_tumor_ssls_staged_sld_rvar,
  bind_rows(tumor_ssls_staged_sld_rvar_right_censored_prior,
            tumor_ssls_staged_sld_rvar_right_censored_posterior,
            tumor_ssls_staged_sld_rvar_uncensored_prior,
            tumor_ssls_staged_sld_rvar_uncensored_posterior))
tar_target(all_tumor_ssls_staged_recist_rvar,
  bind_rows(tumor_ssls_staged_recist_rvar_right_censored_prior,
            tumor_ssls_staged_recist_rvar_right_censored_posterior,
            tumor_ssls_staged_recist_rvar_uncensored_prior,
            tumor_ssls_staged_recist_rvar_uncensored_posterior))
```

---

## 4. Website Pages — Target Mapping

Each page's `_setup.R` block reads from `pub_store`. Target names have no DCO suffix.

| Page | Targets read |
|------|-------------|
| `data-exploration.qmd` | `all_analysis_data`, `all_stan_data` |
| `clinical-endpoints.qmd` | `all_tumor_ssls_orr_rvar`, `all_tumor_ssls_trial_pfs_quant`, `all_tumor_ssls_forecast_target_pfs_n_rvar`, `km_trial_pfs` |
| `survival-analysis.qmd` | `all_tumor_ssls_km_rvar`, `km_trial_pfs` |
| `overall-survival.qmd` | `all_tumor_ssls_km_os_rvar`, `km_trial_os` |
| `parameter-inference.qmd` | `all_tumor_ssls_rates_rvar`, `all_tumor_ssls_coef` |
| `patient-dynamics.qmd` | `all_tumor_ssls_staged_sld_rvar`, `all_tumor_ssls_staged_recist_rvar`, `all_analysis_data` |
| `diagnostics.qmd` | `tumor_ssls_nuts_summary_posterior`, `tumor_ssls_convergence_posterior` |

---

## 5. What is NOT included (compared to sclc)

- No `comparisons/` section (no multiple DCOs yet)
- No `documentation/` section (sclc is canonical)
- No PDL1/ctDNA stratified plots (not in pub data)
- No conditional group plots (no `cond_groups` in pub pipeline)
- No LFO/model-validation page (no LFO run yet)

---

## 6. Build and publish commands

```bash
# Render (from project root)
TAR_PROJECT=publication TAR_RUN=main quarto render quarto/publication/website

# Publish
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/website", server = "az-connect", account = "kmjq089")'
```

---

## Implementation Order

1. Add missing targets to `targets/publication_targets.R`
2. Create `r/publication/target_stores.R`
3. Create `quarto/publication/website/` scaffold (`_quarto.yml`, `_metadata.yml`, `az-theme.scss`, `styles.css`, `_setup.R`)
4. Create `index.qmd`
5. Create `analysis/_metadata.yml` + 7 analysis pages
6. Update `.claude/rules/r-guidelines.md` with publication site build commands
