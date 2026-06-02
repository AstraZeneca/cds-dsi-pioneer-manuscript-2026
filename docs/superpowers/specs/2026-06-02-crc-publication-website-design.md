# CRC Publication Website Design

**Date:** 2026-06-02
**Branch:** karim/crc
**Status:** Approved

## Goal

Add CRC results to the publication website. Currently only SCLC results are
presented. Both diseases should be independent subsites (separate URLs, separate
navigation) that start from the same page template but are free to diverge as
CRC-specific analysis needs emerge.

## Directory Layout

```
quarto/publication/
  sclc/                    ← rename from current website/
    _quarto.yml
    _setup.R               ← sets pub_disease <- "sclc"; sources _setup_common.R
    index.qmd
    az-theme.scss
    styles.css
    images/
    analysis/
      clinical-endpoints.qmd
      data-exploration.qmd
      diagnostics.qmd
      overall-survival.qmd
      parameter-inference.qmd
      patient-dynamics.qmd
      survival-analysis.qmd
      _metadata.yml
  crc/                     ← new subsite, copied from sclc/ as starting template
    _quarto.yml
    _setup.R               ← sets pub_disease <- "crc"; sources _setup_common.R
    index.qmd
    az-theme.scss
    styles.css
    images/
    analysis/
      clinical-endpoints.qmd   ← stub initially
      data-exploration.qmd     ← stub initially
      diagnostics.qmd          ← stub initially
      overall-survival.qmd     ← stub initially
      parameter-inference.qmd  ← stub initially
      patient-dynamics.qmd     ← stub initially
      survival-analysis.qmd    ← stub initially
      _metadata.yml
  _setup_common.R           ← shared: pub_tar_read() helper
```

The existing `quarto/publication/website/` directory is renamed to
`quarto/publication/sclc/` — no content changes, only path change.

## Shared Infrastructure: `_setup_common.R`

All per-disease `.qmd` pages call `pub_tar_read("target_name")` instead of
`tar_read(target_name_sclc, store = pub_store)`. The helper appends the disease
suffix automatically:

```r
# quarto/publication/_setup_common.R
pub_tar_read <- function(name, store = pub_store) {
  targets::tar_read(
    !!rlang::sym(paste0(name, "_", pub_disease)),
    store = store
  )
}
```

This is necessary because `tar_map(disease_map, names = "disease", ...)` in
`publication_targets.R` suffixes every target with `_sclc` or `_crc`. The
helper centralises the suffix logic so pages contain no disease-specific code.

## Per-Subsite `_setup.R`

Each subsite's `_setup.R` sets `pub_disease` before sourcing the rest:

```r
# quarto/publication/sclc/_setup.R
pub_disease <- "sclc"
here::i_am("quarto/publication/sclc/_setup.R")
source(here::here("quarto", "_shared", "_setup.R"))
source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))
source(here::here("r", "publication", "target_stores.R"))
source(here::here("quarto", "publication", "_setup_common.R"))
```

The CRC version is identical except `pub_disease <- "crc"` and the `i_am()`
path.

## `_quarto.yml` Changes

Both subsites get updated titles:
- SCLC: `"PIONEER: Publication — SCLC"`
- CRC: `"PIONEER: Publication — CRC"`

The `r-guidelines.md` rule for publishing is updated to add:
```
quarto render quarto/publication/sclc
quarto render quarto/publication/crc
```

## CRC Stub Pages

CRC pages render immediately with a clear pending notice before the fit is
available:

```markdown
> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

Pages do not attempt `tar_read()` while stubbed. The freeze cache means stub
pages render at zero cost and do not block the SCLC site.

## Page Modifications

All existing `.qmd` pages in the SCLC subsite change their `tar_read()` calls
to the `pub_tar_read()` helper. Note: the current `website/` pages were built
before `tar_map` was introduced and read *unsuffixed* target names that no
longer exist in the store. The migration to `pub_tar_read()` fixes this at the
same time:

```r
# Before (reads unsuffixed name — broken since tar_map was introduced)
orr_data <- targets::tar_read(all_tumor_ssls_orr_rvar, store = pub_store)

# After (reads all_tumor_ssls_orr_rvar_sclc via helper)
orr_data <- pub_tar_read("all_tumor_ssls_orr_rvar")
```

CRC pages start as stubs and adopt the same pattern when activated.

## Deployment

Two `rsconnect::deploySite()` calls, two separate RStudio Connect items:

```r
rsconnect::deploySite(siteDir = "quarto/publication/sclc", server = "az-connect", account = "kmjq089")
rsconnect::deploySite(siteDir = "quarto/publication/crc",  server = "az-connect", account = "kmjq089")
```

The existing publication Connect item (`fdee3ad4-616a-42b4-b258-47dfb14ba50f`)
becomes the SCLC site (path changes, content item ID stays). A new Connect item
is created for CRC when it is deployed for the first time.

## Out of Scope

- Cross-disease comparison views
- Shared navigation between the two subsites
- CRC-specific new analysis pages (added later as disease-specific content diverges)
