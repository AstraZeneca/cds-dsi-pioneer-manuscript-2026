# FlatIron to HISTORICAL Data Pipeline

## Overview

This pipeline transforms FlatIron real-world data (RWD) from the Pluvicto cohort into a standardized format matching the HISTORICAL clinical trial schema. It produces two output datasets:

1. **Patient-level data** (`flatiron_patient_data.csv`) - One row per patient with demographics, treatment info, survival outcomes, and baseline labs
2. **Visit-level data** (`flatiron_visit_data.csv`) - Longitudinal PSA measurements with response classifications

## Input Files

Located in `/mnt/data/Pioneer_data/`:

| File | Description | Records |
|------|-------------|---------|
| `pluvcohort.csv` | Patient-level cohort data from Brad's pipeline | 656 patients |
| `pluvictolonglabs.csv` | Longitudinal lab measurements | 23,630 records (550 patients) |

## Output Files

Located in `/mnt/data/Pioneer_data/`:

| File | Description | Records |
|------|-------------|---------|
| `flatiron_patient_data.csv` | Patient-level data (HISTORICAL format) | 656 patients, 58 columns |
| `flatiron_visit_data.csv` | Visit-level PSA data | ~1,600 visits, 14 columns |

## Usage

```r
# Source the pipeline
source("r/data_pipelines/flatiron_to_historical_pipeline.R")

# Run with default paths
results <- run_flatiron_to_historical_pipeline(
  output_dir = "/mnt/data/Pioneer_data"
)

# Access results
patient_data <- results$patient_data
visit_data <- results$visit_data
```

### Custom paths

```r
results <- run_flatiron_to_historical_pipeline(
  pluvcohort_path = "/path/to/pluvcohort.csv",
  longlabs_path = "/path/to/pluvictolonglabs.csv",
  output_dir = "/path/to/output/"
)
```

## Patient Data Schema

### Core Identification & Treatment

| Column | Description | Source |
|--------|-------------|--------|
| `studyid` | Study identifier | Hardcoded: "FLATIRON_PLUVICTO" |
| `usubjid` | Unique patient ID | `patientid` |
| `trtsdt` | Treatment start date | `indexstartdt` (first Pluvicto dose) |
| `trtedt` | Treatment end date | `indexenddt` |
| `treatment_end_week` | Week of treatment end | Calculated |
| `calendar_week` | Calendar week from first patient | Calculated |
| `arm` | Treatment arm | "Pluvicto monotherapy" or "Pluvicto + other" |

### Demographics

| Column | Description | Source |
|--------|-------------|--------|
| `age` | Age at index | `ageatindex` |
| `age_group` | Categorized age | <18, 18-40, 40-65, 65-75, >75 |
| `sex` | Sex | `gender` |
| `race` | Race | `race` |
| `ecogbl` | Baseline ECOG status | `ecogvalue_index` |

### Survival Outcomes

| Column | Description | Source |
|--------|-------------|--------|
| `death` | Death indicator | `deceased` |
| `death_week` | Week of death | Calculated from `deathdt` |
| `os_time` | Overall survival time (weeks) | `deathdt` or `lkfdt` |
| `os_cnsr` | OS censor flag (1=event, 0=censored) | Derived |
| `progressor` | PSA progression indicator | From visit data (RISE response) |
| `progress_week` | Week of first PSA progression | From visit data |
| `pfs` | PFS lower bound (last clean assessment) | From visit data |
| `pfs_cnsr` | PFS censor flag | Derived |
| `right_censored` | Right censored (no death, no progression) | Derived |

### Disease Characteristics

| Column | Description | Source |
|--------|-------------|--------|
| `histology` | Tumor histology | "Adenocarcinoma" (prostate) |
| `stage` | Disease stage | "De novo" or "Recurrent" |
| `liver_mets` | Liver metastases | `mets_liver` |
| `brain_mets` | Brain metastases | `mets_brain` |
| `bone_mets` | Bone metastases | `mets_bone` |
| `lung_mets` | Lung metastases | `mets_lung` |
| `lymph_mets` | Lymph node metastases | `mets_lymphnode` |
| `visceral_mets` | Visceral metastases summary | Derived |

### Treatment History

| Column | Description | Source |
|--------|-------------|--------|
| `prev_lines` | Prior lines of therapy | `priorlotnum` |
| `chemo_flag` | Prior chemotherapy | `priortaxane` |
| `chemo_naive` | Chemotherapy naive | Inverse of `priortaxane` |
| `first_line` | First-line treatment | `priorlotnum == 0` |
| `prior_arpi` | Prior ARPI | `priorarpi` |
| `prior_taxane_count` | Number of prior taxanes | `priortaxanecount` |
| `prior_arpi_count` | Number of prior ARPIs | `priorarpicount` |

### Baseline Labs

| Column | Description | Source |
|--------|-------------|--------|
| `baseline_albumin` | Albumin (g/dL) | `blrslt_albumin` |
| `baseline_alp` | Alkaline Phosphatase (IU/L) | `blrslt_alp` |
| `baseline_alt` | ALT (IU/L) | `blrslt_alanine` |
| `baseline_ast` | AST (IU/L) | `blrslt_aspar` |
| `baseline_bilirubin` | Bilirubin (mg/dL) | `blrslt_bili` |
| `baseline_gfr` | eGFR (ml/min/1.73m²) | `blrslt_gfr` |
| `baseline_hemoglobin` | Hemoglobin (g/dL) | `blrslt_heme` |
| `baseline_ldh` | LDH (IU/L) | `blrslt_ldh` |
| `baseline_leukocytes` | WBC (×10⁹/L) | `blrslt_leuk` |
| `baseline_neutrophils` | Neutrophils (×10⁹/L) | `blrslt_neut` |
| `baseline_platelets` | Platelets (×10⁹/L) | `blrslt_plat` |
| `baseline_psa` | PSA (ng/mL) | `blrslt_psa` |
| `baseline_nlr` | Neutrophil-to-lymphocyte ratio | Calculated |

## Visit Data Schema

| Column | Description |
|--------|-------------|
| `studyid` | Study identifier |
| `usubjid` | Patient ID |
| `visitnum` | Visit number (0 = baseline) |
| `visit_date` | Date of PSA measurement |
| `ady` | Analysis day (days from treatment start) |
| `day` | Study day |
| `week` | Study week |
| `psa_value` | PSA value (ng/mL) |
| `psa_pct_change` | Percent change from baseline |
| `response` | PSA response category (see below) |
| `is_baseline` | Baseline indicator |
| `baseline_psa` | Baseline PSA value |
| `bestchngcat` | Best change category |
| `worstchngcat` | Worst change category |

### PSA Response Categories

| Response | Definition |
|----------|------------|
| `BL` | Baseline measurement |
| `PSA90` | ≥90% decline from baseline (deep response) |
| `PSA50` | ≥50% decline from baseline (confirmed response) |
| `DECLINE` | Any decline (<50%) |
| `STABLE` | Change between -0% and +25% |
| `RISE` | ≥25% increase (PSA progression) |

## Key Conventions

### Treatment Start Date
- `trtsdt` = `indexstartdt` = **first Pluvicto dose date** (not line start date)
- This is derived from `startdt_pluv` in the upstream pipeline

### PSA Progression
- In prostate cancer, PSA progression (≥25% rise) is analogous to RECIST PD in solid tumors
- The `progressor` flag is set when any visit has `response == "RISE"`

### Interval Censoring
- `progress_week`: Upper bound (week of PSA rise detection)
- `pfs`: Lower bound (last assessment week before PSA rise)
- `interval_censored`: Gap between bounds (`progress_week - pfs - 1`)

## Pipeline Structure

```
r/data_pipelines/
├── flatiron_to_historical_pipeline.R   # Main pipeline
├── historical_utils.R                   # Utility functions
├── helper_functions.r               # General helpers
└── README.md                        # This file
```

### Main Functions

| Function | Description |
|----------|-------------|
| `run_flatiron_to_historical_pipeline()` | Main execution function |
| `generate_patient_data()` | Creates patient-level data |
| `generate_visit_data_from_longlabs()` | Creates visit-level data from longitudinal labs |
| `update_patient_pfs()` | Updates patient PFS metrics from visit data |

## Related Documentation

- [FLATIRON_TO_HISTORICAL_MAPPING.md](../../docs/FLATIRON_TO_HISTORICAL_MAPPING.md) - Detailed schema mapping
- [HISTORICAL_DATA_LINEAGE.md](../../docs/HISTORICAL_DATA_LINEAGE.md) - Target schema definition
- Brad's pipeline: `brads_code/d*.Rmd` - Upstream data processing

## Changelog

- **2026-02-11**: Initial pipeline creation
  - Patient-level data generation from pluvcohort
  - Visit-level data from pluvictolonglabs.csv
  - PSA response classification (PSA90/PSA50/DECLINE/STABLE/RISE)
  - PFS calculation with interval censoring support
