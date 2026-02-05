# Targets Store Timeline Skill

Use this skill to view timeline summaries of targets stores in the SCLC project.

## Usage

Invoke with `/targets-timeline` followed by optional arguments:
- `/targets-timeline` - Show timeline for current store (based on TAR_BRANCH env var)
- `/targets-timeline <store-name>` - Show timeline for specific store (e.g., `dco3`, `main`, `lfo`)
- `/targets-timeline compare <store1> <store2>` - Compare timelines across stores

## Finding Stores

Store paths follow the pattern:
```
/mnt/data/analysis-results/<username>/sclc/<store-name>/_targets
```

Common stores:
- `main` - Main analysis branch
- `dco3` - DCO3 data cut analysis
- `lfo` - Leave-future-out cross-validation
- `all-trials` - All trials sensitivity analysis

## Key Target Categories

| Category | Pattern | Description |
|----------|---------|-------------|
| Data prep | `all_analysis_data_*` | Prepared analysis datasets |
| Stan data | `all_stan_data_*` | Stan-formatted data structures |
| Model fit | `tumor_ssls_fit_*` | Fitted Stan models |
| Forecasts | `all_tumor_ssls_forecast_*` | Patient-level predictions |
| Conditional forecasts | `all_tumor_ssls_cond_forecast_*` | Subgroup predictions |

## Data Cut Suffixes

| Suffix | Description | Typical Date |
|--------|-------------|--------------|
| `_ctdna_jan26` | January 2026 data cut with ctDNA | Latest |
| `_ctdna_aug` | August 2025 data cut with ctDNA | |
| `_ctdna_apr` | April 2025 data cut with ctDNA | |
| `_no_ctdna_*` | Sensitivity: no ctDNA covariates | |
| `_ctdna_only_*` | Sensitivity: only ctDNA covariates | |

## Timeline Script

```bash
#!/bin/bash
# Show timeline for a targets store

STORE_NAME="${1:-$(Rscript -e 'cat(Sys.getenv("TAR_BRANCH", "main"))' 2>/dev/null)}"
STORE_PATH="/mnt/data/analysis-results/$(Rscript -e 'cat(Sys.getenv("DOMINO_STARTING_USERNAME"))' 2>/dev/null)/sclc/${STORE_NAME}/_targets"

if [ ! -d "$STORE_PATH" ]; then
  echo "Error: Store not found at $STORE_PATH"
  exit 1
fi

echo "# Targets Store Timeline: $STORE_NAME"
echo
echo "Store path: \`$STORE_PATH\`"
echo

# Get timeline using R
Rscript -e "
library(targets)
library(dplyr)

store <- '$STORE_PATH'
meta <- tar_meta(store = store)

# Overall stats
cat('\n## Overview\n\n')
cat('Total targets:', nrow(meta), '\n')
cat('First target:', format(min(meta\$time, na.rm = TRUE), '%Y-%m-%d %H:%M'), '\n')
cat('Last target:', format(max(meta\$time, na.rm = TRUE), '%Y-%m-%d %H:%M'), '\n')

# Data cuts available
cat('\n## Data Cuts Available\n\n')
data_targets <- meta %>%
  filter(grepl('all_analysis_data_ctdna', name)) %>%
  mutate(
    data_cut = gsub('all_analysis_data_ctdna_', '', name),
    data_cut = gsub('all_analysis_data_', '', data_cut)
  ) %>%
  select(data_cut, time) %>%
  arrange(desc(time))

if (nrow(data_targets) > 0) {
  cat('| Data Cut | Timestamp |\n')
  cat('|----------|----------|\n')
  for (i in 1:nrow(data_targets)) {
    cat(sprintf('| %s | %s |\n',
      data_targets\$data_cut[i],
      format(data_targets\$time[i], '%Y-%m-%d %H:%M')
    ))
  }
} else {
  cat('No data prep targets found.\n')
}

# Key milestones for each data cut
cat('\n## Key Milestones\n\n')

# Find unique data cut suffixes
suffixes <- meta %>%
  filter(grepl('ctdna_(jan26|aug|apr)', name)) %>%
  pull(name) %>%
  gsub('.*_(ctdna_[^_]+).*', '\\\\1', .) %>%
  unique()

for (suffix in suffixes) {
  cat(sprintf('\n### %s\n\n', toupper(suffix)))

  milestones <- meta %>%
    filter(grepl(suffix, name)) %>%
    filter(name %in% c(
      paste0('all_analysis_data_', suffix),
      paste0('all_stan_data_', suffix),
      paste0('tumor_ssls_fit_', suffix),
      paste0('all_tumor_ssls_forecast_recist_rvar_', suffix),
      paste0('all_tumor_ssls_cond_forecast_target_pfs_n_rvar_', suffix)
    )) %>%
    mutate(
      milestone = case_when(
        grepl('all_analysis_data', name) ~ 'Data prep',
        grepl('all_stan_data', name) ~ 'Stan data',
        grepl('tumor_ssls_fit', name) ~ 'Model fit',
        grepl('forecast_recist', name) ~ 'Forecasts',
        grepl('cond_forecast', name) ~ 'Conditional forecasts',
        TRUE ~ 'Other'
      )
    ) %>%
    select(milestone, time) %>%
    arrange(time)

  if (nrow(milestones) > 0) {
    cat('| Milestone | Timestamp |\n')
    cat('|-----------|----------|\n')
    for (i in 1:nrow(milestones)) {
      cat(sprintf('| %s | %s |\n',
        milestones\$milestone[i],
        format(milestones\$time[i], '%Y-%m-%d %H:%M')
      ))
    }
  } else {
    cat('No milestones found for this data cut.\n')
  }
}

# Subgroups available (from conditional forecasts)
cat('\n## Subgroups Available\n\n')
cond_targets <- meta %>%
  filter(grepl('all_tumor_ssls_cond_forecast', name)) %>%
  arrange(desc(time)) %>%
  slice(1)

if (nrow(cond_targets) > 0) {
  target_name <- cond_targets\$name[1]
  data <- tar_read_raw(target_name, store = store)

  if ('variable' %in% names(data)) {
    variables <- unique(data\$variable)
    cat('Available subgroup variables:\n')
    for (v in variables) {
      groups <- unique(data[data\$variable == v, 'cond_group_name', drop = TRUE])
      cat(sprintf('- **%s**: %s\n', v, paste(groups, collapse = ', ')))
    }
  }
} else {
  cat('No conditional forecast targets found.\n')
}
" 2>/dev/null
```

## Comparing Stores

To compare what's available across stores:

```bash
#!/bin/bash
# Compare two stores

STORE1="${1:-main}"
STORE2="${2:-dco3}"

echo "# Store Comparison: $STORE1 vs $STORE2"
echo

for store_name in "$STORE1" "$STORE2"; do
  store_path="/mnt/data/analysis-results/$(Rscript -e 'cat(Sys.getenv("DOMINO_STARTING_USERNAME"))' 2>/dev/null)/sclc/${store_name}/_targets"

  echo "## $store_name"
  echo

  Rscript -e "
  library(targets)
  library(dplyr)

  store <- '$store_path'
  meta <- tar_meta(store = store)

  # Data cuts
  data_cuts <- meta %>%
    filter(grepl('all_analysis_data_ctdna', name)) %>%
    mutate(data_cut = gsub('all_analysis_data_ctdna_', '', name)) %>%
    pull(data_cut)

  cat('Data cuts:', paste(data_cuts, collapse = ', '), '\n')

  # Latest timestamp
  cat('Last updated:', format(max(meta\$time, na.rm = TRUE), '%Y-%m-%d %H:%M'), '\n')

  # Key targets count
  cat('Total targets:', nrow(meta), '\n')
  " 2>/dev/null

  echo
done
```

## Implementation Notes

- Use `tar_meta(store = store)` to read target metadata without loading data
- Filter by name patterns to find specific target categories
- Use `tar_read_raw()` or `tar_read()` with explicit store argument to load data
- Never use `tar_config_set()` - always pass explicit `store` argument
- Present timestamps in human-readable format with timezone
- Group by data cut suffix for clarity

## Output Format

Present as markdown with clear sections:
1. **Overview** - Total targets, first/last timestamps
2. **Data Cuts Available** - Which data cuts exist in this store
3. **Key Milestones** - Timeline for each data cut (prep → fit → forecasts)
4. **Subgroups Available** - Which conditioning variables are in conditional forecasts

Keep output concise but informative. Focus on what exists and when it was built.
