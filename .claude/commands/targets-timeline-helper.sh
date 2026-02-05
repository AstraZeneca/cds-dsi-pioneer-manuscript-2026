#!/bin/bash
# Targets Store Timeline Helper Script
# Shows timeline and contents of targets stores

set -e

COMMAND="${1:-status}"
ARG1="${2:-}"
ARG2="${3:-}"

# Get username
USERNAME=$(Rscript -e 'cat(Sys.getenv("DOMINO_STARTING_USERNAME"))' 2>/dev/null || echo "karim_naguib")

# Determine store name
if [ "$COMMAND" = "compare" ]; then
  STORE1="${ARG1:-main}"
  STORE2="${ARG2:-dco3}"
else
  STORE_NAME="${ARG1:-$(Rscript -e 'cat(Sys.getenv("TAR_BRANCH", "main"))' 2>/dev/null)}"
  STORE_PATH="/mnt/data/analysis-results/${USERNAME}/sclc/${STORE_NAME}/_targets"
fi

show_timeline() {
  local store_name=$1
  local store_path="/mnt/data/analysis-results/${USERNAME}/sclc/${store_name}/_targets"

  if [ ! -d "$store_path" ]; then
    echo "❌ Store not found: $store_path"
    return 1
  fi

  echo "# Targets Store Timeline: $store_name"
  echo
  echo "**Store path:** \`$store_path\`"
  echo

  Rscript -e "
  suppressPackageStartupMessages({
    library(targets)
    library(dplyr)
  })

  store <- '$store_path'
  meta <- tar_meta(store = store)

  # Overall stats
  cat('\n## Overview\n\n')
  cat('- **Total targets:**', nrow(meta), '\n')
  cat('- **First target:**', format(min(meta\$time, na.rm = TRUE), '%Y-%m-%d %H:%M'), '\n')
  cat('- **Last target:**', format(max(meta\$time, na.rm = TRUE), '%Y-%m-%d %H:%M'), '\n')

  # Data cuts available
  cat('\n## Data Cuts Available\n\n')
  data_targets <- meta %>%
    filter(grepl('all_analysis_data', name)) %>%
    mutate(
      data_cut = gsub('all_analysis_data_', '', name)
    ) %>%
    select(data_cut, time) %>%
    arrange(desc(time))

  if (nrow(data_targets) > 0) {
    cat('| Data Cut | Last Updated |\n')
    cat('|----------|-------------|\n')
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
  cat('\n## Key Milestones by Data Cut\n')

  # Find unique data cut suffixes
  suffixes <- meta %>%
    filter(grepl('ctdna_(jan26|aug25?|apr25?)', name)) %>%
    pull(name) %>%
    gsub('.*_(ctdna_[a-z0-9]+).*', '\\\\1', .) %>%
    unique() %>%
    sort(decreasing = TRUE)

  for (suffix in suffixes) {
    cat(sprintf('\n### %s\n\n', toupper(suffix)))

    # Define milestone target patterns
    milestone_patterns <- c(
      'all_analysis_data',
      'all_stan_data',
      'tumor_ssls_fit',
      'all_tumor_ssls_forecast_recist_rvar',
      'all_tumor_ssls_cond_forecast_target_pfs_n_rvar'
    )

    milestones <- data.frame()
    for (pattern in milestone_patterns) {
      target_name <- paste0(pattern, '_', suffix)
      target_meta <- meta %>% filter(name == target_name)

      if (nrow(target_meta) > 0) {
        milestone_name <- case_when(
          pattern == 'all_analysis_data' ~ 'Data prep',
          pattern == 'all_stan_data' ~ 'Stan data',
          pattern == 'tumor_ssls_fit' ~ 'Model fit',
          pattern == 'all_tumor_ssls_forecast_recist_rvar' ~ 'Forecasts',
          pattern == 'all_tumor_ssls_cond_forecast_target_pfs_n_rvar' ~ 'Conditional forecasts',
          TRUE ~ pattern
        )
        milestones <- bind_rows(milestones, data.frame(
          milestone = milestone_name,
          time = target_meta\$time[1]
        ))
      }
    }

    if (nrow(milestones) > 0) {
      milestones <- milestones %>% arrange(time)
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
    cat('Reading from target:', target_name, '\n\n')

    tryCatch({
      data <- tar_read_raw(target_name, store = store)

      if ('variable' %in% names(data)) {
        variables <- unique(data\$variable)
        for (v in variables) {
          groups <- unique(data[data\$variable == v, 'cond_group_name', drop = TRUE])
          cat(sprintf('- **%s**: %s\n', v, paste(groups, collapse = ', ')))
        }
      }
    }, error = function(e) {
      cat('Could not read conditional forecast data.\n')
    })
  } else {
    cat('No conditional forecast targets found.\n')
  }
  " 2>&1 | grep -v "^\\[conflicted\\]"
}

compare_stores() {
  local store1=$1
  local store2=$2

  echo "# Store Comparison"
  echo

  for store_name in "$store1" "$store2"; do
    store_path="/mnt/data/analysis-results/${USERNAME}/sclc/${store_name}/_targets"

    echo "## $store_name"
    echo

    if [ ! -d "$store_path" ]; then
      echo "❌ Store not found"
      echo
      continue
    fi

    Rscript -e "
    suppressPackageStartupMessages({
      library(targets)
      library(dplyr)
    })

    store <- '$store_path'
    meta <- tar_meta(store = store)

    # Data cuts
    data_cuts <- meta %>%
      filter(grepl('all_analysis_data_ctdna', name)) %>%
      mutate(data_cut = gsub('all_analysis_data_ctdna_', '', name)) %>%
      arrange(desc(data_cut)) %>%
      pull(data_cut)

    cat('- **Data cuts:**', paste(data_cuts, collapse = ', '), '\n')

    # Latest timestamp
    cat('- **Last updated:**', format(max(meta\$time, na.rm = TRUE), '%Y-%m-%d %H:%M'), '\n')

    # Key targets count
    cat('- **Total targets:**', nrow(meta), '\n')

    # Check for conditional forecasts
    has_cond <- any(grepl('cond_forecast', meta\$name))
    cat('- **Conditional forecasts:**', if(has_cond) '✓' else '✗', '\n')
    " 2>&1 | grep -v "^\\[conflicted\\]"

    echo
  done
}

# Main dispatch
case "$COMMAND" in
  status|"")
    show_timeline "${STORE_NAME}"
    ;;
  compare)
    compare_stores "$STORE1" "$STORE2"
    ;;
  *)
    # Assume it's a store name
    show_timeline "$COMMAND"
    ;;
esac
