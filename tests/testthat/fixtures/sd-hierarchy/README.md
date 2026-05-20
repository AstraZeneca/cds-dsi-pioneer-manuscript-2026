# SD Hierarchy Bit-Exact Fixture

Files:
- `bitexact-reference-stan-data.rds` — a `stan_data` list (as produced by
  `prepare_psa_standalone_stan_data(base_stan_data, fit_data = TRUE)`)
  with `enable_sd_level_intercept_mode_tr` all-zero.
- `bitexact-reference.csv` — the CmdStan CSV output of a fit with this data on
  the `main` branch, seed=42, iter_warmup=50, iter_sampling=50, chains=1.

## Capture procedure (run once, from `main`)

1. Check out `main` in a separate worktree:

   ```bash
   git worktree add /tmp/multilevel-sd-ref main
   cd /tmp/multilevel-sd-ref
   ```

2. Source the data-prep helpers and build `stan_data` for a psa_standalone
   configuration that uses minimal fixture data. Save it:

   ```r
   source("r/pioneer/prepare_analysis_data.R")
   # ... load minimal pioneer fixture data into `base_stan_data` ...
   stan_data <- prepare_psa_standalone_stan_data(base_stan_data, fit_data = TRUE)
   saveRDS(stan_data, "/mnt/code/.worktrees/karim/multilevel-sd/tests/testthat/fixtures/sd-hierarchy/bitexact-reference-stan-data.rds")
   ```

3. Fit and save the CSV:

   ```r
   library(cmdstanr)
   model <- cmdstan_model(
     "stan/psa/pioneer.stan",
     include_paths = c("stan", "stan/psa")
   )
   fit <- model$sample(
     data = stan_data,
     seed = 42,
     iter_warmup = 50,
     iter_sampling = 50,
     chains = 1,
     parallel_chains = 1,
     refresh = 0
   )
   file.copy(
     fit$output_files()[1],
     "/mnt/code/.worktrees/karim/multilevel-sd/tests/testthat/fixtures/sd-hierarchy/bitexact-reference.csv"
   )
   ```

4. Clean up:

   ```bash
   cd /mnt/code/.worktrees/karim/multilevel-sd
   git worktree remove /tmp/multilevel-sd-ref
   git add tests/testthat/fixtures/sd-hierarchy/
   git commit -m "test(fixtures): capture bit-exact reference from main"
   ```

## Gotchas

- Bit-exactness depends on cmdstan version. If you upgrade cmdstan, re-capture.
- The test strips CSV comment lines (timestamps, PIDs) before comparing; only
  the numeric draws are compared.
- If `main` is ahead of the pre-feature state (e.g., other commits have landed),
  capture from the commit immediately BEFORE the first SD-hierarchy change on
  this branch instead of literal `main`.

## Active-path fixture

Also required for `test-sd-subhierarchy-active.R`:

- `minimal-active-stan-data.rds` — a `stan_data` list produced by
  `prepare_psa_standalone_stan_data(base_stan_data, fit_data = TRUE,
   tr_sd_intercept_modes = tibble::tribble(
     ~location_level, ~sub_level, ~mode,
     "patient",       "arm",      "re_cp"
   ))`

### Capture procedure

On THIS branch (where the feature exists):

```r
source("r/pioneer/prepare_analysis_data.R")
# ... load minimal pioneer fixture data into `base_stan_data` ...
active_modes <- tibble::tribble(
  ~location_level, ~sub_level, ~mode,
  "patient",       "arm",      "re_cp"
)
stan_data <- prepare_psa_standalone_stan_data(
  base_stan_data, fit_data = TRUE,
  tr_sd_intercept_modes = active_modes
)
saveRDS(
  stan_data,
  "tests/testthat/fixtures/sd-hierarchy/minimal-active-stan-data.rds"
)
```

Then commit:

```bash
git add tests/testthat/fixtures/sd-hierarchy/minimal-active-stan-data.rds
git commit -m "test(fixtures): capture minimal active-path stan_data for smoke test"
```

No need to separately capture a fit-result fixture — the smoke test runs the fit inline (5 iterations only).
