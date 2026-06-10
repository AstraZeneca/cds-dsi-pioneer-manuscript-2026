# (level, velocity) MS Coupling-Basis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a flag-gated `(level, velocity)` time-varying covariate basis to the tumour multistate hazard, replacing the legacy `(level, log-decrease-rate, log-growth-rate)` 3-feature basis when `enable_ms_velocity_basis = 1`.

**Architecture:** A single boolean Stan flag (`enable_ms_velocity_basis`) selects between two mutually-exclusive bases. When on, `n_time_varying_covar = 2` and feature 2 becomes the central-difference velocity of the latent log-burden trajectory; the surrogate-incompatible bi-exponential rates (features 2–3) are dropped. New transformed-data constants `median_velocity_obs`/`iqr_velocity_obs` (observed per-week log-SLD deltas) standardize the velocity. All coefficient machinery is already count-driven by `n_time_varying_covar`, so only the two builder files, the tumour transformed-data block, one helper function, the flag declaration, and the R config wiring change.

**Tech Stack:** Stan (CmdStan 2.39, modular `#include`), R (targets pipelines, cmdstanr), testthat + `fixed_param` Stan-function tests.

**Spec:** `docs/superpowers/specs/2026-06-09-level-velocity-coupling-design.md`

**Branch:** `karim/laplace` (already the working branch / worktree).

---

## File Structure

**Stan — model code:**
- `stan/modules/multistate/flags.stan` — declare `enable_ms_velocity_basis` (MODIFY).
- `stan/_burden.stanfunctions` — add `standardize_velocity` + a `central_difference` helper for the dense path (MODIFY).
- `stan/modules/tumor/transformed_data.stan` — compute `median_velocity_obs`/`iqr_velocity_obs` (MODIFY).
- `stan/tumor/_tumor_observed_covar_transformed_data.stan` — alias them burden-agnostically (MODIFY).
- `stan/_ms_burden_tv_covar.stan` — dense-path velocity branch (MODIFY).
- `stan/_ms_burden_inline_tv_covar.stan` — inline-path velocity branch, 3 blocks (MODIFY).

**Stan — tests:**
- `tests/testthat/stan/test_velocity_covar_all.stan` — `fixed_param` unit test for `standardize_velocity` + central difference (CREATE).
- `tests/testthat/test-velocity-covar.R` — R wrapper asserting `n_failures == 0` (CREATE).

**R — config wiring:**
- `targets/publication_targets.R` — set `enable_ms_velocity_basis` in `default_stan_data_settings`; switch the tribble to velocity mode (MODIFY).
- `targets/sclc_targets.R` — set `enable_ms_velocity_basis = 0L` at both flag sites (MODIFY).
- `targets/pioneer_targets.R` — set `enable_ms_velocity_basis = 0L` (MODIFY).

**Docs:**
- `quarto/sclc/website/documentation/multistate-specification.qmd` — re-justify the coupling as kinetics (MODIFY).

---

## Important codebase facts (read before starting)

- **Flag plumbing.** A new `int` in `stan/modules/multistate/flags.stan` becomes a required entry in the Stan **data list**. Every model that compiles the multistate module must supply it or CmdStan errors at runtime with "variable not found". The three sites that already set `enable_ms_baseline_trend_01` are the exact sites that must set `enable_ms_velocity_basis`:
  - `targets/publication_targets.R:300` (publication; via tribble)
  - `targets/sclc_targets.R:1000` (sclc full model) and `:2357` (ms-standalone priors)
  - `targets/pioneer_targets.R:602` (pioneer)
- **`n_time_varying_covar` is data**, declared at `stan/modules/multistate/data.stan:68`, set in R from `n_tv_covar`. The coefficient arrays (`time_varying_coef_01/02/03`), their hyperpriors, priors, and R initializers are ALL sized off `n_time_varying_covar` and loop index-agnostically — they need NO edit. Confirmed at `parameters.stan:43-44,75-76,149-150`; `hyperparams.stan:114-116`; `priors.stan:95`; `r/initializers_ms.R:247-267`; `r/priors.R:43-56`.
- **The cap.** `standardize_log_burden` (`stan/_burden.stanfunctions:34-40`) caps `log_burden_normalized` at `fmin(·, 10.0)` BEFORE adding baseline. Velocity differences the capped level values but gets **no cap of its own** (spec §3.1).
- **Two forecast paths.** Dense (`_ms_burden_tv_covar.stan`, process-noise ON, materializes `states_full_grid`) vs inline (`_ms_burden_inline_tv_covar.stan`, process-noise OFF — the publication setting, computes burden analytically). Both must branch on the flag.
- **Publication uses the inline path** (`enable_patient_process_noise_tr = FALSE`, `publication_targets.R:339`). The dense path is exercised by pioneer/process-noise configs, which stay on the legacy basis — but the dense path must still emit velocity correctly when the flag is on (it shares the builder).
- **Stan test idiom.** `tests/testthat/helper-stan.R` provides `test_stan_function(stan_file, data, ...)` which compiles with `include_paths=here::here("stan")` and runs 1 `fixed_param` draw. Test Stan files live in `tests/testthat/stan/`, `#include` the `.stanfunctions` under test, and accumulate `n_failures` in `generated quantities`. The R test asserts `n_failures == 0`. See `test_cif_computation_all.stan` / `test-stan-cif-computation.R`.
- **stanc binary:** `~/.cmdstan/cmdstan-2.39.0/bin/stanc` (2.39 pinned on this branch).
- **Use `fatal_error()` not `reject()`** for transformed-data validation (repo Stan guideline).

---

## Task 1: Declare the `enable_ms_velocity_basis` flag

**Files:**
- Modify: `stan/modules/multistate/flags.stan` (after line 78, the `enable_ms_03_time_varying_cov` block)

- [ ] **Step 1: Add the flag declaration**

In `stan/modules/multistate/flags.stan`, immediately after line 78 (`int<lower=0, upper=1> enable_ms_03_time_varying_cov;`), add:

```stan

// --- Time-Varying Covariate Basis Selector ---
// Selects the burden-coupling feature set for the time-varying covariates.
//   0 (default): legacy 3-feature basis [level, log_decrease_rate, log_growth_rate]
//                (n_time_varying_covar = 3). Features 2-3 are bi-exponential SSM
//                component parameters.
//   1          : (level, velocity) 2-feature basis [level, velocity]
//                (n_time_varying_covar = 2). Feature 2 = central-difference of the
//                latent log-burden trajectory. Required for the Laplace surrogate
//                (the bi-exponential rates have no analog under the quadratic
//                surrogate; velocity = d/dw of log-burden does). See
//                docs/superpowers/specs/2026-06-09-level-velocity-coupling-design.md
// The R config MUST set this together with a matching n_time_varying_covar:
//   (enable_ms_velocity_basis, n_time_varying_covar) in {(1,2), (0,3)}.
int<lower=0, upper=1> enable_ms_velocity_basis;
```

- [ ] **Step 2: Verify the flag does not yet break any model (it will — that is expected until R wiring lands)**

This flag is now a required data variable. We will NOT compile a full model yet (R must supply the flag first, Task 6+). Instead, verify the file still parses as a fragment by checking it has no syntax typo:

Run: `grep -c "enable_ms_velocity_basis" stan/modules/multistate/flags.stan`
Expected: `1`

- [ ] **Step 3: Commit**

```bash
git add stan/modules/multistate/flags.stan
git commit -m "feat(laplace): declare enable_ms_velocity_basis flag"
```

---

## Task 2: Add `standardize_velocity` + `central_difference` helpers (TDD)

**Files:**
- Create: `tests/testthat/stan/test_velocity_covar_all.stan`
- Create: `tests/testthat/test-velocity-covar.R`
- Modify: `stan/_burden.stanfunctions` (append after `standardize_log_burden`, line 41)

- [ ] **Step 1: Write the failing test (Stan side)**

Create `tests/testthat/stan/test_velocity_covar_all.stan`:

```stan
// Tests for standardize_velocity() and central_difference_row() (_burden.stanfunctions)
//
// standardize_velocity(v, median, iqr) = (v - median) / iqr
// central_difference_row(g): forward at col 1, central interior, backward at last col.
//   v[1]      = g[2] - g[1]
//   v[k]      = (g[k+1] - g[k-1]) / 2   for 1 < k < n
//   v[n]      = g[n] - g[n-1]
//
// Output: n_failures == 0. R test asserts n_failures == 0.

functions {
  #include "_burden.stanfunctions"
}

data {
  int<lower=0> dummy;
}

generated quantities {
  int n_failures = 0;
  real tol = 1e-10;

  // --- standardize_velocity ---
  if (abs(standardize_velocity(3.0, 1.0, 2.0) - 1.0) > tol) {
    print("FAIL standardize_velocity basic"); n_failures += 1;
  }
  if (abs(standardize_velocity(1.0, 1.0, 2.0) - 0.0) > tol) {
    print("FAIL standardize_velocity centered-at-median"); n_failures += 1;
  }

  // --- central_difference_row on a known parabola g(w) = w^2, w = 1..5 ---
  // g = [1, 4, 9, 16, 25]
  // v[1] = 4-1 = 3 (forward)
  // v[2] = (9-1)/2 = 4
  // v[3] = (16-4)/2 = 6
  // v[4] = (25-9)/2 = 8
  // v[5] = 25-16 = 9 (backward)
  {
    row_vector[5] g = [1, 4, 9, 16, 25];
    row_vector[5] v = central_difference_row(g);
    row_vector[5] expected = [3, 4, 6, 8, 9];
    for (k in 1:5) {
      if (abs(v[k] - expected[k]) > tol) {
        print("FAIL central_difference_row[", k, "]=", v[k], " expected ", expected[k]);
        n_failures += 1;
      }
    }
  }

  // --- central_difference equals analytic derivative for a quadratic (interior) ---
  // For g(w) = b0 + b1*w + b2*w^2, central diff at interior = b1 + 2*b2*w exactly.
  // b0=0, b1=1, b2=0.5 → g(w) = w + 0.5 w^2; g'(w) = 1 + w.
  // w=1..5: g = [1.5, 4, 7.5, 12, 17.5]; interior v should equal 1+w = [_,3,4,5,_].
  {
    row_vector[5] g = [1.5, 4, 7.5, 12, 17.5];
    row_vector[5] v = central_difference_row(g);
    if (abs(v[2] - 3.0) > tol) { print("FAIL quad interior v[2]=", v[2]); n_failures += 1; }
    if (abs(v[3] - 4.0) > tol) { print("FAIL quad interior v[3]=", v[3]); n_failures += 1; }
    if (abs(v[4] - 5.0) > tol) { print("FAIL quad interior v[4]=", v[4]); n_failures += 1; }
  }

  // --- single-column edge case: n=1 → velocity 0 (no neighbour) ---
  {
    row_vector[1] g = [7];
    row_vector[1] v = central_difference_row(g);
    if (abs(v[1] - 0.0) > tol) { print("FAIL n=1 v[1]=", v[1]); n_failures += 1; }
  }
}
```

- [ ] **Step 2: Write the failing R test**

Create `tests/testthat/test-velocity-covar.R`:

```r
library(testthat)
library(cmdstanr)
library(here)

test_that("standardize_velocity and central_difference_row produce expected values (n_failures == 0)", {
  stan_file <- here("tests", "testthat", "stan", "test_velocity_covar_all.stan")

  fit <- test_stan_function(
    stan_file = stan_file,
    data      = list(dummy = 0L),
    seed      = 42L
  )

  draws_df   <- posterior::as_draws_df(fit$draws())
  n_failures <- as.integer(draws_df$n_failures[[1]])

  expect_equal(n_failures, 0L, label = "n_failures")
})
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-velocity-covar.R")'`
Expected: FAIL — compilation error, `central_difference_row`/`standardize_velocity` not defined in `_burden.stanfunctions`.

- [ ] **Step 4: Implement the helpers**

In `stan/_burden.stanfunctions`, after the closing `}` of `standardize_log_burden` (line 41), append:

```stan

/**
 * Standardize a velocity (rate of change of log-burden) for use as a
 * time-varying covariate. Unlike standardize_log_burden there is NO cap:
 * velocity inherits the level cap already applied to the trajectory it is
 * differenced from, and capping velocity would silently force it to zero at
 * the ceiling (consecutive capped values are equal), suppressing gradient
 * information HMC needs. See the (level, velocity) coupling-basis spec §3.1.
 *
 * @param velocity_raw        Raw velocity (per-week d/dw of log-burden)
 * @param median_velocity_obs Population median of observed per-week log-burden deltas
 * @param iqr_velocity_obs    Population IQR of observed per-week log-burden deltas
 * @return Standardized velocity
 */
real standardize_velocity(real velocity_raw,
                          real median_velocity_obs,
                          real iqr_velocity_obs) {
  return (velocity_raw - median_velocity_obs) / iqr_velocity_obs;
}

/**
 * Central-difference derivative of a log-burden trajectory sampled on a unit
 * (weekly) grid. Forward difference at the first column, central in the
 * interior, backward at the last column. For a single-column input the
 * velocity is 0 (no neighbour to difference against).
 *
 * Central difference is algebraically identical to the analytic derivative of
 * a quadratic at interior points, matching the Laplace surrogate's b1 + 2*b2*w.
 *
 * @param g  row vector of log-burden values g(w) on a unit grid
 * @return row vector of velocities v(w), same length as g
 */
row_vector central_difference_row(row_vector g) {
  int n = num_elements(g);
  row_vector[n] v;
  if (n == 1) {
    v[1] = 0.0;
    return v;
  }
  v[1] = g[2] - g[1];                       // forward
  for (k in 2:(n - 1))
    v[k] = (g[k + 1] - g[k - 1]) / 2.0;     // central
  v[n] = g[n] - g[n - 1];                   // backward
  return v;
}
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-velocity-covar.R")'`
Expected: PASS (`n_failures == 0`).

- [ ] **Step 6: Commit**

```bash
git add stan/_burden.stanfunctions tests/testthat/stan/test_velocity_covar_all.stan tests/testthat/test-velocity-covar.R
git commit -m "feat(laplace): add standardize_velocity + central_difference_row helpers (TDD)"
```

---

## Task 3: Compute observed-velocity standardization constants

**Files:**
- Modify: `stan/modules/tumor/transformed_data.stan` (after line 73, the `median_log_sld_obs` block)

This computes `median_velocity_obs`/`iqr_velocity_obs` from observed per-week log-SLD deltas between consecutive **measured** visits. It is unconditional (like the level constants) so the constant exists whenever the tumour model includes the multistate module, regardless of basis.

- [ ] **Step 1: Add the velocity-constants block**

In `stan/modules/tumor/transformed_data.stan`, immediately after line 73 (the closing `}` of the `median_log_sld_obs` block, before the `// Normalize SLD by baseline` comment at line 75), insert:

```stan

// ============================================================================
// VELOCITY NORMALIZATION CONSTANTS (for the (level, velocity) coupling basis)
// ============================================================================
// Robust median/IQR of OBSERVED per-week log(SLD) deltas between consecutive
// measured visits. Anchors the scale of the latent velocity covariate
// (observed-for-scale, latent-for-signal — mirrors median_log_sld_obs). A visit
// is "measured" when sum_tumor_size > 0; deltas are normalized per-week to match
// the latent weekly-grid central difference (visits are irregularly spaced).
// Always computed; consumed only when enable_ms_velocity_basis = 1.

real median_velocity_obs;
real iqr_velocity_obs;
{
  // Upper bound on inter-visit pairs = total visits (each patient contributes
  // at most n_visits - 1 pairs; sum over patients <= sum(n_patient_visits)).
  vector[sum(n_patient_visits)] vel_all_obs;
  int n_deltas = 0;

  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    // Walk consecutive MEASURED visits; difference adjacent measured pairs.
    int prev_v = 0;  // 0 = no previous measured visit yet
    for (v in visit_start:visit_end) {
      if (sum_tumor_size[v] > 0) {
        if (prev_v > 0) {
          int dw = t_patient_visits[v] - t_patient_visits[prev_v];
          if (dw > 0) {  // guard against duplicate-week visits (dw=0)
            n_deltas += 1;
            vel_all_obs[n_deltas] =
              (log(sum_tumor_size[v]) - log(sum_tumor_size[prev_v])) * 1.0 / dw;
          }
        }
        prev_v = v;
      }
    }
  }

  if (n_deltas < 2) {
    // Degenerate: cannot form a robust IQR. Fall back to unit scale so the
    // standardized velocity is just centered raw velocity; flag loudly.
    print("WARNING: only ", n_deltas, " observed inter-visit velocity deltas; ",
          "median_velocity_obs/iqr_velocity_obs fall back to (0, 1).");
    median_velocity_obs = 0.0;
    iqr_velocity_obs = 1.0;
  } else {
    array[3] real q = quantile(vel_all_obs[1:n_deltas], {0.25, 0.5, 0.75});
    median_velocity_obs = q[2];
    iqr_velocity_obs = q[3] - q[1];
    if (iqr_velocity_obs <= 0)
      iqr_velocity_obs = 1.0;  // all deltas equal — avoid divide-by-zero
  }
}
```

- [ ] **Step 2: Syntax-check the tumour model compiles (flag still missing from data → expect a DIFFERENT error)**

Run: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan 2>&1 | head -20`
Expected: NO error mentioning `median_velocity_obs`, `iqr_velocity_obs`, `get_pos`, `quantile`, or `t_patient_visits`. (An error about `enable_ms_velocity_basis` being an unknown identifier is acceptable here only if it appears — but at this task the flag is declared (Task 1) and not yet *read*, so the model should parse cleanly. If `stanc` passes entirely, even better.)

- [ ] **Step 3: Commit**

```bash
git add stan/modules/tumor/transformed_data.stan
git commit -m "feat(laplace): compute observed per-week velocity median/iqr constants"
```

---

## Task 4: Alias velocity constants burden-agnostically

**Files:**
- Modify: `stan/tumor/_tumor_observed_covar_transformed_data.stan` (after the existing `median_log_burden_obs` aliasing, near line 48)

- [ ] **Step 1: Add the aliases**

In `stan/tumor/_tumor_observed_covar_transformed_data.stan`, find the lines:

```stan
real median_log_burden_obs = median_log_sld_obs;
real iqr_log_burden_obs = iqr_log_sld_obs;
```

Immediately after them (before the `#include "../_observed_covar_transformed_data.stan"` line), add:

```stan
// Velocity standardization aliases for the (level, velocity) coupling basis.
// The burden-agnostic names are read by _ms_burden_tv_covar.stan /
// _ms_burden_inline_tv_covar.stan when enable_ms_velocity_basis = 1. PSA does
// not alias these (it stays on the legacy 3-feature basis).
real median_velocity_burden_obs = median_velocity_obs;
real iqr_velocity_burden_obs = iqr_velocity_obs;
```

NOTE: the burden-agnostic alias names are `median_velocity_burden_obs` /
`iqr_velocity_burden_obs` (distinct from the tumour-side `median_velocity_obs`).
The builder files (Tasks 5–6) reference the `*_burden_obs` names.

- [ ] **Step 2: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan 2>&1 | head -20`
Expected: no error about `median_velocity_obs` / `median_velocity_burden_obs`.

- [ ] **Step 3: Commit**

```bash
git add stan/tumor/_tumor_observed_covar_transformed_data.stan
git commit -m "feat(laplace): alias velocity constants burden-agnostically (tumour)"
```

---

## Task 5: Dense-path velocity branch

**Files:**
- Modify: `stan/_ms_burden_tv_covar.stan:66-82` (feature 2/3 block)

The dense path materializes `states_full_grid` and computes feature 1 (level) into `ms_time_varying_covar_01[1][j]`. Features 2–3 (the rates) are at lines 66–82. When the flag is on we replace them with a single velocity feature; when off, the legacy rates stay.

- [ ] **Step 1: Replace the feature 2/3 block with a flag branch**

In `stan/_ms_burden_tv_covar.stan`, replace lines 66–82 (the two `// Feature 2` / `// Feature 3` blocks) with:

```stan
      if (enable_ms_velocity_basis) {
        // Feature 2 (velocity mode): central-difference velocity of the
        // standardized log-burden trajectory. We difference the SAME capped
        // level values used in feature 1 (log_burden_normalized is already
        // fmin(·,10)-capped), then standardize with the velocity constants.
        // No additional cap on velocity (spec §3.1).
        row_vector[max_all_t] vel_raw = central_difference_row(log_burden_normalized);
        row_vector[max_all_t] vel_std;
        for (w in 1:max_all_t)
          vel_std[w] = standardize_velocity(vel_raw[w], median_velocity_burden_obs, iqr_velocity_burden_obs);
        ms_time_varying_covar_01[2][j] = vel_std;
      } else {
        // Feature 2 (legacy): log(decrease rate)
        if (n_time_varying_covar >= 2) {
          if (enable_any_process_noise_tr) {
            ms_time_varying_covar_01[2][j] = patient_log_decrease_rate[j, states_start_col:states_end_col];
          } else {
            ms_time_varying_covar_01[2][j] = rep_row_vector(patient_log_decrease_rate[j, 1], max_all_t);
          }
        }
        // Feature 3 (legacy): log(growth rate)
        if (n_time_varying_covar >= 3) {
          if (enable_any_process_noise_tr) {
            ms_time_varying_covar_01[3][j] = patient_log_growth_rate[j, states_start_col:states_end_col];
          } else {
            ms_time_varying_covar_01[3][j] = rep_row_vector(patient_log_growth_rate[j, 1], max_all_t);
          }
        }
      }
```

NOTE on the velocity unit: `log_burden_normalized` is `log(burden/baseline)`; its central difference over the weekly grid is the per-week velocity of log-burden. The baseline is a per-patient constant, so it cancels in the difference — differencing the normalized trajectory gives the same velocity as differencing the absolute log-burden. This matches the observed per-week delta scale in Task 3.

- [ ] **Step 2: Syntax-check (flag still not supplied by R — see Task 6 note)**

Run: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan 2>&1 | head -20`
Expected: clean parse (flag is declared in Task 1; `central_difference_row`/`standardize_velocity` exist from Task 2; `median_velocity_burden_obs` exists from Task 4).

- [ ] **Step 3: Commit**

```bash
git add stan/_ms_burden_tv_covar.stan
git commit -m "feat(laplace): dense-path (level, velocity) branch"
```

---

## Task 6: Inline-path velocity branch (3 blocks)

**Files:**
- Modify: `stan/_ms_burden_inline_tv_covar.stan` (0→2 block lines 38–55; 0→3 block lines 79–96; 0→1 block lines 117–125)

The inline path has no `states_full_grid`. For velocity we evaluate the analytic burden `g(t) = log_sum_exp(s_d(t), s_g(t))` at three points `t-1, t, t+1` and central-difference. The level cap must be applied to each `g` value before differencing, matching the dense path (which differences `fmin(·,10)`-capped values). We add a small local helper inline via a `let`-style block.

**Define the capped analytic burden once per patient.** In all three blocks the burden at week `t` is:
```
g(t) = fmin(log_sum_exp(base_d - rate_d*t, base_g + rate_g*t), 10.0)
```
(`fmin(·,10.0)` mirrors `standardize_log_burden`'s cap; note the inline level feature currently passes the *uncapped* `log_sum_exp` into `standardize_log_burden`, which caps internally — so capping here reproduces the same level value the difference is taken over.)

- [ ] **Step 1: Replace the 0→2 feature 2/3 block**

In `stan/_ms_burden_inline_tv_covar.stan`, replace lines 46–55 (the `if (n_time_varying_covar >= 2)` and `if (n_time_varying_covar >= 3)` blocks inside the 0→2 loop) with:

```stan
          if (enable_ms_velocity_basis) {
            // Velocity: central difference of capped analytic log-burden at each t.
            for (t in 1:T_max) {
              real g_prev = fmin(log_sum_exp(base_d - rate_d * (t - 1), base_g + rate_g * (t - 1)), 10.0);
              real g_next = fmin(log_sum_exp(base_d - rate_d * (t + 1), base_g + rate_g * (t + 1)), 10.0);
              real vel = (g_next - g_prev) / 2.0;  // central difference, unit grid
              log_cond_surv_02[j, t] += time_varying_coef_02[2]
                * standardize_velocity(vel, median_velocity_burden_obs, iqr_velocity_burden_obs);
            }
          } else {
            if (n_time_varying_covar >= 2) {
              real contrib = time_varying_coef_02[2] * patient_log_decrease_rate[j, 1];
              for (t in 1:T_max)
                log_cond_surv_02[j, t] += contrib;
            }
            if (n_time_varying_covar >= 3) {
              real contrib = time_varying_coef_02[3] * patient_log_growth_rate[j, 1];
              for (t in 1:T_max)
                log_cond_surv_02[j, t] += contrib;
            }
          }
```

NOTE: velocity at the grid edges uses `t-1=0` / `t+1=T_max+1`, i.e. a symmetric difference straddling the event window rather than clamped forward/backward. This is intentional and harmless — the analytic burden is defined for all `t` (no array bound), so the central formula is exact-in-form everywhere; only the dense path (bounded by `states_full_grid` columns) needs forward/backward at the array edges.

- [ ] **Step 2: Replace the 0→3 feature 2/3 block**

In the 0→3 loop, replace lines 87–96 (the `if (n_time_varying_covar >= 2)` / `>= 3` blocks) with the same structure, substituting `03` for `02` and `log_cond_surv_03`:

```stan
          if (enable_ms_velocity_basis) {
            for (t in 1:T_max) {
              real g_prev = fmin(log_sum_exp(base_d - rate_d * (t - 1), base_g + rate_g * (t - 1)), 10.0);
              real g_next = fmin(log_sum_exp(base_d - rate_d * (t + 1), base_g + rate_g * (t + 1)), 10.0);
              real vel = (g_next - g_prev) / 2.0;
              log_cond_surv_03[j, t] += time_varying_coef_03[2]
                * standardize_velocity(vel, median_velocity_burden_obs, iqr_velocity_burden_obs);
            }
          } else {
            if (n_time_varying_covar >= 2) {
              real contrib = time_varying_coef_03[2] * patient_log_decrease_rate[j, 1];
              for (t in 1:T_max)
                log_cond_surv_03[j, t] += contrib;
            }
            if (n_time_varying_covar >= 3) {
              real contrib = time_varying_coef_03[3] * patient_log_growth_rate[j, 1];
              for (t in 1:T_max)
                log_cond_surv_03[j, t] += contrib;
            }
          }
```

- [ ] **Step 3: Extend the 0→1 visit-gated block to add velocity**

The 0→1 block (lines 117–125) currently adds ONLY feature 1 (level) at visit weeks — it never read features 2/3. In velocity mode we add the velocity contribution at each visit week. Replace lines 117–126 (the `for (v in v_start:v_end)` loop body) with:

```stan
        for (v in v_start:v_end) {
          int wk = t_patient_visits[v];
          if (wk >= 1 && wk <= max_all_t) {
            real s_d = base_d - rate_d * wk;
            real s_g = base_g + rate_g * wk;
            real burden_covar = standardize_log_burden(
              log_sum_exp(s_d, s_g), log_bburden, median_log_burden_obs, iqr_log_burden_obs);
            log_cond_surv_01[j, wk] += time_varying_coef_01[1] * burden_covar;

            if (enable_ms_velocity_basis) {
              // Velocity at visit week wk: central difference of capped analytic
              // log-burden over a unit (weekly) grid straddling wk.
              real g_prev = fmin(log_sum_exp(base_d - rate_d * (wk - 1), base_g + rate_g * (wk - 1)), 10.0);
              real g_next = fmin(log_sum_exp(base_d - rate_d * (wk + 1), base_g + rate_g * (wk + 1)), 10.0);
              real vel = (g_next - g_prev) / 2.0;
              log_cond_surv_01[j, wk] += time_varying_coef_01[2]
                * standardize_velocity(vel, median_velocity_burden_obs, iqr_velocity_burden_obs);
            }
          }
        }
```

NOTE: velocity is a per-week instantaneous rate evaluated AT the visit week `wk`, using a ±1-week straddle — independent of the (irregular) gap to the next visit. This matches the per-week observed-delta scale (Task 3) and the dense-path weekly central difference (Task 5).

- [ ] **Step 4: Syntax-check the inline path via a model that uses it**

The publication model (`sf-ssm-log-space.stan`) compiles the inline path. Run:
`~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan 2>&1 | head -20`
Expected: clean parse.

- [ ] **Step 5: Also syntax-check the LFO variants that include the builder**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan 2>&1 | head -10
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo-endpoints.stan 2>&1 | head -10
```
Expected: clean parse for both.

- [ ] **Step 6: Confirm pioneer (legacy, flag off) still parses**

Run: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan 2>&1 | head -20`
Expected: clean parse. (Pioneer reads the same builder; the `else` branch is the unchanged legacy code. The velocity branch references `median_velocity_burden_obs`, which the PSA adapter does NOT alias — but that name only appears inside `if (enable_ms_velocity_basis)`, which pioneer compiles. **Stan compiles all branches regardless of flag value**, so the name must resolve. See Step 7.)

- [ ] **Step 7: Resolve the PSA missing-alias problem if Step 6 fails**

If Step 6 errors with `median_velocity_burden_obs` not found: the cleanest fix consistent with the spec (PSA stays legacy, no velocity constants) is to have the PSA adapter declare the aliases too, set to placeholder `(0, 1)` since they are never read when the flag is off. In `stan/psa/_psa_observed_covar_transformed_data.stan`, after the existing `iqr_log_burden_obs` alias, add:

```stan
// Placeholder velocity aliases so _ms_burden_*_tv_covar.stan parse cleanly for
// PSA. PSA stays on the legacy basis (enable_ms_velocity_basis = 0), so these
// are never read — the velocity branch is compiled but not executed.
real median_velocity_burden_obs = 0.0;
real iqr_velocity_burden_obs = 1.0;
```

Re-run Step 6; expect clean parse. (If Step 6 already passed, skip this step.)

- [ ] **Step 8: Commit**

```bash
git add stan/_ms_burden_inline_tv_covar.stan stan/psa/_psa_observed_covar_transformed_data.stan
git commit -m "feat(laplace): inline-path (level, velocity) branch (0->1/0->2/0->3)"
```

---

## Task 7: Config-error guard (matched flag + count)

**Files:**
- Modify: `stan/modules/multistate/transformed_data.stan` (append a guard near the top, after existing declarations)

A defensive `fatal_error` ensures `(enable_ms_velocity_basis, n_time_varying_covar)` is one of the two legal pairs, so a misconfigured R side fails loudly instead of silently mis-sizing.

- [ ] **Step 1: Find a safe insertion point**

Run: `grep -n "enable_ms_pop_time_varying_cov\|n_time_varying_covar\|^}" stan/modules/multistate/transformed_data.stan | head`
Expected: shows the transformed_data block structure. Insert the guard inside the block, after the opening declarations (any point where `enable_ms_velocity_basis` and `n_time_varying_covar` are in scope — both are data, so anywhere in the block works).

- [ ] **Step 2: Add the guard**

Near the top of the executable portion of `stan/modules/multistate/transformed_data.stan` (after variable declarations, before/independent of other logic), add:

```stan
// Guard: the (level, velocity) basis is a pure 2-vs-3 feature mode switch.
// enable_ms_velocity_basis must be paired with the matching feature count, or
// the coefficient arrays mis-size relative to what the builders emit.
if (enable_ms_pop_time_varying_cov) {
  if (enable_ms_velocity_basis && n_time_varying_covar != 2)
    fatal_error("enable_ms_velocity_basis=1 requires n_time_varying_covar=2 (level, velocity); got ",
                n_time_varying_covar);
  if (!enable_ms_velocity_basis && n_time_varying_covar == 2)
    fatal_error("n_time_varying_covar=2 with enable_ms_velocity_basis=0 is ambiguous; ",
                "set enable_ms_velocity_basis=1 for the (level, velocity) basis or use the 3-feature legacy basis.");
}
```

NOTE: the guard only fires under `enable_ms_pop_time_varying_cov` (when TV covariates are actually used). The visit-gated observed mode (`n_time_varying_covar` may be 1 internally) is unaffected because that path sets `enable_ms_velocity_basis = 0` and `n_tv_covar != 2` legacy. If the guard proves too strict in the standalone-prior path (`enable_ms_pop_time_varying_cov = 0`), it is already gated off there.

- [ ] **Step 3: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan 2>&1 | head -10`
Expected: clean parse.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/multistate/transformed_data.stan
git commit -m "feat(laplace): guard matched (velocity_basis, n_tv_covar) pair"
```

---

## Task 8: R config wiring — set the flag everywhere the model is built

**Files:**
- Modify: `targets/sclc_targets.R:1000` and `:2357`
- Modify: `targets/pioneer_targets.R:602`
- Modify: `targets/publication_targets.R` (`default_stan_data_settings` ~line 300, and the tribble ~line 230)

- [ ] **Step 1: sclc full model — flag OFF (legacy)**

In `targets/sclc_targets.R`, find line 1000 (`enable_ms_baseline_trend_01 = 0L,`). Immediately after it, add:

```r
          # (level, velocity) coupling basis: OFF for sclc (legacy
          # 3-feature [level, dec_rate, gro_rate] basis). Only the publication
          # model opts into the velocity basis (Laplace surrogate prerequisite).
          enable_ms_velocity_basis = 0L,
```

- [ ] **Step 2: sclc ms-standalone priors — flag OFF**

In `targets/sclc_targets.R`, find line 2357 (`enable_ms_baseline_trend_01 = 0L,` in the standalone-prior block). Immediately after it, add:

```r
            enable_ms_velocity_basis = 0L,
```

- [ ] **Step 3: pioneer — flag OFF**

In `targets/pioneer_targets.R`, find line 602 (`enable_ms_baseline_trend_01 = 0L,`). Immediately after it, add:

```r
        # (level, velocity) basis: OFF for pioneer (PSA stays on the legacy
        # 3-feature basis; velocity is not migrated for the PSA marker).
        enable_ms_velocity_basis = 0L,
```

- [ ] **Step 4: publication default_stan_data_settings — flag ON**

In `targets/publication_targets.R`, find line 300 (`enable_ms_baseline_trend_01 = enable_trend,`). Immediately after it, add:

```r
        # (level, velocity) coupling basis ON for the publication model. Matches
        # n_tv_covar = 2L in the tribble. Required for the Laplace surrogate: the
        # bi-exponential rates have no analog under the quadratic surrogate, but
        # velocity = d/dw of log-burden = b1 + 2*b2*w is linear in the
        # marginalized latents (keeps the survival term log-concave).
        enable_ms_velocity_basis = 1L,
```

- [ ] **Step 5: publication tribble — switch n_tv_covar 3L → 2L and cold-start**

In `targets/publication_targets.R`, find the tribble (~lines 230–243). Change:

```r
      n_tv_covar     = 3L,
```
to:
```r
      n_tv_covar     = 2L,   # (level, velocity) basis — see enable_ms_velocity_basis
```

Verify `warmstart = FALSE` is already set in that tribble (it is, per the surrogate cold-start note at `publication_targets.R:234-242`). If `warmstart` is `TRUE`, set it to `FALSE` — the 3→2 dimension change invalidates any saved inv-metric.

- [ ] **Step 6: Update the tribble comment block**

The comment at `publication_targets.R:213` says `"full" : n_time_varying_covar = 3 (standardised burden + log decrease rate + log growth rate)`. Update it to:

```r
  #   - "full"        : n_time_varying_covar = 2 (standardised burden LEVEL +
  #                     VELOCITY = central-difference of log-burden) under the
  #                     (level, velocity) basis (enable_ms_velocity_basis = 1L),
  #                     PLUS the hierarchical 0->1 baseline log-time trend.
```

- [ ] **Step 7: Verify the R files parse**

Run:
```bash
Rscript -e 'invisible(parse("targets/sclc_targets.R")); invisible(parse("targets/pioneer_targets.R")); invisible(parse("targets/publication_targets.R")); cat("parse OK\n")'
```
Expected: `parse OK`

- [ ] **Step 8: Commit**

```bash
git add targets/sclc_targets.R targets/pioneer_targets.R targets/publication_targets.R
git commit -m "feat(laplace): wire enable_ms_velocity_basis (publication ON, others OFF)"
```

---

## Task 9: Landing-gate validation — compile both modes + short MCMC

**Files:** none (validation only). Uses the targets pipeline.

- [ ] **Step 1: Run the unit tests (velocity helpers)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-velocity-covar.R")'`
Expected: PASS, `n_failures == 0`.

- [ ] **Step 2: Compile the publication model in VELOCITY mode (flag=1, n_tv=2)**

Build the publication stan-data + a model compile. The fastest check is to compile via cmdstanr with the publication data settings. Run the targets that assemble data and compile (NOT a full sample):

Run (MCP targets tool preferred, or):
```bash
Rscript -e 'library(targets); tar_make(names = matches("all_stan_data|default_stan_data_settings"), store = "/mnt/data/analysis-results/karim_naguib/publication/<RUN>/_targets")'
```
Then compile the model file directly:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: exit 0, no diagnostics. **Confirm the model BUILDS** (the `fatal_error` guard does not fire because data is (1, 2)).

NOTE: confirm the actual publication store path with `ls /mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/publication/` and the active `TAR_BRANCH`/run before running. Do not assume the run name.

- [ ] **Step 3: Compile in LEGACY mode (flag=0, n_tv=3) to confirm no regression**

Temporarily flip the publication tribble back to `n_tv_covar = 3L` and `enable_ms_velocity_basis = 0L` in a scratch check — OR rely on the sclc model (already flag=0, n_tv=3) compiling:

Run: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan` (sclc data settings)
Expected: clean parse in both modes. (The branch is compile-time identical; both `if`/`else` arms are always compiled. This step confirms the legacy arm is untouched.)

- [ ] **Step 4: Short MCMC in velocity mode — convergence + identification**

Launch a SHORT publication posterior fit in velocity mode (use the `pioneer-toolkit:start-job` skill / `/start-job`, or a local short run with reduced iterations). Target: a few hundred warmup + sampling iterations, 2 chains.

Acceptance:
- All `time_varying_coef_01[*]` and `time_varying_coef_03[*]` have `Rhat < 1.01`, adequate ESS.
- Zero divergences (or a count consistent with the legacy model's baseline).
- `time_varying_coef_01[2]` (velocity coef) posterior is **shifted off the prior** (not stuck at the `Normal(0, σ)` hyperprior mean) — confirms the feature is identified.

Run diagnostics via the `bayesian-toolkit:stan-diagnose` skill on the resulting fit.

- [ ] **Step 5: Scale sanity — standardized velocity ≈ unit IQR**

From the short fit (or a `generated_quantities` extraction), confirm the standardized velocity feature has IQR roughly of order 1 across patients/time. An IQR ≫ 1 or ≪ 1 signals a units bug (e.g. a missing per-week normalization in Task 3).

A quick check: extract `median_velocity_obs` / `iqr_velocity_obs` from the compiled model's transformed data (they are fixed) and confirm `iqr_velocity_obs > 0` and is on the scale of weekly log-SLD changes (typically 0.01–0.2 per week).

- [ ] **Step 6: Record the gate result**

Append a short PASS/FAIL note (with Rhat/ESS/divergence numbers and the velocity-coef posterior summary) to `docs/superpowers/specs/2026-06-09-laplace-surrogate-findings-and-blockers.md`.

- [ ] **Step 7: Commit the findings note**

```bash
git add docs/superpowers/specs/2026-06-09-laplace-surrogate-findings-and-blockers.md
git commit -m "docs(laplace): record (level, velocity) landing-gate validation result"
```

---

## Task 10: Documentation — re-justify the coupling as kinetics

**Files:**
- Modify: `quarto/sclc/website/documentation/multistate-specification.qmd`

- [ ] **Step 1: Locate the time-varying covariate description**

Run: `grep -n "time-varying\|decrease rate\|growth rate\|log SLD\|covariate" quarto/sclc/website/documentation/multistate-specification.qmd | head`
Expected: the section describing the 3-feature coupling.

- [ ] **Step 2: Add a paragraph describing the (level, velocity) basis**

Add a paragraph (respecting the repo's automatic section-numbering and cross-reference conventions — no hardcoded section numbers) explaining:
- The publication model couples the burden-driven transitions (0→1, 0→3) to the tumour trajectory via two time-varying features: the standardized log-burden **level** $g(w)$ and its **velocity** $g'(w)$ (the per-week rate of change, a central difference of the latent log-burden trajectory).
- This replaces the earlier (level, decrease-rate, growth-rate) basis: the bi-exponential rates are state-space component parameters with no analog under the Laplace surrogate's quadratic burden, whereas velocity is a derivative of the *output* trajectory, computable identically for forecast and backgrounded patients.
- Clinical reading: the hazard responds to where the burden is and which way it is moving (tumour kinetics).

- [ ] **Step 3: Verify the qmd has no obvious syntax breakage**

Run: `grep -c "velocity" quarto/sclc/website/documentation/multistate-specification.qmd`
Expected: ≥ 1 (the new paragraph is present). (Full `quarto render` is out of scope for this gate; the comparative-refit follow-up will re-render.)

- [ ] **Step 4: Commit**

```bash
git add quarto/sclc/website/documentation/multistate-specification.qmd
git commit -m "docs(laplace): document (level, velocity) kinetics coupling in MS spec"
```

---

## Task 11 (follow-up, NON-BLOCKING): Comparative refit

**Files:** none (analysis). Tracked here so it is not forgotten; does NOT block the landing of Tasks 1–10.

- [ ] **Step 1: Run both publication fits**

Run a full publication posterior fit in velocity mode (Tasks 1–8) and a second in legacy mode (tribble flipped to `n_tv_covar = 3L`, `enable_ms_velocity_basis = 0L`). Use the `pioneer-toolkit:start-job` skill for both.

- [ ] **Step 2: Compare PFS/OS KM fit quality**

Compare the two fits' predicted PFS/OS KM against observed (use the existing `km_est` comparison machinery — see `.claude/rules/sclc.md` for the `km_est` convention). 

**Acceptance:** the velocity basis must fit PFS/OS *at least as well* as the 3-feature basis before it is declared the publication default.

- [ ] **Step 3: Write up the result**

Document the comparison (KM overlays, coefficient interpretation, divergence/ESS comparison) in the findings doc and decide whether the velocity basis becomes the default.

---

## Self-Review Notes

- **Spec coverage:** §2 flag → Task 1; §3.1 velocity definition → Tasks 2,5,6; §3.3 standardization → Tasks 3,4; §4 coefficient machinery (no edit needed) → documented in "codebase facts"; §5 R wiring → Task 8; §2 config-error guard → Task 7; §6 validation (gate) → Task 9; §6.4 comparative refit → Task 11; clinical doc → Task 10. All spec sections map to a task.
- **Velocity unit consistency:** dense path differences `log_burden_normalized` (per-week, baseline cancels); inline path differences the capped analytic `log_sum_exp` (per-week ±1 straddle); Task 3 observed constants are per-week deltas. All three are per-week log-burden velocities — consistent scale. Verified in Step notes of Tasks 3, 5, 6.
- **Alias name consistency:** tumour-side constants are `median_velocity_obs`/`iqr_velocity_obs` (Task 3); burden-agnostic aliases read by the builders are `median_velocity_burden_obs`/`iqr_velocity_burden_obs` (Task 4, referenced in Tasks 5, 6, 6.7). Consistent across tasks.
- **Compile-all-branches caveat:** Stan compiles both `if`/`else` arms regardless of flag value, so the velocity branch's identifiers must resolve even for legacy/PSA models — handled by Task 6 Step 7 (PSA placeholder aliases). This is the subtlest correctness point.
```
