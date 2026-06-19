# Static Initial-Fraction Compartment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a third *static* (zero-growth) compartment to the tumor state-space model's initial-state split, so that stable-disease patients can be represented as flat rather than being structurally forced into an exponentially-growing compartment.

**Architecture:** The static compartment is constant in log-space (rate ≡ 0), so it is NOT a propagating state column. It is added as a single per-patient scalar `static_log_level` (= `log(π_static)`) inside the one combine function `calc_log_burden_mean`, where all SLD values are assembled via `log_sum_exp`. The baseline split becomes a nested conditional logit: `init_logit_loc` keeps its exact current meaning (decrease vs. rest), and a new population scalar `init_logit_static_loc_pop` splits the remainder into static vs. growth. With the feature flag off (default), `static_log_level = negative_infinity()` reproduces the current 2-component model exactly because `exp(-∞) = 0`.

**Tech Stack:** Stan (CmdStanR 2.38), modular `#include` architecture, R (tidyverse), testthat with `fixed_param` Stan-function unit tests.

**Spec:** `docs/superpowers/specs/2026-06-12-static-init-fraction-design.md`

**Key design facts (verified against current code):**
- Single combine funnel: `stan/modules/state_space/sf.stanfunctions:187` `calc_log_burden_mean(matrix, real)`.
- Call sites of `calc_log_burden_mean`: `sf.stanfunctions:870, 889, 1009, 1013, 1128, 1135, 1142`; `stan/tumor/sf-ssls-lfo.stan:213, 220, 227`.
- PSA uses a SEPARATE function `calc_log_psa_mean` (`stan/modules/psa/psa.stanfunctions:24`) — **out of scope**, untouched.
- `init` module 7 files at `stan/modules/init/`. New parameter follows the documented 8-step "Adding Module Parameters" checklist (`.claude/rules/stan-guidelines.md`).
- R wiring: `r/priors.R:243-244,347-348` (priors), `r/initializers.R:197-220` (init), `targets/publication_targets.R` (flag).
- Stan-function unit-test harness: `tests/testthat/helper-stan.R::test_stan_function(stan_file, data)` compiles a `tests/testthat/stan/<name>_all.stan` fixture with `include_paths = stan/` and runs 1 `fixed_param` draw. Driver lives in `tests/testthat/test-stan-<name>.R`.

---

## File Structure

**Stan — combine function (the math):**
- Modify `stan/modules/state_space/sf.stanfunctions` — add 3-arg `calc_log_burden_mean` overload; keep 2-arg as a thin wrapper; thread `static_log_level` through the two RNG functions and their call sites.

**Stan — init module (the parameter):**
- Modify `stan/modules/init/flags.stan` — `enable_static_init`
- Modify `stan/modules/init/hyperparams.stan` — prior mean/sd for the static logit
- Modify `stan/modules/init/parameters.stan` — `init_logit_static_loc_pop` (size-gated)
- Modify `stan/modules/init/transformed_parameters.stan` — nested-logit → `init_log_static_patient`
- Modify `stan/modules/init/priors.stan` — prior on the static logit
- Modify `stan/modules/init/generated_quantities.stan` — expose population π fractions

**Stan — wiring the scalar from transformed params into gen-quants/likelihood:**
- Modify `stan/modules/state_space/generated_quantities.stan` — pass `init_log_static_patient` to the all-patients RNG
- Modify `stan/tumor/sf-ssm-log-space.stan` — pass static level at the likelihood/mean call sites
- Modify `stan/tumor/sf-ssls-lfo.stan` — same for the LFO model

**R wiring:**
- Modify `r/priors.R` — default `init_logit_static_loc_pop_mean/sd`
- Modify `r/initializers.R` — init value for `init_logit_static_loc_pop`
- Modify `targets/publication_targets.R` — `enable_static_init` flag in stan-data assembly

**Tests:**
- Create `tests/testthat/stan/test_calc_log_burden_mean_all.stan` — fixture exercising 2-arg, 3-arg with `-inf`, and 3-arg with finite static
- Create `tests/testthat/test-stan-calc-log-burden-mean.R` — driver with expected values
- Create `tests/testthat/test-static-init-fraction.R` — R-level test that the nested-logit simplex sums to 1 and recovers the 2-way split when the flag is off

**Docs:**
- Modify `quarto/sclc/website/documentation/tumor-dynamics-specification.qmd` — document the static compartment (spec-as-contract rule)

---

## Task 1: Backward-compatible 3-arg `calc_log_burden_mean` overload

The heart of the change. We add a 3-arg overload taking `static_log_level`, and make the existing 2-arg signature delegate to it with `negative_infinity()`. Because `log_sum_exp(x, -∞) = x`, all existing callers are unaffected at the numeric level.

**Files:**
- Modify: `stan/modules/state_space/sf.stanfunctions:187-190`
- Create: `tests/testthat/stan/test_calc_log_burden_mean_all.stan`
- Create: `tests/testthat/test-stan-calc-log-burden-mean.R`

- [ ] **Step 1: Write the failing test fixture (Stan)**

Create `tests/testthat/stan/test_calc_log_burden_mean_all.stan`:

```stan
functions {
  #include "modules/state_space/sf.stanfunctions"
}

data {
  int<lower=1> n_visits;
  matrix[n_visits, 2] patient_states;  // log-space states (decrease, growth)
  real baseline;
  real static_log_level;               // log(pi_static); -inf disables
}

generated quantities {
  // Two-arg overload (must equal current behavior)
  vector[n_visits] out_2arg = calc_log_burden_mean(patient_states, baseline);
  // Three-arg overload with the supplied static level
  vector[n_visits] out_3arg = calc_log_burden_mean(patient_states, baseline, static_log_level);
}
```

- [ ] **Step 2: Write the failing test driver (R)**

Create `tests/testthat/test-stan-calc-log-burden-mean.R`:

```r
library(testthat)
library(here)

# Hand-computed expectations:
# state row = (log_decrease, log_growth); baseline multiplies inside log.
# 2-arg:  log(baseline) + log(exp(ld) + exp(lg))
# 3-arg:  log(baseline) + log(exp(ld) + exp(lg) + exp(static_log_level))
test_that("calc_log_burden_mean: 3-arg with -inf static equals 2-arg", {
  ld <- c(log(40), log(20))
  lg <- c(log(10), log(50))
  baseline <- 100
  data <- list(
    n_visits = 2L,
    patient_states = cbind(ld, lg),
    baseline = baseline,
    static_log_level = -Inf
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_calc_log_burden_mean_all.stan"),
    data
  )
  out_2 <- as.numeric(fit$draws("out_2arg", format = "draws_matrix"))
  out_3 <- as.numeric(fit$draws("out_3arg", format = "draws_matrix"))
  expect_equal(out_3, out_2, tolerance = 1e-10)
})

test_that("calc_log_burden_mean: finite static adds a constant compartment", {
  ld <- c(log(40), log(20))
  lg <- c(log(10), log(50))
  baseline <- 100
  static_log_level <- log(30)  # pi_static-scaled constant compartment
  data <- list(
    n_visits = 2L,
    patient_states = cbind(ld, lg),
    baseline = baseline,
    static_log_level = static_log_level
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_calc_log_burden_mean_all.stan"),
    data
  )
  out_3 <- as.numeric(fit$draws("out_3arg", format = "draws_matrix"))
  expected <- log(baseline) + log(exp(ld) + exp(lg) + exp(static_log_level))
  expect_equal(out_3, expected, tolerance = 1e-10)
})
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-calc-log-burden-mean.R")'`
Expected: FAIL — compilation error, no 3-arg overload of `calc_log_burden_mean` exists.

- [ ] **Step 4: Implement the overload**

In `stan/modules/state_space/sf.stanfunctions`, replace the existing function (lines 187-190):

```stan
vector calc_log_burden_mean(matrix patient_states, real baseline) {
  assert_equal(cols(patient_states), 2);
  return to_vector(log_sum_exp(patient_states[, 1], patient_states[, 2])) + log(baseline);
}
```

with the 3-arg overload plus a 2-arg wrapper:

```stan
/**
 * Three-component combine: decrease + growth (state columns) + static (scalar).
 * static_log_level = log(pi_static) in the same log-proportion units as the
 * state columns (before log(baseline) is added). Pass negative_infinity() to
 * disable the static compartment; exp(-inf)=0 recovers the 2-component model.
 */
vector calc_log_burden_mean(matrix patient_states, real baseline, real static_log_level) {
  assert_equal(cols(patient_states), 2);
  vector[rows(patient_states)] lse2 =
    to_vector(log_sum_exp(patient_states[, 1], patient_states[, 2]));
  return log_sum_exp(lse2, rep_vector(static_log_level, rows(patient_states))) + log(baseline);
}

// Backward-compatible 2-arg form: no static compartment.
vector calc_log_burden_mean(matrix patient_states, real baseline) {
  return calc_log_burden_mean(patient_states, baseline, negative_infinity());
}
```

Note: `log_sum_exp(vector, vector)` is elementwise in Stan and returns a vector — confirm against `stanc` in Step 5; if the elementwise vector overload is unavailable, use the explicit form:
`for (v in 1:rows(patient_states)) out[v] = log_sum_exp(lse2[v], static_log_level);`

- [ ] **Step 5: Syntax-check the model**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: no errors. If `log_sum_exp(vector, vector)` is rejected, switch to the explicit loop shown in Step 4 and re-run.

- [ ] **Step 6: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-calc-log-burden-mean.R")'`
Expected: PASS (both tests).

- [ ] **Step 7: Commit**

```bash
git add stan/modules/state_space/sf.stanfunctions \
        tests/testthat/stan/test_calc_log_burden_mean_all.stan \
        tests/testthat/test-stan-calc-log-burden-mean.R
git commit -m "feat(stan): add static compartment to calc_log_burden_mean

Three-arg overload adds a constant log-space compartment; 2-arg form
delegates with negative_infinity() so existing callers are byte-identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: `init` module — flag, parameter, prior, hyperparams

Add the new population scalar `init_logit_static_loc_pop` and its gating flag, following the documented 8-step checklist. The parameter is size-gated to length 0 when the flag is off so it costs nothing.

**Files:**
- Modify: `stan/modules/init/flags.stan:5`
- Modify: `stan/modules/init/hyperparams.stan:6`
- Modify: `stan/modules/init/parameters.stan:7`
- Modify: `stan/modules/init/priors.stan:6`

- [ ] **Step 1: Add the flag**

In `stan/modules/init/flags.stan`, after line 5 (`int<lower=0,upper=1> enable_pop_cov_init;`), add:

```stan
// Static compartment: when 1, baseline burden splits decrease/static/growth
// via a nested conditional logit (init_logit_loc unchanged: decrease vs. rest).
int<lower=0,upper=1> enable_static_init;
```

- [ ] **Step 2: Add the hyperparameters (prior mean/sd)**

In `stan/modules/init/hyperparams.stan`, after line 6 (`real<lower=0> init_logit_loc_pop_sd;`), add:

```stan
// Static-vs-growth split prior (used only when enable_static_init = 1).
// Always declared (cost-free real) so stan-data shape is flag-independent.
real init_logit_static_loc_pop_mean;
real<lower=0> init_logit_static_loc_pop_sd;
```

- [ ] **Step 3: Add the parameter (size-gated)**

In `stan/modules/init/parameters.stan`, after line 7 (`real init_logit_loc_pop;`), add:

```stan
// Static-vs-growth logit among the non-decreasing fraction.
// Length 1 when enabled, 0 when disabled (no sampling cost when off).
array[enable_static_init ? 1 : 0] real init_logit_static_loc_pop;
```

- [ ] **Step 4: Add the prior (guarded)**

In `stan/modules/init/priors.stan`, after line 6 (`init_logit_loc_pop ~ normal(...)`), add:

```stan
// Static-vs-growth split prior
if (enable_static_init) {
  init_logit_static_loc_pop[1] ~ normal(init_logit_static_loc_pop_mean, init_logit_static_loc_pop_sd);
}
```

- [ ] **Step 5: Syntax-check the model**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: no errors. (Compiles even though `init_log_static_patient` is not yet produced — that is Task 3.)

- [ ] **Step 6: Commit**

```bash
git add stan/modules/init/flags.stan stan/modules/init/hyperparams.stan \
        stan/modules/init/parameters.stan stan/modules/init/priors.stan
git commit -m "feat(stan/init): add enable_static_init flag and static logit parameter

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: `init` transformed parameters — nested-logit simplex

Produce `init_log_static_patient` (per-patient `log(π_static) + log(B₀)` is NOT added here — the baseline is applied inside `calc_log_burden_mean`; here we produce `log(π_static)` only, matching that `init_log_growth_patient` is also a log-proportion). Re-split the existing growth fraction so the three proportions form a simplex.

**Files:**
- Modify: `stan/modules/init/transformed_parameters.stan:102-103`

- [ ] **Step 1: Replace the 2-way split with the gated nested split**

In `stan/modules/init/transformed_parameters.stan`, replace lines 102-103:

```stan
vector[n_forecast_patients] init_log_decrease_patient = log_inv_logit(init_logit_loc_patient);
vector[n_forecast_patients] init_log_growth_patient   = log1m_inv_logit(init_logit_loc_patient);
```

with:

```stan
// Decrease keeps its exact meaning: log(pi_decrease) = log_inv_logit(loc).
vector[n_forecast_patients] init_log_decrease_patient = log_inv_logit(init_logit_loc_patient);

// "rest" = 1 - pi_decrease, in log space: log1m_inv_logit(loc).
// 2-way (flag off): all of "rest" is growth; static is absent (size 0).
// 3-way (flag on):  rest splits static/growth via init_logit_static_loc_pop.
//   log(pi_growth) = log(rest) + log1m_inv_logit(static)
//   log(pi_static) = log(rest) + log_inv_logit(static)
vector[n_forecast_patients] init_log_rest_patient = log1m_inv_logit(init_logit_loc_patient);
vector[enable_static_init ? n_forecast_patients : 0] init_log_static_patient;
vector[n_forecast_patients] init_log_growth_patient;
if (enable_static_init) {
  real ls = init_logit_static_loc_pop[1];
  init_log_static_patient = init_log_rest_patient + log_inv_logit(ls);
  init_log_growth_patient = init_log_rest_patient + log1m_inv_logit(ls);
} else {
  init_log_growth_patient = init_log_rest_patient;
}
```

- [ ] **Step 2: Syntax-check the model**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: no errors.

- [ ] **Step 3: Write an R test for the simplex math (failing first)**

Create `tests/testthat/test-static-init-fraction.R`:

```r
library(testthat)

# Mirror the Stan nested-logit math in R to lock the contract.
nested_split <- function(logit_loc, logit_static = NULL) {
  pi_decrease <- plogis(logit_loc)
  rest <- 1 - pi_decrease
  if (is.null(logit_static)) {
    return(c(decrease = pi_decrease, static = 0, growth = rest))
  }
  pi_static <- rest * plogis(logit_static)
  pi_growth <- rest * (1 - plogis(logit_static))
  c(decrease = pi_decrease, static = pi_static, growth = pi_growth)
}

test_that("2-way split (flag off) puts all non-decrease mass in growth", {
  p <- nested_split(0.3)
  expect_equal(unname(p["static"]), 0)
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that("3-way split is a valid simplex and recovers 2-way as static logit -> -Inf", {
  p3 <- nested_split(0.3, logit_static = 0.5)
  expect_equal(sum(p3), 1, tolerance = 1e-12)
  expect_true(all(p3 >= 0))
  p_recover <- nested_split(0.3, logit_static = -Inf)
  expect_equal(unname(p_recover["static"]), 0, tolerance = 1e-12)
  expect_equal(unname(p_recover["growth"]), unname(nested_split(0.3)["growth"]), tolerance = 1e-12)
})
```

- [ ] **Step 4: Run the R test**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-static-init-fraction.R")'`
Expected: PASS (this test documents the contract; it is pure R and should pass immediately — it guards future edits).

- [ ] **Step 5: Commit**

```bash
git add stan/modules/init/transformed_parameters.stan tests/testthat/test-static-init-fraction.R
git commit -m "feat(stan/init): nested-logit 3-way simplex for static compartment

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Thread `static_log_level` through the RNG functions

The two generator functions must accept and forward the static level. `generate_patient_states_rng` (single patient: scalar) and `generate_all_patients_states_with_means_rng` (all patients: vector indexed `[p]`). Both currently mirror `baseline_obs_value` / `baseline_obs_per_patient`.

**Files:**
- Modify: `stan/modules/state_space/sf.stanfunctions:848-892` (`generate_patient_states_rng`)
- Modify: `stan/modules/state_space/sf.stanfunctions:975-1014` (`generate_patient_states_with_means_rng`)
- Modify: `stan/modules/state_space/sf.stanfunctions:1061-1167` (`generate_all_patients_states_with_means_rng`)

- [ ] **Step 1: Add `static_log_level` scalar arg to `generate_patient_states_rng`**

In the signature of `generate_patient_states_rng` (around line 848-865), add a parameter `real static_log_level` immediately after `real baseline_obs_value`. At the two `calc_log_burden_mean` call sites inside it (`:870`, `:889`), add `, static_log_level` as the third argument:

```stan
rep_log_obs[2:] = to_vector(normal_rng(
  calc_log_burden_mean(patient_states[2:], baseline_obs_value, static_log_level), measure_sd));
...
forecast_log_obs = to_vector(normal_rng(
  calc_log_burden_mean(forecast_states[2:], baseline_obs_value, static_log_level), measure_sd));
```

- [ ] **Step 2: Forward through `generate_patient_states_with_means_rng`**

In `generate_patient_states_with_means_rng` (around line 975-1014), add `real static_log_level` after `real baseline_obs_value` in the signature. Forward it into the nested `generate_patient_states_rng(...)` call (around line 995-1005) and into the two `calc_log_burden_mean` calls (`:1009`, `:1013`):

```stan
generate_patient_states_rng(
  patient_states, forecast_time,
  patient_log_decrease_rate, patient_log_growth_rate,
  baseline_obs_value, static_log_level,
  forecast_growth_lag, forecast_growth_transition,
  forecast_process_noise, measure_sd);
...
vector[n_patient_visits] rep_mean_patient_log_obs =
  calc_log_burden_mean(patient_states, baseline_obs_value, static_log_level);
vector[forecast_size] forecast_mean_patient_log_obs =
  forecast_size > 0
    ? calc_log_burden_mean(forecast_patient_states, baseline_obs_value, static_log_level)
    : zeros_vector(0);
```

- [ ] **Step 3: Add `static_log_level_per_patient` vector to the all-patients function**

In `generate_all_patients_states_with_means_rng` (signature around line 1061), add `vector static_log_level_per_patient` immediately after `vector baseline_obs_per_patient`. Inside the per-patient loop, at the three `calc_log_burden_mean` call sites (`:1128`, `:1135`, `:1142`), pass `static_log_level_per_patient[p]`:

```stan
calc_log_burden_mean(patient_states[2:], baseline_obs_per_patient[p], static_log_level_per_patient[p])
...
calc_log_burden_mean(patient_states, baseline_obs_per_patient[p], static_log_level_per_patient[p])
...
calc_log_burden_mean(temp_forecast_patient_states, baseline_obs_per_patient[p], static_log_level_per_patient[p])
```

- [ ] **Step 4: Syntax-check (will fail at callers — expected)**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: errors at the *callers* of these functions (gen-quants + LFO) because their call signatures don't yet pass the new arg. That is fixed in Task 5. If there are errors *inside* `sf.stanfunctions` itself, fix those now.

- [ ] **Step 5: Commit**

```bash
git add stan/modules/state_space/sf.stanfunctions
git commit -m "feat(stan): thread static_log_level through state RNG generators

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Wire the static level at the model call sites (full + LFO)

Build a per-patient static-level vector in the model scope and pass it to the generators and direct `calc_log_burden_mean` calls. When the flag is off, the vector is all `negative_infinity()` (or we pass the scalar form), recovering current behavior.

**Files:**
- Modify: `stan/modules/state_space/generated_quantities.stan:14-89`
- Modify: `stan/tumor/sf-ssm-log-space.stan` (likelihood/mean call sites — locate via grep below)
- Modify: `stan/tumor/sf-ssls-lfo.stan:213, 220, 227`

- [ ] **Step 1: Construct the per-patient static-level vector in gen-quants**

In `stan/modules/state_space/generated_quantities.stan`, before the `profile("gen_quant_trajectories")` block (after line 12), add:

```stan
// Per-patient static log-level for the combine function.
// Flag off OR no static patients -> -inf (static compartment absent).
// init_log_static_patient is indexed forecast-locally (size n_forecast_patients);
// map to the unified patient index p via forecast_patient_idx, matching states_full_grid.
vector[n_patients] static_log_level_per_patient = rep_vector(negative_infinity(), n_patients);
if (enable_static_init) {
  for (j in 1:n_forecast_patients) {
    static_log_level_per_patient[forecast_patient_idx[j]] = init_log_static_patient[j];
  }
}
```

NOTE FOR IMPLEMENTER: verify `n_patients` is the correct length for `baseline_obs_per_patient` here (it is the array the all-patients RNG indexes by `[p]`). If `baseline_obs_per_patient` has a different length symbol, match it exactly.

- [ ] **Step 2: Pass the vector to the all-patients RNG**

In the `enable_patient_process_noise_tr` branch (around line 19-32), add `static_log_level_per_patient` immediately after `baseline_obs_per_patient` in the `generate_all_patients_states_with_means_rng(...)` argument list.

- [ ] **Step 3: Pass the scalar in the on-the-fly branch**

In the `else` branch's `generate_patient_states_with_means_rng(...)` call (around line 69-79), add the per-patient scalar after `baseline_obs_per_patient[p]`:

```stan
generate_patient_states_with_means_rng(
  states[state_start:state_end],
  forecast_time,
  patient_log_decrease_rate[j, 1],
  patient_log_growth_rate[j, 1],
  baseline_obs_per_patient[p],
  static_log_level_per_patient[p],   // NEW
  negative_infinity(),
  1.0,
  rep_matrix(0.0, forecast_size, 2),
  measure_sd_obs
);
```

- [ ] **Step 4: Locate and fix likelihood-side call sites in the full model**

Run: `grep -n "calc_log_burden_mean\|generate_patient_states\|generate_all_patients" stan/tumor/sf-ssm-log-space.stan`
For each call site, add the appropriate static argument: `static_log_level_per_patient[p]` (per-patient scalar) or the vector (all-patients), consistent with whether the call is inside a per-patient loop. If a likelihood mean uses `states` directly with `calc_log_burden_mean`, pass `static_log_level_per_patient[p]`.

- [ ] **Step 5: Fix the LFO model call sites**

In `stan/tumor/sf-ssls-lfo.stan`, at lines 213, 220, 227, add the static argument. Determine the in-scope per-patient static level the same way (build a `static_log_level_per_patient` vector in the LFO model's transformed-parameters/gen-quants scope mirroring Step 1). Each becomes:

```stan
calc_log_burden_mean(patient_states[2:], sum_tumor_size[visit_start], static_log_level_per_patient[p]),
...
calc_log_burden_mean(patient_states, sum_tumor_size[visit_start], static_log_level_per_patient[p]);
...
calc_log_burden_mean(forecast_patient_states, sum_tumor_size[visit_start], static_log_level_per_patient[p]);
```

NOTE FOR IMPLEMENTER: confirm the patient index symbol in scope at these LFO lines (it may be a local loop variable, not `p`); use whatever indexes `init_log_static_patient` correctly.

- [ ] **Step 6: Syntax-check BOTH models**

Run:
```
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: no errors for either.

- [ ] **Step 7: Commit**

```bash
git add stan/modules/state_space/generated_quantities.stan \
        stan/tumor/sf-ssm-log-space.stan stan/tumor/sf-ssls-lfo.stan
git commit -m "feat(stan): wire per-patient static level into full and LFO call sites

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Expose population fractions in generated quantities

For diagnostics and the recoverability check, output the three population proportions.

**Files:**
- Modify: `stan/modules/init/generated_quantities.stan`

- [ ] **Step 1: Add the population simplex outputs**

In `stan/modules/init/generated_quantities.stan`, after the existing `init_coef_pop` block, add:

```stan
// Population-level baseline composition (diagnostics).
// pi_decrease = inv_logit(loc_pop); rest = 1 - pi_decrease.
real init_pi_decrease_pop = inv_logit(init_logit_loc_pop);
real init_pi_static_pop = enable_static_init
  ? (1 - init_pi_decrease_pop) * inv_logit(init_logit_static_loc_pop[1])
  : 0.0;
real init_pi_growth_pop = (1 - init_pi_decrease_pop) - init_pi_static_pop;
```

- [ ] **Step 2: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: no errors.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/init/generated_quantities.stan
git commit -m "feat(stan/init): expose population decrease/static/growth fractions

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: R-side priors, initializers, and stan-data flag

Provide defaults so the model can be assembled and fit. Default the flag OFF so all existing pipelines are unchanged.

**Files:**
- Modify: `r/priors.R:243-244, 347-348`
- Modify: `r/initializers.R:197-220`
- Modify: `targets/publication_targets.R`

- [ ] **Step 1: Add prior defaults in `r/priors.R`**

After line 244 (`init_logit_loc_pop_sd <- 0.8 ...`), add:

```r
# Static-vs-growth split prior. Centered so static is a modest minority of the
# non-decreasing fraction at baseline; weakly informative.
init_logit_static_loc_pop_mean <- -0.85  # inv_logit(-0.85) ~ 0.30 of the "rest"
init_logit_static_loc_pop_sd <- 0.8
```

Then in the returned list (after line 348, `init_logit_loc_pop_sd = init_logit_loc_pop_sd,`), add:

```r
    init_logit_static_loc_pop_mean = init_logit_static_loc_pop_mean,
    init_logit_static_loc_pop_sd = init_logit_static_loc_pop_sd,
```

- [ ] **Step 2: Add the initializer in `r/initializers.R`**

In the block that sets `init_vals$init_logit_loc_pop` (around line 197-220), after `init_vals$init_logit_loc_pop <- init_logit_loc_pop`, add:

```r
    if (isTRUE(stan_data$enable_static_init == 1L)) {
      init_vals$init_logit_static_loc_pop <- as.array(rnorm(
        1,
        stan_data$init_logit_static_loc_pop_mean,
        stan_data$init_logit_static_loc_pop_sd
      ))
    }
```

NOTE FOR IMPLEMENTER: `as.array` because the Stan parameter is `array[1] real` (length-1), not a scalar; CmdStanR needs a 1-element array for it.

- [ ] **Step 3: Set the flag default in stan-data assembly**

In `targets/publication_targets.R`, locate where `enable_pop_cov_init` is set in the stan-data list (grep: `grep -n "enable_pop_cov_init\|enable_static_init" targets/publication_targets.R`). Add alongside it:

```r
    enable_static_init = 0L,  # static compartment OFF by default; flip to 1L per variant
```

If `enable_pop_cov_init` is not set in this file, set `enable_static_init = 0L` wherever the init flags are assembled for the publication stan data.

- [ ] **Step 4: Verify R sources without error**

Run: `Rscript -e 'source("r/priors.R"); source("r/initializers.R"); cat("OK\n")'`
Expected: `OK` (no parse/source errors).

- [ ] **Step 5: Run the existing testthat suite (backward-compat gate)**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: all pre-existing tests pass (flag defaults to 0, so nothing changes), plus the new tests from Tasks 1 and 3.

- [ ] **Step 6: Commit**

```bash
git add r/priors.R r/initializers.R targets/publication_targets.R
git commit -m "feat(r): default priors, initializer, and OFF flag for static init

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8: End-to-end compile + backward-compat fit smoke test

Confirm the full model compiles via CmdStanR (not just `stanc`) and that a tiny fit with the flag OFF behaves identically to baseline, and with the flag ON runs without error.

**Files:**
- Create: `tests/testthat/test-static-init-end-to-end.R` (skipped by default if no fixture data; otherwise minimal)

- [ ] **Step 1: Compile the full model through CmdStanR**

Run:
```r
Rscript -e 'library(cmdstanr); m <- cmdstan_model("stan/tumor/sf-ssm-log-space.stan", include_paths=c("stan","stan/tumor"), compile=TRUE); cat("compiled:", m$exe_file(), "\n")'
```
Expected: prints a compiled exe path, no errors.

- [ ] **Step 2: Confirm `enable_static_init` appears in the model's data**

Run:
```bash
grep -n "enable_static_init\|init_logit_static_loc_pop" stan/tumor/sf-ssm-log-space.hpp | head
```
Expected: the new symbols appear in the regenerated `.hpp` (confirms the include wired through). If the `.hpp` is stale, recompile in Step 1 first.

- [ ] **Step 3: Run the full testthat suite once more**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: green.

- [ ] **Step 4: Commit (if any fixture/test added)**

```bash
git add tests/testthat/test-static-init-end-to-end.R 2>/dev/null || true
git commit -m "test: end-to-end compile + backward-compat smoke for static init" --allow-empty
```

---

## Task 9: Update the model specification documentation

The spec-as-contract rule (`.claude/rules/r-guidelines.md`) requires the tumor-dynamics spec page to reflect Stan model changes.

**Files:**
- Modify: `quarto/sclc/website/documentation/tumor-dynamics-specification.qmd`

- [ ] **Step 1: Document the static compartment**

In `quarto/sclc/website/documentation/tumor-dynamics-specification.qmd`, find the section describing the two-component decrease/growth split (grep: `grep -n "growth\|decrease\|log_sum_exp\|init_logit_loc" quarto/sclc/website/documentation/tumor-dynamics-specification.qmd`). Add a subsection documenting:
- The nested conditional-logit 3-way split (decrease / static / growth), with the equations from the spec.
- That the static compartment is constant in log-space (rate ≡ 0) and enters via `calc_log_burden_mean`.
- That `enable_static_init = 0` recovers the 2-component model exactly.
- The motivation: relieving the structural forcing that makes stable-disease patients extrapolate to unbounded SLD.

- [ ] **Step 2: Commit**

```bash
git add quarto/sclc/website/documentation/tumor-dynamics-specification.qmd
git commit -m "docs(spec): document static initial-fraction compartment

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 10: Recoverability simulation (the gating validation)

The novel-structure gate from the spec. Simulate from a known 3-way split, fit, confirm the static fraction is recovered. This is exploratory analysis, not a unit test — it produces a short report, not a green/red assertion (though we add a loose posterior-coverage check).

**Files:**
- Create: `r/process_noise/recoverability_sim.R`

- [ ] **Step 1: Write the simulation script**

Create `r/process_noise/recoverability_sim.R`. It must:
1. Simulate `N` patients with a known population split (e.g. π_decrease ≈ 0.4, π_static ≈ 0.3, π_growth ≈ 0.3), realistic visit schedule, and measurement noise — reusing the model's own generative form (`exp(state_decrease) + π_static·B₀ + exp(state_growth)`).
2. Assemble stan-data with `enable_static_init = 1L` and the priors from `r/priors.R`.
3. Fit `stan/tumor/sf-ssm-log-space.stan` with a short run (e.g. 2 chains, 500 warmup / 500 sampling).
4. Extract `init_pi_static_pop`, `init_pi_decrease_pop`, `init_pi_growth_pop` and compare posterior to the true values.
5. Print a summary table (true vs. posterior mean/CI) and the divergence count.

The acceptance reading: the static fraction's posterior 90% CI should cover the true value and not collapse to 0, at the simulated follow-up length. Document the result inline in the script header as a comment after the first run.

NOTE FOR IMPLEMENTER: do NOT run any data-preparation pipeline; this is synthetic data generated in-script. Use `cmdstanr` directly, store path is irrelevant (no targets).

- [ ] **Step 2: Run the simulation**

Run: `Rscript r/process_noise/recoverability_sim.R`
Expected: prints the true-vs-recovered table; static fraction recovered (CI covers truth). If the static fraction is NOT recoverable, STOP and report — this is the gate that decides whether the approach is viable before any real fit.

- [ ] **Step 3: Commit**

```bash
git add r/process_noise/recoverability_sim.R
git commit -m "test(sim): recoverability check for static init fraction

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:**
- Spec §"Mathematical Parameterization" (nested conditional logits) → Task 3 ✓
- Spec §"Hierarchy scope (population intercept only)" → Task 2 (size-gated `array[1]` scalar) ✓
- Spec §A "Combine function (extra scalar argument)" → Task 1 ✓
- Spec §B "Threading the static scalar to call sites" → Tasks 4 + 5 ✓
- Spec §C "init module 8-step" → Tasks 2, 3, 6 (flags/hyperparams/parameters/transformed/priors/gen-quants) ✓
- Spec §D "R-side wiring" → Task 7 ✓
- Spec §E "Tests" → Tasks 1, 3, 8 (backward-compat via flag-off suite) ✓
- Spec §"Validation: recoverability sim" → Task 10 ✓; "backward-compat" → Task 7 Step 5 + Task 8; "explosion-reduction measurement" → deferred to post-implementation real-data run (out of plan scope, noted in spec) ; "syntax/compile gate" → every task's stanc step ✓

**Placeholder scan:** Three steps carry explicit "NOTE FOR IMPLEMENTER" verification asks (index symbol at LFO call sites; `n_patients` length symbol; `as.array` for length-1 param). These are real verification instructions, not vague placeholders — each names the exact thing to confirm and what to do. The `log_sum_exp(vector,vector)` uncertainty in Task 1 has an explicit fallback (the per-row loop). No "TBD"/"handle edge cases"/"similar to Task N".

**Type consistency:**
- `init_logit_static_loc_pop` is `array[enable_static_init ? 1 : 0] real` in parameters (Task 2), indexed `[1]` in priors (Task 2), transformed-params (Task 3), and gen-quants (Task 6) — consistent. R initializer uses `as.array(...)` length-1 (Task 7) — consistent.
- `init_log_static_patient` is `vector[enable_static_init ? n_forecast_patients : 0]` (Task 3), indexed forecast-locally `[j]` and mapped to `static_log_level_per_patient[p]` via `forecast_patient_idx` (Task 5) — consistent with how `states_full_grid` is indexed.
- `static_log_level` scalar arg name is identical across `calc_log_burden_mean`, `generate_patient_states_rng`, `generate_patient_states_with_means_rng` (Tasks 1, 4); `static_log_level_per_patient` vector arg in `generate_all_patients_states_with_means_rng` (Task 4) and the model scope (Task 5) — consistent.

**Open verification carried into execution (flagged, not hidden):** the exact patient-index symbol at `sf-ssls-lfo.stan:213/220/227` and the precise length symbol for `baseline_obs_per_patient` — both must be confirmed at implementation time and have explicit NOTE steps.
