# Gompertz Growth-Rate Decay Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the unbounded linear-in-time growth accumulation in the two-component log-space Stein-Fojo SSM with Gompertz growth-rate decay (`g(t)=growth_rate·e^{−κt}`), so long-horizon SLD forecasts plateau instead of exploding, gated by a master flag `enable_gr_decay` that defaults OFF (bit-identical to today's model).

**Architecture:** A new **lean pop-level** optional Stan module `stan/modules/gr_decay/` supplies a per-patient decay rate κ_i = exp(pop intercept + optional QR covariate slopes). A shared `growth_warp(t,κ)=φ(t)=(1−e^{−κt})/κ` helper warps growth-time at every growth-accumulation site (four branches in `state_space/transformed_parameters.stan` + two LFO forecast paths). All warp sites use φ-differences so constant-rate trajectories telescope to `growth_rate·φ(t)` identically across branches. Off-path takes a literal `t` branch — no `exp`, no helper call.

**Tech Stack:** Stan (CmdStanR 2.38, modular `#include`), R (tidyverse, testthat), targets pipeline.

**Spec:** `docs/superpowers/specs/2026-06-12-gompertz-growth-decay-design.md` (lean pop-level amendment, 2026-06-12).

---

## File Structure

**New Stan module — `stan/modules/gr_decay/` (6 files, NO `transformed_data.stan`):**
- `flags.stan` — `enable_gr_decay` (master) + `enable_pop_cov_gr_decay` (inner covariate gate)
- `hyperparams.stan` — `gr_decay_log_loc_pop_mean`/`_sd`, `gr_decay_coef_qr_pop_mean`/`_sd` (data, unconditional)
- `parameters.stan` — `gr_decay_log_loc_pop` (array size `enable_gr_decay?1:0`), `gr_decay_coef_qr_pop` (size `(enable_gr_decay&&enable_pop_cov_gr_decay)?n_covar:0`)
- `transformed_parameters.stan` — assemble `gr_decay_kappa` (size `n_forecast_patients`, zeros when off)
- `priors.stan` — `if (enable_gr_decay)` { intercept + guarded covariate priors }
- `generated_quantities.stan` — `gr_decay_kappa_pop` + implied plateau offset (diagnostics)

**New helper:** `growth_warp()` scalar + row_vector overloads in `stan/modules/state_space/sf.stanfunctions`.

**Modified Stan:**
- `stan/tumor/sf-ssm-log-space.stan` — `#include` the 6 gr_decay fragments
- `stan/tumor/sf-ssls-lfo.stan` — `#include` gr_decay fragments + warp the LFO forecast call (L205)
- `stan/modules/state_space/transformed_parameters.stan` — warp 4 branches (B1@74, B2@119, B3@145, B4@184)
- `stan/tumor/_lfo_endpoints_generated_quantities.stan` — warp the cutoff-compact forecast call (L201)

**Modified R:**
- `r/priors.R` — κ prior defaults inside `get_tumor_priors`'s returned `lst(...)`
- `r/sclc/initializers_fixed.R` — conditional `gr_decay_log_loc_pop` init
- `targets/publication_targets.R` — set `enable_gr_decay = 1L`, `enable_pop_cov_gr_decay = 0L`

**New tests:**
- `tests/testthat/stan/test_growth_warp.stan` + `tests/testthat/test-stan-growth-warp.R`
- `tests/testthat/stan/test_gr_decay_warp_branches.stan` + `tests/testthat/test-stan-gr-decay-branches.R`

**New deliverables (not commits):**
- `r/process_noise/recoverability_sim_gompertz.R` — adversarial κ-threshold sim
- prior-sensitivity narrative (run output captured in docs)

---

## Conventions every task must honour

- Native pipe `|>`, tidyverse, no `%>%`, no `install.packages()`, no hardcoded subject IDs.
- Stan: `expm1`, NOT deprecated `fabs` (use `abs`); built-in zero constructors; NCP.
- The off-path (`enable_gr_decay = 0`) is a **literal `? :` branch guard** at each growth site: `enable_gr_decay ? warp : linear`. It must produce bit-identical output to the current model.
- Compile syntax-check command (fast): `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
- Full test run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
- **Bootstrap on a fresh worktree** (empty renv): `R_PROFILE_USER=/dev/null Rscript --vanilla -e 'source("renv/activate.R"); renv::restore(prompt=FALSE)'`

---

## Task 1: `growth_warp()` helper (scalar + row_vector overloads)

**Files:**
- Modify: `stan/modules/state_space/sf.stanfunctions` (insert after `get_growth_lag_factor` overloads, i.e. after line 170, before `calc_log_burden_mean`)
- Create test: `tests/testthat/stan/test_growth_warp.stan`
- Create test: `tests/testthat/test-stan-growth-warp.R`

- [ ] **Step 1: Write the failing Stan test model**

Create `tests/testthat/stan/test_growth_warp.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
  #include "full_model.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "_burden.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
}

data {
  int<lower=1> n_t;
  vector[n_t] t;          // elapsed times (weeks)
  real kappa;             // decay rate
}

generated quantities {
  // Scalar warp at each time
  vector[n_t] phi;
  for (k in 1:n_t) phi[k] = growth_warp(t[k], kappa);

  // Row-vector overload must agree with scalar elementwise
  row_vector[n_t] phi_rv = growth_warp(to_row_vector(t), kappa);

  // Telescoping identity: sum of phi-differences from 0 == phi(t_n)
  // (build differences against a prepended 0)
  real phi_last = growth_warp(t[n_t], kappa);
  real telescoped = 0;
  {
    real prev = 0;  // phi(0) = 0
    for (k in 1:n_t) {
      real cur = growth_warp(t[k], kappa);
      telescoped += (cur - prev);
      prev = cur;
    }
  }
}
```

- [ ] **Step 2: Run the test model to verify it fails (helper undefined)**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor tests/testthat/stan/test_growth_warp.stan`
Expected: FAIL — `Identifier 'growth_warp' not in scope.`

- [ ] **Step 3: Implement the helper**

In `stan/modules/state_space/sf.stanfunctions`, immediately after the second `get_growth_lag_factor` overload (line 170), insert:

```stan
/**
 * Gompertz growth-time warp: phi(t) = integral_0^t exp(-kappa*s) ds = (1 - exp(-kappa*t)) / kappa.
 * Limits: kappa->0 => phi(t)->t (recovers linear growth exactly);
 *         t->inf   => phi(t)->1/kappa (growth arm plateaus).
 * The small-|kappa*t| branch (2-term Taylor) guards the 0/0 division AND keeps
 * d phi / d kappa finite for autodiff. expm1 fixes numerator precision only.
 */
real growth_warp(real t, real kappa) {
  if (abs(kappa * t) < 1e-8) {
    return t * (1 - 0.5 * kappa * t);   // phi ~= t - kappa*t^2/2
  }
  return -expm1(-kappa * t) / kappa;
}

// Elementwise row_vector overload (constant kappa) for the grid/cumsum branches.
row_vector growth_warp(row_vector t, real kappa) {
  int n = num_elements(t);
  row_vector[n] out;
  for (k in 1:n) out[k] = growth_warp(t[k], kappa);
  return out;
}
```

- [ ] **Step 4: Write the R test asserting the math**

Create `tests/testthat/test-stan-growth-warp.R`:

```r
library(testthat)
library(here)

phi_ref <- function(t, k) if (abs(k * t) < 1e-12) t else (1 - exp(-k * t)) / k

test_that("growth_warp matches phi and row_vector overload agrees with scalar", {
  data <- list(n_t = 5L, t = c(1, 4, 12, 52, 200), kappa = 0.05)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  phi <- vapply(1:5, function(k) get_stan_val(draws, "phi", k), numeric(1))
  phi_rv <- vapply(1:5, function(k) get_stan_val(draws, "phi_rv", k), numeric(1))
  expect_equal(phi, vapply(data$t, phi_ref, numeric(1), k = data$kappa), tolerance = 1e-6)
  expect_equal(phi, phi_rv, tolerance = 1e-10)
})

test_that("growth_warp recovers linear growth as kappa -> 0", {
  data <- list(n_t = 3L, t = c(1, 12, 52), kappa = 1e-10)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  phi <- vapply(1:3, function(k) get_stan_val(draws, "phi", k), numeric(1))
  expect_equal(phi, data$t, tolerance = 1e-6)
})

test_that("growth_warp plateaus at 1/kappa as t grows", {
  k <- 0.1
  data <- list(n_t = 1L, t = array(5000), kappa = k)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  expect_equal(get_stan_val(draws, "phi", 1), 1 / k, tolerance = 1e-4)
})

test_that("phi-differences telescope to phi(t_n)", {
  data <- list(n_t = 6L, t = c(2, 5, 9, 20, 60, 150), kappa = 0.03)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  expect_equal(
    get_stan_val(draws, "telescoped"),
    get_stan_val(draws, "phi_last"),
    tolerance = 1e-10
  )
})
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-growth-warp.R")'`
Expected: PASS (4 tests).

- [ ] **Step 6: Commit**

```bash
git add stan/modules/state_space/sf.stanfunctions tests/testthat/stan/test_growth_warp.stan tests/testthat/test-stan-growth-warp.R
git commit -m "feat(stan): add growth_warp() Gompertz time-warp helper with tests"
```

---

## Task 2: gr_decay module — flags + hyperparams (data layer)

**Files:**
- Create: `stan/modules/gr_decay/flags.stan`
- Create: `stan/modules/gr_decay/hyperparams.stan`

- [ ] **Step 1: Write flags.stan**

```stan
// gr_decay/flags.stan
// Gompertz growth-rate decay module — gating flags (optional module, propensity-style).
// enable_gr_decay is the MASTER gate on all gr_decay parameters/transforms/priors.

int<lower=0,upper=1> enable_gr_decay;          // master: warp growth-time when 1
int<lower=0,upper=1> enable_pop_cov_gr_decay;  // pop-level baseline-covariate slopes on log(kappa)
```

- [ ] **Step 2: Write hyperparams.stan**

```stan
// gr_decay/hyperparams.stan
// Hyperparameter (data) declarations for the Gompertz decay (gr_decay) module.
// Always declared with flag-independent shape, unused when off (cf. init/hyperparams.stan).

// Population intercept prior on log(kappa). Mean is NEGATIVE (e.g. log(0.02) ~ -3.9),
// so this is unconstrained real (NOT <lower=0>).
real gr_decay_log_loc_pop_mean;
real<lower=0> gr_decay_log_loc_pop_sd;

// QR-space covariate coefficient hyperparameters (applied in model block when enabled).
vector[n_covar] gr_decay_coef_qr_pop_mean;
vector<lower=0>[n_covar] gr_decay_coef_qr_pop_sd;
```

- [ ] **Step 3: Verify syntax in isolation is deferred to Task 8** (these fragments only compile inside the full model, which doesn't `#include` them until Task 8). No standalone gate here.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/gr_decay/flags.stan stan/modules/gr_decay/hyperparams.stan
git commit -m "feat(stan/gr_decay): add flags and hyperparams (data layer)"
```

---

## Task 3: gr_decay module — parameters

**Files:**
- Create: `stan/modules/gr_decay/parameters.stan`

- [ ] **Step 1: Write parameters.stan**

```stan
// gr_decay/parameters.stan
// Lean pop-level parameters for log(kappa). Master-flag gated: sized 0 when off
// (an off model samples/stores nothing), mirroring propensity/parameters.stan.

// Population intercept on log scale => kappa = exp(.) > 0. Array size 1 when on, 0 when off
// (mirrors init_logit_static_loc_pop).
array[enable_gr_decay ? 1 : 0] real gr_decay_log_loc_pop;

// Population covariate coefficients (QR space). Inner gate: 0 unless both the master
// flag and the covariate flag are on.
vector[(enable_gr_decay && enable_pop_cov_gr_decay) ? n_covar : 0] gr_decay_coef_qr_pop;
```

- [ ] **Step 2: Commit**

```bash
git add stan/modules/gr_decay/parameters.stan
git commit -m "feat(stan/gr_decay): add pop-level parameters (intercept + optional QR coefs)"
```

---

## Task 4: gr_decay module — transformed_parameters (per-patient κ)

**Files:**
- Create: `stan/modules/gr_decay/transformed_parameters.stan`

- [ ] **Step 1: Write transformed_parameters.stan**

```stan
// gr_decay/transformed_parameters.stan
// Per-patient decay rate kappa_i = exp(pop intercept + optional QR covariate slopes).
// Sizing: ALWAYS vector[n_forecast_patients], filled inside the flag guard and defaulting
// to zeros when off (the value is never consumed when off — warp sites take the t branch).
// Indexed by the forecast-local patient index j (matches frac_log_growth_patient, init_*).

vector[n_forecast_patients] gr_decay_log_loc_patient = zeros_vector(n_forecast_patients);
vector[n_forecast_patients] gr_decay_kappa = zeros_vector(n_forecast_patients);

if (enable_gr_decay) {
  vector[n_forecast_patients] gr_decay_linpred_pop = enable_pop_cov_gr_decay
    ? (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop)
    : zeros_vector(n_forecast_patients);
  gr_decay_log_loc_patient = rep_vector(gr_decay_log_loc_pop[1], n_forecast_patients)
    + gr_decay_linpred_pop;          // pop-only reduction: == gr_decay_log_loc_pop[1] when covariates off
  gr_decay_kappa = exp(gr_decay_log_loc_patient);
}
```

- [ ] **Step 2: Commit**

```bash
git add stan/modules/gr_decay/transformed_parameters.stan
git commit -m "feat(stan/gr_decay): compute per-patient kappa via pop-level linear predictor"
```

---

## Task 5: gr_decay module — priors + generated_quantities

**Files:**
- Create: `stan/modules/gr_decay/priors.stan`
- Create: `stan/modules/gr_decay/generated_quantities.stan`

- [ ] **Step 1: Write priors.stan**

```stan
// gr_decay/priors.stan — flag-gated weakly-informative prior on log(kappa).
// Entire block omitted when off (zero sampling cost).

if (enable_gr_decay) {
  gr_decay_log_loc_pop[1] ~ normal(gr_decay_log_loc_pop_mean, gr_decay_log_loc_pop_sd);
  if (enable_pop_cov_gr_decay) {
    gr_decay_coef_qr_pop ~ normal(gr_decay_coef_qr_pop_mean, gr_decay_coef_qr_pop_sd);
  }
}
```

- [ ] **Step 2: Write generated_quantities.stan**

```stan
// gr_decay/generated_quantities.stan — diagnostics for prior-vs-posterior contraction.
// pop_log_growth_rate is in scope (declared earlier in the model's GQ block).

real gr_decay_kappa_pop = enable_gr_decay ? exp(gr_decay_log_loc_pop[1]) : 0.0;
// Implied long-horizon plateau OFFSET on the log-growth arm: growth_rate / kappa.
real gr_decay_plateau_offset = enable_gr_decay
  ? exp(pop_log_growth_rate) / gr_decay_kappa_pop
  : positive_infinity();
```

- [ ] **Step 3: Commit**

```bash
git add stan/modules/gr_decay/priors.stan stan/modules/gr_decay/generated_quantities.stan
git commit -m "feat(stan/gr_decay): add prior and diagnostic generated quantities"
```

---

## Task 6: Warp the 4 growth-accumulation branches in state_space/transformed_parameters.stan

**Files:**
- Modify: `stan/modules/state_space/transformed_parameters.stan` (B1 @74, B2 @119–120, B3 @145, B4 @183–184)

> The state_space fragment is included AFTER frac/init in the model, and gr_decay's transformed_parameters must be included BEFORE state_space (wired in Task 8). `gr_decay_kappa` is therefore in scope here. `time_since_first_visit` is `row_vector[max_t_width]` with values `1..max_t_width`; *elapsed* time at grid position k is `time_since_first_visit[k] - 1`.

- [ ] **Step 1: Warp B1 — patient-level process noise (line 74)**

Replace line 74:
```stan
          states_full_grid[2][j, 2:] = to_row_vector(init_log_growth_patient[j] + cumulative_sum(patient_growth_rate[j, :(max_t_width-1)]));
```
with:
```stan
          if (enable_gr_decay) {
            // phi-differences over elapsed time so cumsum telescopes to growth_rate*phi(t).
            row_vector[max_t_width] warp = growth_warp(time_since_first_visit - 1, gr_decay_kappa[j]);
            row_vector[max_t_width-1] warp_diff = warp[2:] - warp[:(max_t_width-1)];
            states_full_grid[2][j, 2:] = to_row_vector(init_log_growth_patient[j]
              + cumulative_sum(patient_growth_rate[j, :(max_t_width-1)] .* warp_diff));
          } else {
            states_full_grid[2][j, 2:] = to_row_vector(init_log_growth_patient[j] + cumulative_sum(patient_growth_rate[j, :(max_t_width-1)]));
          }
```

- [ ] **Step 2: Warp B2 — pop-level process noise (lines 119–120), structural per-patient rewrite**

Replace lines 119–120:
```stan
        states_full_grid[2][, 2:] = rep_matrix(init_log_growth_patient, max_t_width - 1)
          + baseline_growth_rate * pop_cumsum_exp[:(max_t_width - 1)];
```
with:
```stan
        if (enable_gr_decay) {
          // Outer-product factorization breaks: the decay weight depends on BOTH patient
          // (via kappa_i) and time. Per-patient phi-difference cumsum (mirrors B1).
          for (i in 1:n_forecast_patients) {
            row_vector[max_t_width] warp = growth_warp(time_since_first_visit - 1, gr_decay_kappa[i]);
            row_vector[max_t_width-1] warp_diff = warp[2:] - warp[:(max_t_width-1)];
            states_full_grid[2][i, 2:] = init_log_growth_patient[i]
              + cumulative_sum(baseline_growth_rate[i] * (pop_exp[:(max_t_width-1)] .* warp_diff));
          }
        } else {
          states_full_grid[2][, 2:] = rep_matrix(init_log_growth_patient, max_t_width - 1)
            + baseline_growth_rate * pop_cumsum_exp[:(max_t_width - 1)];
        }
```

> Note: `pop_exp` (the per-step `exp(pop_noise)`, declared at L106) replaces `pop_cumsum_exp` because under the warp we re-accumulate per patient. With κ→0 the warp_diff → 1 and this reduces to `baseline_growth_rate[i] * cumulative_sum(pop_exp)` = the original outer-product row — confirm in the cross-branch test (Task 7).

- [ ] **Step 3: Warp B3 — no-noise full grid (line 145)**

Replace line 145:
```stan
      states_full_grid[2] = init_log_growth_patient * ones_row_vector(max_t_width) + patient_growth_rate[, 1] * (time_since_first_visit - 1);
```
with:
```stan
      if (enable_gr_decay) {
        // Per-patient warp of elapsed time (time_since_first_visit - 1), since kappa_i differs.
        for (i in 1:n_forecast_patients) {
          states_full_grid[2][i] = init_log_growth_patient[i]
            + patient_growth_rate[i, 1] * growth_warp(time_since_first_visit - 1, gr_decay_kappa[i]);
        }
      } else {
        states_full_grid[2] = init_log_growth_patient * ones_row_vector(max_t_width) + patient_growth_rate[, 1] * (time_since_first_visit - 1);
      }
```

- [ ] **Step 4: Warp B4 — no-noise direct visit (line 184)**

Replace line 184:
```stan
        states[state_start:state_end, 2] = init_log_growth_patient[j] + patient_growth_rate[j, 1] * dt;
```
with:
```stan
        if (enable_gr_decay) {
          vector[data_end - data_start + 1] dt_warp;
          for (v in 1:rows(dt)) dt_warp[v] = growth_warp(dt[v], gr_decay_kappa[j]);
          states[state_start:state_end, 2] = init_log_growth_patient[j] + patient_growth_rate[j, 1] * dt_warp;
        } else {
          states[state_start:state_end, 2] = init_log_growth_patient[j] + patient_growth_rate[j, 1] * dt;
        }
```

- [ ] **Step 5: Compile-gate the full model**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: clean parse (after Task 8 wires the includes). If run before Task 8, expect `gr_decay_kappa not in scope` — that is the signal to do Task 8 next. (Recommended execution order keeps Task 8 immediately after this one; the compile-gate is asserted at the end of Task 8.)

- [ ] **Step 6: Commit**

```bash
git add stan/modules/state_space/transformed_parameters.stan
git commit -m "feat(stan): warp all four growth-accumulation branches for Gompertz decay"
```

---

## Task 7: Cross-branch invariance test (constant-rate telescoping)

**Files:**
- Create: `tests/testthat/stan/test_gr_decay_warp_branches.stan`
- Create: `tests/testthat/test-stan-gr-decay-branches.R`

This is the invariance the naive Riemann sum would break: a constant-rate patient must get **identical** `growth_rate·φ(t)` whether computed by the cumsum form (B1/B2) or the direct form (B3/B4).

- [ ] **Step 1: Write the test Stan model (isolated, exercises both forms)**

Create `tests/testthat/stan/test_gr_decay_warp_branches.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
  #include "full_model.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "_burden.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
}

data {
  int<lower=1> n_t;        // grid width
  real growth_rate;        // constant rate
  real kappa;
  real init_g;
}

generated quantities {
  row_vector[n_t] elapsed = linspaced_row_vector(n_t, 0, n_t - 1);

  // Direct form (B3/B4): init + rate * phi(elapsed)
  row_vector[n_t] direct;
  for (k in 1:n_t) direct[k] = init_g + growth_rate * growth_warp(elapsed[k], kappa);

  // Cumsum form (B1/B2): init + cumsum(rate * (phi(t_k)-phi(t_{k-1})))
  row_vector[n_t] cumsum_form;
  cumsum_form[1] = init_g;
  {
    row_vector[n_t] warp = growth_warp(elapsed, kappa);
    row_vector[n_t-1] warp_diff = warp[2:] - warp[:(n_t-1)];
    cumsum_form[2:] = init_g + cumulative_sum(rep_row_vector(growth_rate, n_t-1) .* warp_diff);
  }
}
```

- [ ] **Step 2: Run to verify it fails (model uncompiled / helper just added)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-gr-decay-branches.R")'`
Expected: FAIL until Step 3's R file exists; if helper missing, compile error.

- [ ] **Step 3: Write the R test**

Create `tests/testthat/test-stan-gr-decay-branches.R`:

```r
library(testthat)
library(here)

test_that("cumsum (phi-difference) form equals direct form for constant rate", {
  for (k in c(1e-10, 0.01, 0.05, 0.2)) {
    data <- list(n_t = 40L, growth_rate = 0.08, kappa = k, init_g = log(3))
    fit <- test_stan_function(
      here::here("tests/testthat/stan/test_gr_decay_warp_branches.stan"), data
    )
    draws <- posterior::as_draws_df(fit$draws())
    direct <- vapply(1:40, function(i) get_stan_val(draws, "direct", i), numeric(1))
    cform  <- vapply(1:40, function(i) get_stan_val(draws, "cumsum_form", i), numeric(1))
    expect_equal(direct, cform, tolerance = 1e-9,
                 info = paste("kappa =", k))
  }
})

test_that("kappa -> 0 recovers the linear growth arm exactly", {
  data <- list(n_t = 20L, growth_rate = 0.08, kappa = 1e-12, init_g = log(3))
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_warp_branches.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  direct <- vapply(1:20, function(i) get_stan_val(draws, "direct", i), numeric(1))
  expected <- data$init_g + data$growth_rate * (0:19)
  expect_equal(direct, expected, tolerance = 1e-6)
})
```

- [ ] **Step 4: Run to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-gr-decay-branches.R")'`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/stan/test_gr_decay_warp_branches.stan tests/testthat/test-stan-gr-decay-branches.R
git commit -m "test(stan): cross-branch invariance for Gompertz phi-difference telescoping"
```

---

## Task 8: Wire gr_decay includes into both tumor models + warp both LFO forecast paths

**Files:**
- Modify: `stan/tumor/sf-ssm-log-space.stan` (data/parameters/transformed parameters/model/GQ includes)
- Modify: `stan/tumor/sf-ssls-lfo.stan` (same includes + warp the LFO forecast call at L205)
- Modify: `stan/tumor/_lfo_endpoints_generated_quantities.stan` (warp the cutoff-compact forecast call at L201)

- [ ] **Step 1: Add gr_decay includes to `sf-ssm-log-space.stan`**

In the `data` block, after `#include "modules/init/flags.stan"` (line 32):
```stan
  #include "modules/gr_decay/flags.stan"
  #include "modules/gr_decay/hyperparams.stan"
```
In the `parameters` block, after `#include "modules/init/parameters.stan"` (line 63):
```stan
  #include "modules/gr_decay/parameters.stan"
```
In the `transformed parameters` block, **before** `#include "modules/state_space/transformed_parameters.stan"` (i.e. after the init line 69):
```stan
  #include "modules/gr_decay/transformed_parameters.stan"
```
In the `model` block, after `#include "modules/init/priors.stan"` (line 86):
```stan
  #include "modules/gr_decay/priors.stan"
```
In the `generated quantities` block, after the `pop_log_growth_rate` declaration (line 129) so it is in scope:
```stan
  #include "modules/gr_decay/generated_quantities.stan"
```

- [ ] **Step 2: Mirror the same five includes into `sf-ssls-lfo.stan`**

Identical insertions at the matching locations (data after L31, parameters after L65, transformed parameters after L71 and before the state_space include L72, model after L87). For GQ: `sf-ssls-lfo.stan` has no `pop_log_growth_rate` in its GQ — so include only the prior/params/tp; **omit `gr_decay/generated_quantities.stan`** from the LFO model (its diagnostics belong to the full model). Confirm by reading the LFO GQ block.

- [ ] **Step 3: Warp the LFO forecast call in `sf-ssls-lfo.stan` (L205)**

The forecast loop variable is `i` (the original patient ID). κ for that patient is `gr_decay_kappa[forecast_local]`, but the LFO GQ does not carry a forecast-local index here. Mirror the `static_log_level_per_patient` pattern (L152–157): build a per-patient κ array indexed by the unified patient id.

After the `static_log_level_per_patient` block (L157), add:
```stan
  // Per-patient kappa for the forecast path (mirror of static_log_level_per_patient).
  // 0.0 sentinel => growth_warp(t,0)=t (no attenuation) when gr_decay is off.
  vector[n_patients] gr_decay_kappa_per_patient = zeros_vector(n_patients);
  if (enable_gr_decay) {
    for (j in 1:n_forecast_patients) {
      gr_decay_kappa_per_patient[forecast_patient_idx[j]] = gr_decay_kappa[j];
    }
  }
```

Then replace the forecast call (L205–214). **Key design point:** the decrease arm must NOT be warped — only the growth rate decays. So we cannot warp a shared time axis; instead we apply a per-step Gompertz factor to the *growth rate alone*, which is exactly what the helper's `time_varying_factor` already multiplies (`sf.stanfunctions:225,262`). To preserve the cross-branch φ-difference invariance, that factor is set to the **exact φ-difference per step** (`Δφ/Δt`), not a naive Riemann `e^{−κt}`. This needs one new helper overload `sf_log_space_trajectory_ncp_decay` (Step 4) that takes a precomputed factor vector instead of `(growth_lag, transition_rate)`.

Replace:
```stan
          (full_forecast_expected, full_forecast) = sf_log_space_trajectory_ncp(
            patient_states[visit_size],  // Last observed state as initial
            forecast_time,
            exp(patient_log_decrease_rate[i, 1]),
            exp(patient_log_growth_rate[i, 1]),
            negative_infinity(),  // growth lag disabled
            1.0,  // growth transition
            rep_matrix(0.0, n_oos_visits, 2),  // No process noise
            0  // No debug
          );
```
with:
```stan
          real kappa_i = gr_decay_kappa_per_patient[i];
          // Per-step Gompertz factor on the GROWTH rate only = exact phi-difference / dt,
          // so the forecast telescopes to growth_rate*phi(t) and matches the in-sample branches.
          vector[size(forecast_time)] tv_factor;
          for (t in 1:size(forecast_time)) {
            if (t == 1 || !enable_gr_decay) {
              tv_factor[t] = 1.0;
            } else {
              real e_hi = forecast_time[t] - forecast_time[1];
              real e_lo = forecast_time[t - 1] - forecast_time[1];
              real dphi = growth_warp(e_hi, kappa_i) - growth_warp(e_lo, kappa_i);
              real dt   = forecast_time[t] - forecast_time[t - 1];
              tv_factor[t] = dt > 0 ? dphi / dt : 1.0;
            }
          }
          (full_forecast_expected, full_forecast) = sf_log_space_trajectory_ncp_decay(
            patient_states[visit_size],
            forecast_time,
            exp(patient_log_decrease_rate[i, 1]),
            exp(patient_log_growth_rate[i, 1]),
            tv_factor,
            rep_matrix(0.0, n_oos_visits, 2),
            0
          );
```

> When `enable_gr_decay = 0`, `tv_factor` is all 1.0 → identical to the current disabled-lag call (backward-compat). With κ>0, the increment `tv_factor[t]·growth_rate·Δt = growth_rate·Δφ` is the **exact** φ-difference, so the forecast matches the in-sample warp (no Riemann gap).

- [ ] **Step 4: Add the `sf_log_space_trajectory_ncp_decay` helper**

In `stan/modules/state_space/sf.stanfunctions`, after the existing `sf_log_space_trajectory_ncp(... debug)` overload (line 270), add:
```stan
/**
 * Trajectory NCP variant taking a precomputed per-visit time_varying_factor (Gompertz
 * decay weight on the GROWTH rate only). decrease_rate is unscaled. With factor==1 this
 * is identical to the disabled-lag path; with factor=exp(-kappa*(t-t0)) the growth
 * increments telescope to growth_rate*(phi(t)-phi(t')).
 */
tuple(matrix, matrix) sf_log_space_trajectory_ncp_decay(
  row_vector x0, array[] real times,
  real decrease_rate, real growth_rate,
  vector time_varying_factor,
  matrix process_noise, int debug
) {
  int T = size(times);
  matrix[T, 2] x;
  matrix[T, 2] expected_x;
  x[1] = x0;
  expected_x[1] = x0;
  for (t in 2:T) {
    (expected_x[t], x[t]) = sf_log_space_transition_ncp(
      x[t - 1], times[t], times[t - 1],
      decrease_rate, time_varying_factor[t] * growth_rate, process_noise[t - 1]);
  }
  return (expected_x, x);
}
```

> The `tv_factor` built in Step 3 supplies the exact φ-difference per step, so this helper's `time_varying_factor[t] * growth_rate * Δt` accumulates `growth_rate · Δφ` — telescoping to `growth_rate · φ(t)`, identical to the in-sample branches. The plateau test (Task 11 Step 3) confirms the match.

- [ ] **Step 5: Apply the identical treatment to `_lfo_endpoints_generated_quantities.stan` (L201)**

That call uses `generate_patient_states_with_means_rng(...)`, which internally calls `sf_log_space_trajectory_ncp(... growth_lag, transition ...)`. Add a parallel decay-aware path: build `cutoff_gr_decay_kappa_per_patient` (mirroring `cutoff_static_log_level_per_patient` at L132–137) and pass an exact-φ-difference `tv_factor` into a new `generate_patient_states_with_means_decay_rng` overload, OR — simpler and lower-risk — replace the `forecast_time` axis passed to that function with the same exact-φ-difference factor applied only to the growth arm. Implement the same helper-overload approach as Step 4 for consistency. Read L120–223 of this file first and follow the `cutoff_static_*` pattern exactly.

- [ ] **Step 6: Compile-gate both models**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: both parse cleanly, no errors.

- [ ] **Step 7: Commit**

```bash
git add stan/tumor/sf-ssm-log-space.stan stan/tumor/sf-ssls-lfo.stan stan/tumor/_lfo_endpoints_generated_quantities.stan stan/modules/state_space/sf.stanfunctions
git commit -m "feat(stan): wire gr_decay module + warp both LFO forecast paths (exact phi-differences)"
```

---

## Task 9: R wiring — priors, initializer, target flags

**Files:**
- Modify: `r/priors.R` (inside the `lst(...)` returned by `get_tumor_priors`, ~line 350)
- Modify: `r/sclc/initializers_fixed.R` (inside `create_tumor_ssls_initializer_fixed`, ~line 333)
- Modify: `targets/publication_targets.R` (the default_stan_data_settings list, ~line 389)

- [ ] **Step 1: Add κ prior defaults to `get_tumor_priors`**

In `r/priors.R`, near the top of `get_tumor_priors` (after the other `*_pop_mean` scalars, ~line 248), add:
```r
  # Gompertz growth-rate decay: weakly-informative prior on log(kappa).
  # Centered at log(0.02) /week (growth-rate half-life ~35 wk; plateau over a multi-year
  # horizon). sd=0.75 => 95% prior kappa in ~[0.0045, 0.087], wide enough that data can
  # contract toward 0 (linear arms) or strong decay (bent arms). See spec Identifiability.
  gr_decay_log_loc_pop_mean <- log(0.02)
  gr_decay_log_loc_pop_sd <- 0.75
```
Then inside the returned `lst(...)` (after the `init_*` entries, ~line 359), add:
```r
    # Gompertz decay module hyperparams
    gr_decay_log_loc_pop_mean = gr_decay_log_loc_pop_mean,
    gr_decay_log_loc_pop_sd = gr_decay_log_loc_pop_sd,
    gr_decay_coef_qr_pop_mean = as.array(rep(0, n_covar)),
    gr_decay_coef_qr_pop_sd = as.array(rep(1, n_covar)),
```

- [ ] **Step 2: Verify the priors are reachable**

Run:
```bash
Rscript -e 'source("r/priors.R"); p <- get_tumor_priors(list(n_levels=2L, n_covar=0L), list(coef_mean=numeric(0), coef_sd=numeric(0))); stopifnot(p$gr_decay_log_loc_pop_mean == log(0.02), p$gr_decay_log_loc_pop_sd == 0.75); cat("OK\n")'
```
Expected: `OK`.

- [ ] **Step 3: Add the initializer**

In `r/sclc/initializers_fixed.R`, inside the `tumor_init <- tibble::lst(...)` block (alongside `tr_loc_pop = -2.0` etc., ~line 333), add a conditional entry. Since `lst()` evaluates expressions, gate on the stan_data flag (in scope via `with(stan_data, {...})`):
```r
        # Gompertz decay intercept: drawn from prior when enabled, length-0 when off.
        gr_decay_log_loc_pop = if (isTRUE(enable_gr_decay == 1L)) {
          as.array(rnorm(1, gr_decay_log_loc_pop_mean, gr_decay_log_loc_pop_sd))
        } else {
          numeric(0)
        },
        gr_decay_coef_qr_pop = if (isTRUE(enable_gr_decay == 1L) && isTRUE(enable_pop_cov_gr_decay == 1L) && n_covar > 0) {
          as.array(rep(0, n_covar))
        } else {
          numeric(0)
        },
```

> Confirm `gr_decay_log_loc_pop_mean`/`_sd` are in scope inside the initializer. If the initializer does not receive priors, hardcode the same numeric defaults (`log(0.02)`, `0.75`) with a comment cross-referencing `r/priors.R`. Read the function's argument list first.

- [ ] **Step 4: Set the flags in `publication_targets.R`**

In `targets/publication_targets.R`, in the `default_stan_data_settings` list, right after `enable_static_init = 0L,` (line 389), add:
```r
        enable_gr_decay = 1L,           # activate Gompertz growth-rate decay on this branch
        enable_pop_cov_gr_decay = 0L,   # pop-intercept only (no baseline-covariate slopes yet)
```
Leave `sclc_targets.R` and `pioneer_targets.R` unset so the data-prep/stan-data default (`0L`) flows there — confirm those pipelines default unset flags to 0 (read `prepare_tumor_stan_data`). If they do NOT default, add `enable_gr_decay = 0L` explicitly to each.

- [ ] **Step 5: Verify targets file sources cleanly**

Run: `Rscript -e 'source("targets/publication_targets.R")' 2>&1 | tail -5` (expect no parse error; target graph errors about missing stores are unrelated).

- [ ] **Step 6: Commit**

```bash
git add r/priors.R r/sclc/initializers_fixed.R targets/publication_targets.R
git commit -m "feat(r): wire gr_decay priors, initializer, and publication flags"
```

---

## Task 10: Backward-compat gate + R hierarchy-reduction test

**Files:**
- Create: `tests/testthat/test-gr-decay-pop-reduction.R`

- [ ] **Step 1: Backward-compat — full suite green with flag OFF**

The defaults leave `enable_gr_decay` unset/0 everywhere except `publication_targets.R`. Run the full suite:

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: PASS at the same level as base commit `378c8a5e` (the pre-existing `arm`-column / `conflicted::is_null` failures are unrelated — confirm they are identical on base, do not attempt to fix them here).

- [ ] **Step 2: Write the pop-only reduction test**

This asserts the spec's "hierarchy reduces correctly" claim: with `enable_pop_cov_gr_decay = 0`, every patient's `gr_decay_log_loc_patient` equals `gr_decay_log_loc_pop[1]`. Test against the gr_decay transformed_parameters logic in isolation via a tiny model.

Create `tests/testthat/stan/test_gr_decay_reduction.stan`:
```stan
data {
  int<lower=1> n_forecast_patients;
  int<lower=0> n_covar;
  int<lower=0,upper=1> enable_gr_decay;
  int<lower=0,upper=1> enable_pop_cov_gr_decay;
  array[n_forecast_patients] int forecast_patient_idx;
  matrix[n_forecast_patients, n_covar] Q_covar_design_matrix;
  array[enable_gr_decay ? 1 : 0] real gr_decay_log_loc_pop;
  vector[(enable_gr_decay && enable_pop_cov_gr_decay) ? n_covar : 0] gr_decay_coef_qr_pop;
}
generated quantities {
  vector[n_forecast_patients] gr_decay_log_loc_patient = zeros_vector(n_forecast_patients);
  if (enable_gr_decay) {
    vector[n_forecast_patients] linpred = enable_pop_cov_gr_decay
      ? (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop)
      : zeros_vector(n_forecast_patients);
    gr_decay_log_loc_patient = rep_vector(gr_decay_log_loc_pop[1], n_forecast_patients) + linpred;
  }
}
```

Create `tests/testthat/test-gr-decay-pop-reduction.R`:
```r
library(testthat)
library(here)

test_that("pop-only config: kappa constant across patients", {
  np <- 6L
  data <- list(
    n_forecast_patients = np, n_covar = 2L,
    enable_gr_decay = 1L, enable_pop_cov_gr_decay = 0L,
    forecast_patient_idx = 1:np,
    Q_covar_design_matrix = matrix(rnorm(np * 2), np, 2),
    gr_decay_log_loc_pop = array(log(0.02)),
    gr_decay_coef_qr_pop = numeric(0)
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_reduction.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  vals <- vapply(1:np, function(i) get_stan_val(draws, "gr_decay_log_loc_patient", i), numeric(1))
  expect_equal(vals, rep(log(0.02), np), tolerance = 1e-9)
})

test_that("flag off: kappa log-loc is zero (value unused downstream)", {
  np <- 4L
  data <- list(
    n_forecast_patients = np, n_covar = 0L,
    enable_gr_decay = 0L, enable_pop_cov_gr_decay = 0L,
    forecast_patient_idx = 1:np,
    Q_covar_design_matrix = matrix(0, np, 0),
    gr_decay_log_loc_pop = numeric(0),
    gr_decay_coef_qr_pop = numeric(0)
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_reduction.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  vals <- vapply(1:np, function(i) get_stan_val(draws, "gr_decay_log_loc_patient", i), numeric(1))
  expect_equal(vals, rep(0, np), tolerance = 1e-12)
})
```

- [ ] **Step 3: Run the new test**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-gr-decay-pop-reduction.R")'`
Expected: PASS (2 tests).

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/stan/test_gr_decay_reduction.stan tests/testthat/test-gr-decay-pop-reduction.R
git commit -m "test: gr_decay pop-only reduction + backward-compat gate"
```

---

## Task 11: End-to-end compile + short sampling smoke test

**Files:** none (verification only)

- [ ] **Step 1: Compile both tumor models via CmdStanR**

Run:
```bash
Rscript -e 'library(cmdstanr); m <- cmdstan_model("stan/tumor/sf-ssm-log-space.stan", include_paths="stan", compile=TRUE); cat("compiled OK\n")'
```
Expected: `compiled OK`.

- [ ] **Step 2: Short smoke fit with `enable_gr_decay = 1` on synthetic data**

Use the smallest available test stan_data fixture (or the recoverability sim's data generator from Task 12) with `enable_gr_decay = 1L`, `enable_pop_cov_gr_decay = 0L`, 1 chain, 50 warmup / 50 sampling. Assert: (a) `gr_decay_log_loc_pop` is in `fit$draws()`, (b) zero post-warmup exceptions, (c) `gr_decay_kappa_pop` finite and positive.

- [ ] **Step 3: Plateau assertion**

From the smoke fit, generate a long-horizon forecast and assert the mean log-growth arm approaches `init_log_growth + growth_rate / kappa` (within tolerance) as horizon → large, and does NOT exceed it. This is the LFO plateau check from the spec's Testing section.

- [ ] **Step 4: No commit** (verification gate). If any step fails, return to the relevant task.

---

## Task 12 (deliverable, not a code gate): Adversarial recoverability sim

**File:** Create `r/process_noise/recoverability_sim_gompertz.R`

Per the spec's "Recoverability simulation" section — this is **diagnostic**, characterizing the κ-detection threshold under empirical follow-up, NOT a pass/fail gate.

- [ ] **Step 1: Implement the sim**
  - Simulate N≈300 patients; draw follow-up length and nadir timing from the **empirical SCLC/sclc distribution** (most censored while declining or ≤~7wk post-nadir; minority ≥12wk). Cross-reference `[[lfo-forecast-regrowth-bias]]`.
  - Per-patient `growth_rate` heterogeneity; realistic measurement noise.
  - Sweep κ_true ∈ {0, 0.01, 0.02, 0.05, 0.1} (spanning prior center → clearly decelerating).
  - Assemble stan_data with `enable_gr_decay = 1L`, priors from `get_tumor_priors`; fit `sf-ssm-log-space.stan` short (2 chains, 500/500) per κ_true.
- [ ] **Step 2: Report per κ_true** — posterior mean/CI for `gr_decay_kappa_pop`, prior-vs-posterior contraction (does it move off the prior?), divergences.
- [ ] **Step 3: Capture the output table** into the prior-sensitivity narrative (Task 13). No git gate; this is run interactively.

---

## Task 13 (deliverable): Prior-sensitivity analysis + doc page

**Files:**
- Create/append: model-spec documentation (the tumor-dynamics specification page) — document Gompertz growth-rate decay, the warp, the weakly-informative prior, and the identifiability precondition.

- [ ] **Step 1: Prior-sensitivity** — vary `gr_decay_log_loc_pop_mean`/`_sd`; report (a) forecast-SLD plateau and (b) in-window fit (log-lik / LOO) at each setting. Make the prior's influence explicit (the spec requires this as the honest companion to a weakly-identified parameter).
- [ ] **Step 2: Update the tumor-dynamics specification page** (per r-guidelines: "When changing Stan code, update the relevant specification page"). Document the φ(t) warp, the four-branch wiring, the LFO plateau, and the `enable_gr_decay`/`enable_pop_cov_gr_decay` flags.
- [ ] **Step 3: Commit docs**

```bash
git add quarto/.../tumor-dynamics-specification.qmd docs/  # exact path per project
git commit -m "docs: document Gompertz growth-rate decay and prior-sensitivity"
```

---

## Execution Order & Notes

1. **T1** (helper + tests) — foundational, fully testable in isolation.
2. **T2–T5** (gr_decay module files) — independent file creations; compile only inside the model.
3. **T6** (warp 4 branches) then **T8** (wire includes) — T6's compile-gate is asserted at the end of T8 (the includes make `gr_decay_kappa` reachable).
4. **T7** (cross-branch invariance) — can run right after T1 (uses only the helper).
5. **T9** (R wiring) → **T10** (backward-compat + reduction tests) → **T11** (e2e compile/smoke/plateau).
6. **T12** (sim) + **T13** (prior-sensitivity + docs) — deliverables, no code gate.

**Critical correctness invariants (the make-or-break):**
- κ→0 ⟹ φ(t)→t at *every* site → flag-off bit-identical (asserted by backward-compat suite + reduction test).
- Constant-rate patient gets identical `growth_rate·φ(t)` across B1/B2/B3/B4 and both LFO paths (asserted by the cross-branch test; the LFO uses the exact φ-difference `tv_factor`, NOT a naive Riemann `e^{−κt}·Δt`).
- The **decrease** arm is never warped — only the growth rate is scaled by the Gompertz factor (this is why the LFO path uses `time_varying_factor` on `growth_rate` alone, not a warped shared time axis).
- `gr_decay_log_loc_pop_mean` is **negative** (`log(0.02)`) ⟹ hyperparam is unconstrained `real`, never `<lower=0>`.

**Pre-existing test failures** (unrelated, identical on base `378c8a5e`): `arm` column and `conflicted::is_null` errors. Do not fix here; just confirm parity.
