# Gompertz κ RE Hierarchy (trial-arm + patient) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the gr_decay module's pooled population-scalar κ to a full level-indexed RE hierarchy (population intercept + per-level NONE/FE/RE/RE_CP intercepts + per-level covariate slopes), mirroring the `frac` module exactly, so κ can vary by trial-arm and patient. Covariate-slope machinery is wired but defaults OFF.

**Architecture:** κ's linear predictor on log scale becomes `gr_decay_log_loc_pop + linpred_pop(QR cov) + Σ_lv level_intercepts[lv] + Σ_lv level_slopes[lv]`, exactly like `frac_logit_loc_patient`. Uses the shared `compute_level_module_flags` / `split_cp_ncp_pos` helpers (hierarchy.stanfunctions) for flat-index gathers. Downstream warp sites are UNCHANGED — they consume κ only via the forecast-local `vector[n_forecast_patients] gr_decay_kappa`, which this still produces.

**Tech Stack:** Stan (modular #include), CmdStanR, R (tidyverse), targets, testthat.

**Key invariants (must hold at every step):**
- `enable_gr_decay=0` ⟹ zero gr_decay params, bit-identical to pre-gr_decay model.
- `enable_level_intercept_gr_decay = c(none, none)` AND `enable_pop_cov_gr_decay=0` ⟹ reduces EXACTLY to the current pop-only scalar κ (the just-fitted job #1952 config). This is the backward-compat anchor.
- κ is consumed downstream ONLY via `gr_decay_kappa[j]`. No warp site (`sf.stanfunctions` overloads, `state_space/transformed_parameters.stan`, `sf-ssls-lfo.stan`, `_lfo_endpoints_generated_quantities.stan`) changes.
- Mirror `frac`, NOT `tr` (no process-noise, no SD sub-hierarchy, no nu/Student-t unless frac's loop already carries it — gr_decay does NOT add Student-t).

**Reference files (read before implementing):**
- `stan/modules/frac/{flags,hyperparams,parameters,transformed_data,transformed_parameters,priors}.stan` — the exact template.
- `stan/hierarchy.stanfunctions` — `compute_level_module_flags`, `split_cp_ncp_pos`.
- `stan/_hierarchy_transformed_data.stan` — LEVEL_MODE_* constants, `n_forecast_groups_per_level`.
- `r/priors.R:344-355` (frac block), `r/sclc/initializers_fixed.R:296-372` (frac raw-level construction).
- `targets/publication_targets.R:372-396` (frac/init level config).

---

## File Structure

- **Modify** `stan/modules/gr_decay/flags.stan` — add `enable_level_intercept_gr_decay[n_levels]`, `enable_level_cov_gr_decay[n_levels]`.
- **Modify** `stan/modules/gr_decay/hyperparams.stan` — add per-level SD hyperparams (intercept sd, fe sd, slope sd).
- **Create** `stan/modules/gr_decay/transformed_data.stan` — flat-index + bucket-position computation (copy frac, rename).
- **Modify** `stan/modules/gr_decay/parameters.stan` — add level SD raw, raw/cp intercept, slope SD, raw/cp slope params.
- **Modify** `stan/modules/gr_decay/transformed_parameters.stan` — assemble `gr_decay_log_loc_patient` from pop + level intercepts + slopes; `gr_decay_kappa = exp(...)`.
- **Modify** `stan/modules/gr_decay/priors.stan` — frac-style unified level loop.
- **Modify** `stan/tumor/sf-ssm-log-space.stan` AND `stan/tumor/sf-ssls-lfo.stan` — add `#include "modules/gr_decay/transformed_data.stan"` after the frac transformed_data include.
- **Modify** `r/priors.R` — add gr_decay per-level SD hyperparam defaults.
- **Modify** `r/sclc/initializers_fixed.R` — add gr_decay raw-level / cp / slope inits.
- **Modify** `targets/publication_targets.R` — set `enable_level_intercept_gr_decay = c(trial_arm=re, patient=re)`, `enable_level_cov_gr_decay = c(trial_arm=FALSE, patient=FALSE)`.
- **Create** `tests/testthat/stan/test_gr_decay_hierarchy.stan` + `tests/testthat/test-stan-gr-decay-hierarchy.R` — reduction-to-pop-scalar test.

---

### Task 1: Flags — add per-level intercept and cov flag arrays

**Files:**
- Modify: `stan/modules/gr_decay/flags.stan`

- [ ] **Step 1: Add the flag arrays** (mirror `frac/flags.stan` lines 11-12)

Append to `stan/modules/gr_decay/flags.stan` after the existing two flags:

```stan
// Per-level intercept modes on log(kappa): index 1..n_levels.
// 0=NONE, 1=FE, 2=RE, 3=RE_GP (rejected here), 4=RE_CP. Mirrors frac.
array[n_levels] int<lower=0,upper=4> enable_level_intercept_gr_decay;
// Per-level covariate-slope flags on log(kappa). Default all 0 (off).
array[n_levels] int<lower=0,upper=1> enable_level_cov_gr_decay;
```

- [ ] **Step 2: Syntax check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: FAIL — because parameters.stan/transformed_data.stan reference `n_re_levels_gr_decay_intercept` etc. that don't exist yet. (This is the failing-first signal; the flags themselves parse.) If the ONLY errors are about not-yet-defined gr_decay hierarchy symbols, proceed. If there's a syntax error in flags.stan itself, fix it.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/gr_decay/flags.stan
git commit -m "feat(stan/gr_decay): per-level intercept + cov flag arrays

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Hyperparams — per-level SD hyperparameters

**Files:**
- Modify: `stan/modules/gr_decay/hyperparams.stan`

- [ ] **Step 1: Add per-level SD hyperparams** (mirror `frac/hyperparams.stan`)

Append to `stan/modules/gr_decay/hyperparams.stan` (keep the existing pop-intercept-mean/sd and pop-cov-mean/sd lines):

```stan
// Hierarchical intercept prior scale hyperparameters — one per level (RE/RE_CP).
array[n_levels] real<lower=0> gr_decay_sd_level_intercept_sd;

// Fixed-effect SD for FE levels (mode=1); unused for RE/disabled levels.
array[n_levels] real<lower=0> gr_decay_fe_sd_level_intercept;

// Hierarchical slope SD hyperpriors — one row_vector per level.
array[n_levels] row_vector<lower=0>[n_covar] gr_decay_sd_level_slope_sd;
```

- [ ] **Step 2: Syntax check** (same command as Task 1 Step 2)

Expected: still FAIL on not-yet-defined parameters; hyperparams parse OK.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/gr_decay/hyperparams.stan
git commit -m "feat(stan/gr_decay): per-level SD hyperparameters

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: transformed_data — flat indices + bucket positions

**Files:**
- Create: `stan/modules/gr_decay/transformed_data.stan`

- [ ] **Step 1: Create the file** (copy `frac/transformed_data.stan`, replace `frac`→`gr_decay`, `enable_level_cov_frac`→`enable_level_cov_gr_decay`, `enable_level_intercept_frac`→`enable_level_intercept_gr_decay`)

Create `stan/modules/gr_decay/transformed_data.stan`:

```stan
// gr_decay/transformed_data.stan
// Flat indices + bucket positions for the level-indexed log(kappa) hierarchy.
// Mirrors frac/transformed_data.stan. Computed unconditionally (cheap; all counts
// collapse to 0 when every level mode is NONE, so the pop-only path is bit-exact).

int n_enabled_groups_gr_decay_intercept;
array[n_levels + 1] int enabled_level_pos_gr_decay_intercept;
array[n_patients, n_levels] int patient_gr_decay_intercept_flat_idx;
int n_enabled_groups_gr_decay_slope;
array[n_levels + 1] int enabled_level_pos_gr_decay_slope;
array[n_patients, n_levels] int patient_gr_decay_slope_flat_idx;
(n_enabled_groups_gr_decay_intercept, enabled_level_pos_gr_decay_intercept, patient_gr_decay_intercept_flat_idx,
 n_enabled_groups_gr_decay_slope, enabled_level_pos_gr_decay_slope, patient_gr_decay_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_gr_decay, enable_level_cov_gr_decay);

int n_re_levels_gr_decay_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE ||
      enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) n_re_levels_gr_decay_intercept += 1;

int n_raw_groups_gr_decay_intercept;
int n_cp_groups_gr_decay_intercept;
array[n_levels + 1] int raw_level_pos_gr_decay_intercept;
array[n_levels + 1] int cp_level_pos_gr_decay_intercept;
(n_raw_groups_gr_decay_intercept, raw_level_pos_gr_decay_intercept,
 n_cp_groups_gr_decay_intercept,  cp_level_pos_gr_decay_intercept) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_level_intercept_gr_decay);

array[n_levels] int gr_decay_slope_mode;
for (lv in 1:n_levels) {
  gr_decay_slope_mode[lv] = enable_level_cov_gr_decay[lv] ? enable_level_intercept_gr_decay[lv] : 0;
}
int n_raw_groups_gr_decay_slope;
int n_cp_groups_gr_decay_slope;
array[n_levels + 1] int raw_level_pos_gr_decay_slope;
array[n_levels + 1] int cp_level_pos_gr_decay_slope;
(n_raw_groups_gr_decay_slope, raw_level_pos_gr_decay_slope,
 n_cp_groups_gr_decay_slope,  cp_level_pos_gr_decay_slope) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, gr_decay_slope_mode);
```

- [ ] **Step 2: Wire the include into both models**

In `stan/tumor/sf-ssm-log-space.stan`, add after the `#include "modules/frac/transformed_data.stan"` line (around line 52):
```stan
  #include "modules/gr_decay/transformed_data.stan"
```
In `stan/tumor/sf-ssls-lfo.stan`, add after its `#include "modules/frac/transformed_data.stan"` line (around line 54):
```stan
  #include "modules/gr_decay/transformed_data.stan"
```

- [ ] **Step 3: Syntax check** (same command)

Expected: still FAIL on not-yet-defined parameters (parameters.stan not updated). transformed_data parses.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/gr_decay/transformed_data.stan stan/tumor/sf-ssm-log-space.stan stan/tumor/sf-ssls-lfo.stan
git commit -m "feat(stan/gr_decay): transformed_data flat indices for level hierarchy

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: parameters — level SD, raw/cp intercept, slope params

**Files:**
- Modify: `stan/modules/gr_decay/parameters.stan`

- [ ] **Step 1: Replace the file body** (keep the existing pop intercept + pop-cov coef; ADD frac-style level params)

The existing `gr_decay/parameters.stan` has `gr_decay_log_loc_pop` (array size enable_gr_decay?1:0) and `gr_decay_coef_qr_pop`. KEEP both. Append the level structure, mirroring `frac/parameters.stan` lines 17-32 but gated so they vanish when the master flag is off:

```stan
// ===== UNIFIED LEVEL STRUCTURE (mirrors frac) =====
// All sized by enabled-group counts from transformed_data; collapse to 0 when
// every level mode is NONE (pop-only path) or when enable_gr_decay is off
// (counts are 0 because the flag arrays are passed as all-NONE in that case —
// see r/priors.R / initializer, which set the flag arrays to 0 when off).

// Intercept SD free params — one per RE/RE_CP level.
array[n_re_levels_gr_decay_intercept] real<lower=0> gr_decay_sd_level_intercept_raw;

// Raw NCP and CP intercept draws — split by mode.
vector[n_raw_groups_gr_decay_intercept] gr_decay_raw_level_intercept;
vector[n_cp_groups_gr_decay_intercept]  gr_decay_cp_level_intercept;

// Slope SD hyperparameters — one vector per level.
array[n_levels] vector<lower=0>[n_covar] gr_decay_sd_level_slope;

// Raw NCP and CP slope draws — split by mode.
matrix[n_raw_groups_gr_decay_slope, n_covar] gr_decay_raw_level_slope;
matrix[n_cp_groups_gr_decay_slope,  n_covar] gr_decay_cp_level_slope;
```

**IMPORTANT backward-compat note:** When `enable_gr_decay=0`, the R side MUST pass `enable_level_intercept_gr_decay = rep(0, n_levels)` and `enable_level_cov_gr_decay = rep(0, n_levels)` so all these counts are 0 and the level params vanish — preserving the "off = zero params" invariant. `gr_decay_sd_level_slope` is `array[n_levels]` of `vector[n_covar]`; when `n_covar=0` it is an array of length-0 vectors (harmless). This matches frac (frac always declares it `array[n_levels]`).

- [ ] **Step 2: Syntax check** (same command)

Expected: still FAIL — transformed_parameters.stan and priors.stan reference assembled quantities not yet built. parameters.stan parses.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/gr_decay/parameters.stan
git commit -m "feat(stan/gr_decay): level-indexed intercept + slope parameters

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: transformed_parameters — assemble log_loc_patient + kappa

**Files:**
- Modify: `stan/modules/gr_decay/transformed_parameters.stan`

- [ ] **Step 1: Rewrite to frac-style assembly** (replace the current pop-only body)

Replace `stan/modules/gr_decay/transformed_parameters.stan` entirely with the frac-style assembly. Note the differences from frac: (a) gr_decay's pop intercept is an `array[..]?1:0` so use `gr_decay_log_loc_pop[1]` guarded by enable_gr_decay; (b) final transform is `exp` (kappa>0) not logit-split; (c) when off, set kappa to zeros_vector (the value is never consumed when off — warp sites take the t-branch).

```stan
// gr_decay/transformed_parameters.stan
// Per-patient decay rate kappa_i = exp(log-scale linear predictor):
//   gr_decay_log_loc_pop + pop-cov + Σ_lv level intercepts + Σ_lv level slopes
// Mirrors frac/transformed_parameters.stan. Indexed by forecast-local patient j.
// When enable_gr_decay is off, kappa is zeros (never consumed — warp takes t-branch).

vector[n_forecast_patients] gr_decay_log_loc_patient = zeros_vector(n_forecast_patients);
vector[n_forecast_patients] gr_decay_kappa = zeros_vector(n_forecast_patients);

if (enable_gr_decay) {
  // Population covariate effects
  vector[n_forecast_patients] gr_decay_linpred_pop = enable_pop_cov_gr_decay ?
    (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop) : zeros_vector(n_forecast_patients);

  // ===== SD EXPANSION (frac-style) =====
  array[n_levels] real<lower=0> gr_decay_sd_level_intercept;
  {
    int sd_idx = 0;
    for (lv in 1:n_levels) {
      if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_FE) {
        gr_decay_sd_level_intercept[lv] = gr_decay_fe_sd_level_intercept[lv];
      } else if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE ||
                 enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) {
        sd_idx += 1;
        gr_decay_sd_level_intercept[lv] = gr_decay_sd_level_intercept_raw[sd_idx];
      } else {
        gr_decay_sd_level_intercept[lv] = 0.0;
      }
    }
  }

  // ===== INTERCEPT EFFECTS =====
  vector[n_enabled_groups_gr_decay_intercept] gr_decay_scaled_level_intercept;
  for (lv in 1:n_levels) {
    int mode = enable_level_intercept_gr_decay[lv];
    if (mode == LEVEL_MODE_NONE) continue;
    int e_lo, e_hi;
    (e_lo, e_hi) = get_pos(enabled_level_pos_gr_decay_intercept, lv);
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_gr_decay_intercept[lv];
      int c_hi = cp_level_pos_gr_decay_intercept[lv + 1] - 1;
      gr_decay_scaled_level_intercept[e_lo:e_hi] = gr_decay_cp_level_intercept[c_lo:c_hi];
    } else {
      int r_lo = raw_level_pos_gr_decay_intercept[lv];
      int r_hi = raw_level_pos_gr_decay_intercept[lv + 1] - 1;
      gr_decay_scaled_level_intercept[e_lo:e_hi] =
        gr_decay_sd_level_intercept[lv] * gr_decay_raw_level_intercept[r_lo:r_hi];
    }
  }

  vector[n_forecast_patients] gr_decay_linpred_level_intercepts = zeros_vector(n_forecast_patients);
  for (lv in 1:n_levels) {
    if (enable_level_intercept_gr_decay[lv]) {
      gr_decay_linpred_level_intercepts += gr_decay_scaled_level_intercept[patient_gr_decay_intercept_flat_idx[forecast_patient_idx, lv]];
    }
  }

  // ===== COVARIATE SLOPE EFFECTS =====
  matrix[n_enabled_groups_gr_decay_slope, n_covar] gr_decay_scaled_level_slope;
  if (n_covar > 0 && n_enabled_groups_gr_decay_slope > 0) {
    for (lv in 1:n_levels) {
      if (!enable_level_cov_gr_decay[lv]) continue;
      int mode = enable_level_intercept_gr_decay[lv];
      int e_lo, e_hi;
      (e_lo, e_hi) = get_pos(enabled_level_pos_gr_decay_slope, lv);
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_gr_decay_slope[lv];
        int c_hi = cp_level_pos_gr_decay_slope[lv + 1] - 1;
        gr_decay_scaled_level_slope[e_lo:e_hi, :] = gr_decay_cp_level_slope[c_lo:c_hi, :];
      } else {
        int r_lo = raw_level_pos_gr_decay_slope[lv];
        int r_hi = raw_level_pos_gr_decay_slope[lv + 1] - 1;
        gr_decay_scaled_level_slope[e_lo:e_hi, :] =
          gr_decay_raw_level_slope[r_lo:r_hi, :] .*
          rep_matrix(gr_decay_sd_level_slope[lv]', r_hi - r_lo + 1);
      }
    }
  }

  vector[n_forecast_patients] gr_decay_linpred_level_slopes = zeros_vector(n_forecast_patients);
  if (n_covar > 0) {
    for (lv in 1:n_levels) {
      if (enable_level_cov_gr_decay[lv]) {
        gr_decay_linpred_level_slopes += rows_dot_product(
          Q_covar_design_matrix[forecast_patient_idx, :],
          gr_decay_scaled_level_slope[patient_gr_decay_slope_flat_idx[forecast_patient_idx, lv], :]
        );
      }
    }
  }

  // ===== FINAL LINEAR PREDICTOR =====
  gr_decay_log_loc_patient = rep_vector(gr_decay_log_loc_pop[1], n_forecast_patients)
    + gr_decay_linpred_pop
    + gr_decay_linpred_level_intercepts
    + gr_decay_linpred_level_slopes;
  gr_decay_kappa = exp(gr_decay_log_loc_patient);
}
```

- [ ] **Step 2: Syntax check** (same command)

Expected: still FAIL on priors.stan referencing the new params (until Task 6) — OR may now pass through transformed_parameters cleanly. If it fails ONLY in priors.stan, proceed. If it fails inside transformed_parameters.stan, fix it (likely a typo in a position-array name).

- [ ] **Step 3: Commit**

```bash
git add stan/modules/gr_decay/transformed_parameters.stan
git commit -m "feat(stan/gr_decay): assemble level-indexed log(kappa) + per-patient kappa

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: priors — frac-style unified level loop

**Files:**
- Modify: `stan/modules/gr_decay/priors.stan`

- [ ] **Step 1: Rewrite to frac-style** (the whole block stays inside `if (enable_gr_decay)`)

Replace `stan/modules/gr_decay/priors.stan` with (note: gr_decay does NOT use Student-t, so omit those branches; keep pop intercept + pop cov, then add the level loop):

```stan
// gr_decay/priors.stan — flag-gated priors on log(kappa) with level hierarchy.
// Mirrors frac/priors.stan (minus Student-t). Entire block omitted when off.

if (enable_gr_decay) {
  gr_decay_log_loc_pop[1] ~ normal(gr_decay_log_loc_pop_mean, gr_decay_log_loc_pop_sd);
  if (enable_pop_cov_gr_decay) {
    gr_decay_coef_qr_pop ~ normal(gr_decay_coef_qr_pop_mean, gr_decay_coef_qr_pop_sd);
  }

  int sd_idx = 0;
  for (lv in 1:n_levels) {
    if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE ||
        enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      gr_decay_sd_level_intercept_raw[sd_idx] ~ normal(0, gr_decay_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      gr_decay_sd_level_slope[lv] ~ normal(0, gr_decay_sd_level_slope_sd[lv]);
    }

    // Intercept effects
    if (enable_level_intercept_gr_decay[lv]) {
      if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_gr_decay_intercept[lv];
        int c_hi = cp_level_pos_gr_decay_intercept[lv + 1] - 1;
        gr_decay_cp_level_intercept[c_lo:c_hi] ~ normal(0, gr_decay_sd_level_intercept_raw[sd_idx]);
      } else {
        int r_lo = raw_level_pos_gr_decay_intercept[lv];
        int r_hi = raw_level_pos_gr_decay_intercept[lv + 1] - 1;
        gr_decay_raw_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }

    // Slope effects
    if (enable_level_cov_gr_decay[lv] && n_covar > 0) {
      int mode = enable_level_intercept_gr_decay[lv];
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_gr_decay_slope[lv];
        int c_hi = cp_level_pos_gr_decay_slope[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_covar) {
            gr_decay_cp_level_slope[c_lo:c_hi, k] ~ normal(0, gr_decay_sd_level_slope[lv, k]);
          }
        }
      } else {
        int r_lo = raw_level_pos_gr_decay_slope[lv];
        int r_hi = raw_level_pos_gr_decay_slope[lv + 1] - 1;
        to_vector(gr_decay_raw_level_slope[r_lo:r_hi, :]) ~ std_normal();
      }
    }
  }
}
```

- [ ] **Step 2: Full syntax check — BOTH models**

Run:
```
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: BOTH PASS (exit 0, no output). This is the first point where the full Stan code is consistent.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/gr_decay/priors.stan
git commit -m "feat(stan/gr_decay): frac-style level-hierarchy priors

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 7: R priors — per-level SD hyperparameter defaults

**Files:**
- Modify: `r/priors.R`

- [ ] **Step 1: Add gr_decay level hyperparams to the returned `lst(...)`**

In `r/priors.R`, find the existing gr_decay block in the returned list (currently lines ~374-377: `gr_decay_log_loc_pop_mean`, `_sd`, `gr_decay_coef_qr_pop_mean`, `_sd`). Add immediately after them:

```r
    gr_decay_sd_level_intercept_sd = rep(0.25, n_levels),
    gr_decay_fe_sd_level_intercept = rep(0, n_levels),
    gr_decay_sd_level_slope_sd = list(
      trial = rep(0.05, n_covar),
      patient = rep(0.03, n_covar)
    ),
```

(Mirror frac's values: `0.25` intercept SD prior scale, slope SD scales `0.05`/`0.03`. These are weakly-informative — log-scale SD of 0.25 ⟹ arm/patient κ multipliers mostly within exp(±0.5) ≈ [0.6, 1.65].)

- [ ] **Step 2: Verify the list is reachable** (priors.R uses `tibble::lst` / rlang)

Run: `Rscript -e 'suppressMessages({library(rlang); library(tibble)}); source("r/priors.R"); cat("SOURCED OK\n")'`
Expected: `SOURCED OK` (no parse error). If `get_tumor_priors` requires args to actually call it, just confirming source parses is enough here.

- [ ] **Step 3: Commit**

```bash
git add r/priors.R
git commit -m "feat(r): gr_decay per-level SD hyperparameter defaults

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 8: R initializer — gr_decay raw-level / cp / slope inits

**Files:**
- Modify: `r/sclc/initializers_fixed.R`

- [ ] **Step 1: Compute gr_decay group counts** (near the frac counts, ~line 296-313)

In the `with(stan_data)` block, after the `n_*_groups_frac_*` lines, add:

```r
      n_raw_groups_gr_decay_intercept <- sum(n_groups_per_level[enable_level_intercept_gr_decay %in% c(1L, 2L, 3L)])
      n_cp_groups_gr_decay_intercept  <- sum(n_groups_per_level[enable_level_intercept_gr_decay == 4L])
      n_raw_groups_gr_decay_slope <- sum(n_groups_per_level[(enable_level_intercept_gr_decay %in% c(1L, 2L, 3L)) & (enable_level_cov_gr_decay == 1L)])
      n_cp_groups_gr_decay_slope  <- sum(n_groups_per_level[(enable_level_intercept_gr_decay == 4L) & (enable_level_cov_gr_decay == 1L)])
```

After the `frac_sd_level`/`frac_raw_level`/`frac_cp_level_intercept` lines (~318-328), add:

```r
      gr_decay_sd_level <- rep(0.30, n_levels)
      gr_decay_raw_level <- rnorm(n_raw_groups_gr_decay_intercept, sd = 0.2)
      gr_decay_cp_level_intercept <- rnorm(n_cp_groups_gr_decay_intercept, 0, gr_decay_sd_level[enable_level_intercept_gr_decay == 4L] * 0.2)
```

- [ ] **Step 2: Add the init entries to `tumor_init <- tibble::lst(...)`**

Find the existing gr_decay entries in `tumor_init` (`gr_decay_log_loc_pop`, `gr_decay_coef_qr_pop`, ~lines 339-348). Add after them (gate on enable_gr_decay so off ⟹ zero-length, matching the param sizing):

```r
        gr_decay_sd_level_intercept_raw = if (isTRUE(enable_gr_decay == 1L)) {
          as.array(gr_decay_sd_level[enable_level_intercept_gr_decay %in% c(2L, 4L)])
        } else {
          numeric(0)
        },
        gr_decay_raw_level_intercept = if (isTRUE(enable_gr_decay == 1L)) gr_decay_raw_level else numeric(0),
        gr_decay_cp_level_intercept = if (isTRUE(enable_gr_decay == 1L)) gr_decay_cp_level_intercept else numeric(0),
        gr_decay_sd_level_slope = if (isTRUE(enable_gr_decay == 1L) && n_covar > 0) {
          lapply(seq_len(n_levels), function(lv) rep(if (lv == n_levels) 0.03 else 0.05, n_covar))
        } else {
          lapply(seq_len(n_levels), function(lv) numeric(0))
        },
        gr_decay_raw_level_slope = if (isTRUE(enable_gr_decay == 1L) && n_covar > 0) matrix(0, n_raw_groups_gr_decay_slope, n_covar) else matrix(0, 0, max(n_covar, 0)),
        gr_decay_cp_level_slope = if (isTRUE(enable_gr_decay == 1L) && n_covar > 0) matrix(0, n_cp_groups_gr_decay_slope, n_covar) else matrix(0, 0, max(n_covar, 0)),
```

**Note:** `gr_decay_sd_level_slope` is `array[n_levels] vector[n_covar]` in Stan — when n_covar=0 each element is a length-0 vector, so `lapply(..., numeric(0))` over n_levels is correct (NOT zero-length list). Match frac's `tr_sd_level_slope` shape exactly.

- [ ] **Step 3: Source-parse check**

Run: `Rscript -e 'invisible(parse("r/sclc/initializers_fixed.R")); cat("PARSE OK\n")'`
Expected: `PARSE OK`

- [ ] **Step 4: Commit**

```bash
git add r/sclc/initializers_fixed.R
git commit -m "feat(r): gr_decay level-hierarchy initializer values

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 9: Backward-compat — pass all-NONE flag arrays when off; set publication config

**Files:**
- Modify: `targets/publication_targets.R`

- [ ] **Step 1: Add gr_decay level flags to the publication stan_data config** (near lines 390-396, after `enable_pop_cov_gr_decay`)

The flag arrays must ALWAYS be present in stan_data (the Stan data block declares them unconditionally). Add:

```r
        enable_level_intercept_gr_decay = c(
          trial_arm = level_intercept_mode["re"],
          patient   = level_intercept_mode["re"]
        ),
        enable_level_cov_gr_decay = c(trial_arm = FALSE, patient = FALSE),
```

- [ ] **Step 2: Ensure OFF models still get all-NONE arrays**

Search the codebase for where `enable_gr_decay` is set to 0 in OTHER pipelines (sclc_targets.R, pioneer_targets.R, any test fixtures) and where the tumor stan_data is assembled. The Stan data block now REQUIRES `enable_level_intercept_gr_decay[n_levels]` and `enable_level_cov_gr_decay[n_levels]` to always be present.

Run: `grep -rn "enable_gr_decay" targets/ r/ tests/` and for EVERY stan_data assembly that sets `enable_gr_decay`, confirm the two new flag arrays are also set. Where `enable_gr_decay = 0L`, set both to `rep(0L, n_levels)` (or `c(trial_arm=0L, patient=0L)`) so the off-path has zero gr_decay level params. If a central helper assembles tumor stan_data (e.g. `prepare_tumor_stan_data`), add a default there: `enable_level_intercept_gr_decay = rep(0L, n_levels)`, `enable_level_cov_gr_decay = rep(0L, n_levels)` unless overridden.

**This is the critical backward-compat task** — if any existing OFF model omits these arrays, it won't supply the required data and will error at run. Document what you found and changed.

- [ ] **Step 3: Parse check**

Run: `Rscript -e 'invisible(parse("targets/publication_targets.R")); cat("PARSE OK\n")'`
Expected: `PARSE OK`

- [ ] **Step 4: Commit**

```bash
git add targets/publication_targets.R [any other files touched in Step 2]
git commit -m "feat(targets): gr_decay trial-arm + patient RE hierarchy config

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 10: Test — reduction to pop-only scalar (backward-compat anchor)

**Files:**
- Create: `tests/testthat/stan/test_gr_decay_hierarchy.stan`
- Create: `tests/testthat/test-stan-gr-decay-hierarchy.R`

- [ ] **Step 1: Write the test model** — a minimal standalone Stan model that includes the gr_decay module fragments with a 2-level hierarchy, exposes `gr_decay_kappa`, and is driven by data so the test can set modes.

Create `tests/testthat/stan/test_gr_decay_hierarchy.stan`. It must `#include` the real module fragments (flags, hyperparams, transformed_data, parameters, transformed_parameters) plus the minimal hierarchy scaffolding they depend on (`hierarchy.stanfunctions`, `_hierarchy_data.stan` subset, `_hierarchy_transformed_data.stan`, `Q_covar_design_matrix`, `forecast_patient_idx`, `n_forecast_patients`, `n_forecast_groups_per_level`). Look at how `test_gr_decay_reduction.stan` scaffolds (it inlines a minimal data block) and EXTEND it to include the real `modules/gr_decay/transformed_data.stan` + `transformed_parameters.stan` rather than inlining. If full module include proves too heavy to scaffold, inline the assembly but keep it faithful.

The model exposes `generated quantities { vector[n_forecast_patients] kappa_out = gr_decay_kappa; }`.

- [ ] **Step 2: Write the R test** — TWO configs, assert reduction.

Create `tests/testthat/test-stan-gr-decay-hierarchy.R`:

```r
test_that("gr_decay kappa reduces to pop scalar when all level modes are NONE", {
  # Config A: pop-only (all level modes NONE) — every patient kappa must equal exp(log_loc_pop)
  # Config B: trial-arm + patient RE with zero raw draws — must ALSO equal exp(log_loc_pop)
  #           (RE with raw=0 contributes exactly 0 on the log scale → same kappa).
  # This proves the hierarchy is a strict generalization: zero RE effect ⟺ pop-only.
  # Use a fixed log_loc_pop (e.g. -3.9), n_covar = 0, 2 levels (trial_arm, patient).
  # Compile with cmdstanr include_paths = c("stan", "stan/tumor"), force_recompile,
  # verify symbols, run with fixed_param or optimize so kappa is deterministic.
  # Assert: max(abs(kappa_A - exp(-3.9))) < 1e-9 ; max(abs(kappa_B - kappa_A)) < 1e-9.
})
```

Fill in the actual cmdstanr compile + run. Mirror the structure of `tests/testthat/test-stan-gr-decay-reduction.R` for the compile/run boilerplate (include_paths, force_recompile, the stale-exe guard). Use `expect_lt` with tolerance 1e-9 for the strictly-zero RE reduction (raw=0 exactly ⟹ no CSV round-trip on the difference) — but if driven through CSV, loosen to 1e-6 per the known precision note.

- [ ] **Step 3: Run the test**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-gr-decay-hierarchy.R")'`
Expected: PASS. If the test model won't compile, debug the scaffolding (most likely a missing data field the module needs). Document the stale-exe guard: `rm` the test exe + `force_recompile=TRUE`.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/stan/test_gr_decay_hierarchy.stan tests/testthat/test-stan-gr-decay-hierarchy.R
git commit -m "test(stan): gr_decay hierarchy reduces to pop scalar at zero RE

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 11: End-to-end compile verification (both production models, g++)

**Files:** none (verification only)

- [ ] **Step 1: Compile BOTH production models with cmdstanr** (catches the stale-exe trap that stanc can't)

```r
Rscript -e '
library(cmdstanr); set_cmdstan_path("~/.cmdstan/cmdstan-2.38.0")
for (m in c("stan/tumor/sf-ssm-log-space.stan", "stan/tumor/sf-ssls-lfo.stan")) {
  exe <- sub("\\.stan$", "", m)
  if (file.exists(exe)) file.remove(exe)            # stale-exe guard
  mod <- cmdstan_model(m, include_paths = c("stan", "stan/tumor"), force_recompile = TRUE)
  cat(m, "COMPILED:", mod$exe_file(), "\n")
}'
```
Expected: both print COMPILED with no g++ errors.

- [ ] **Step 2: Verify gr_decay hierarchy symbols are in the binaries**

```bash
strings stan/tumor/sf-ssm-log-space | grep -c "gr_decay_raw_level_intercept\|gr_decay_sd_level_intercept"
strings stan/tumor/sf-ssls-lfo | grep -c "gr_decay_raw_level_intercept\|gr_decay_sd_level_intercept"
```
Expected: both > 0 (symbols present → the new code is actually compiled in, not a stale binary).

- [ ] **Step 3: Run the FULL gr_decay test suite** (regression — the original 4 tests must still pass)

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "gr.decay")'`
Expected: all gr_decay tests PASS (warp, branches, reduction, plateau, hierarchy). The original 4 tests inline their own gr_decay logic so they should be unaffected; confirm.

- [ ] **Step 4: No commit** (verification task). Record results in the execution summary.

---

## Self-Review

**Spec coverage:**
- Trial-arm RE on κ → `enable_level_intercept_gr_decay = c(trial_arm=re, ...)` (Task 9). ✓
- Patient RE on κ → `c(..., patient=re)` (Task 9). ✓
- Covariate path wired, default off → slope params (Task 4), slope assembly (Task 5), slope priors (Task 6), `enable_level_cov_gr_decay = FALSE` (Task 9). ✓
- Backward-compat (off + all-NONE = current scalar) → reduction test (Task 10), all-NONE-when-off wiring (Task 9 Step 2). ✓
- Downstream warp sites unchanged → confirmed κ consumed only via `gr_decay_kappa`; no task touches sf.stanfunctions / state_space / lfo forecast loops. ✓

**Type consistency:** All position-array names (`raw_level_pos_gr_decay_intercept`, `cp_level_pos_gr_decay_slope`, `enabled_level_pos_gr_decay_intercept`) match between transformed_data (Task 3), parameters (Task 4), transformed_parameters (Task 5), priors (Task 6). `gr_decay_sd_level_slope` shape `array[n_levels] vector[n_covar]` consistent across parameters/priors/initializer.

**Placeholder scan:** Task 10's test .stan scaffolding is described, not fully written — flagged as the one judgment-heavy task; the implementer must extend `test_gr_decay_reduction.stan`'s scaffolding faithfully.

**Risk note (identifiability — NOT a code gate):** patient-level RE on κ is weakly identified; the likely empirical failure mode is `gr_decay_sd_level_intercept[patient]` collapsing toward 0 with divergences. This is a MODELING outcome to diagnose post-fit, not a reason to block the implementation. The reduction test only proves the math is a correct generalization; the fit will tell us if the data support patient-level κ variation.
