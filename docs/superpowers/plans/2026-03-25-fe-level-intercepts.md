# Fixed Effects Level Intercepts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add fixed-effects (FE) intercept mode to the multistate level baseline hazard, allowing trial-level hazards to be estimated independently without hierarchical pooling.

**Architecture:** Extend the existing `enable_ms_level_baseline_hazard` flag from 3-state (0=off, 1=RE intercept, 2=RE GP) to 4-state (0=off, 1=FE intercept, 2=RE intercept, 3=RE GP). FE mode reuses the existing NCP raw parameters but scales them by a fixed data hyperparameter instead of an estimated SD.

**Spec:** `docs/superpowers/specs/2026-03-25-fe-level-intercepts-design.md`

**Tech Stack:** Stan, R (targets pipeline), CmdStan 2.38.0

---

### Task 1: Update Flag Encoding in `flags.stan`

**Files:**
- Modify: `stan/modules/multistate/flags.stan:21-22`

- [ ] **Step 1: Update the flag comment and upper bound**

```stan
// 4-state flag: 0 = no level effect, 1 = fixed-effect intercept, 2 = RE intercept, 3 = RE GP
array[n_levels] int<lower=0, upper=3> enable_ms_level_baseline_hazard;
```

- [ ] **Step 2: Commit**

```bash
git add stan/modules/multistate/flags.stan
git commit -m "feat: extend enable_ms_level_baseline_hazard to 4-state (0=off, 1=FE, 2=RE, 3=GP)"
```

---

### Task 2: Add FE Hyperparameters in `hyperparams.stan`

**Files:**
- Modify: `stan/modules/multistate/hyperparams.stan:38-44`

- [ ] **Step 1: Add FE prior SD fields at end of level hyperparameters section**

Insert after line 92 (after the last 3→2 level GP hyperparameter `log_lambda_gp_32_s_level_rho_beta`), before the covariate section (line 94):

```stan
// --- Fixed-Effect Prior SD for Level Intercepts (used when flag == 1) ---
array[n_levels] real<lower=0> fe_log_lambda_gp_01_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_02_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_12_s_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_12_t_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_03_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_32_s_level_intercept_sd;
```

- [ ] **Step 2: Commit**

```bash
git add stan/modules/multistate/hyperparams.stan
git commit -m "feat: add fe_level_intercept_sd data hyperparameters for FE mode"
```

---

### Task 3: Update GP Threshold in `transformed_data.stan`

**Files:**
- Modify: `stan/modules/multistate/transformed_data.stan:58-63`

- [ ] **Step 1: Change GP detection from `== 2` to `== 3`**

Update line 62:
```stan
  ms_level_baseline_is_gp[lv] = (enable_ms_level_baseline_hazard[lv] == 3) ? 1 : 0;
```

Update the comment at line 58-59:
```stan
// For intercept-only modes (FE=1, RE=2), we don't need eta vectors.
// GP-only arrays track which levels use full GP (mode==3).
```

- [ ] **Step 2: Add `any_re_level` flag for conditional parameter sizing**

Insert after line 74 (after the `n_gp_groups_ms_baseline_*` declarations):

```stan
// --- Any-RE Flag (gates conditional sizing of estimated level_intercept_sd) ---
int any_re_level = max(to_array_1d(enable_ms_level_baseline_hazard)) >= 2 ? 1 : 0;
```

This is safe because `enable_ms_level_baseline_hazard` always has length `n_levels >= 1`.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/multistate/transformed_data.stan
git commit -m "feat: update GP threshold == 2 to == 3, add any_re_level flag"
```

---

### Task 4: Conditional Parameter Sizing in `parameters.stan`

**Files:**
- Modify: `stan/modules/multistate/parameters.stan:19,49,78,91,112,129`

- [ ] **Step 1: Gate `level_intercept_sd` on `any_re_level`**

For each transition, change the `level_intercept_sd` sizing. Example for 0→1 (line 19):

Before:
```stan
array[enable_ms_01 ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_intercept_sd;
```

After:
```stan
array[enable_ms_01 && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_intercept_sd;
```

Apply the same pattern to all 6 transitions:
- Line 19: `log_lambda_gp_01_level_intercept_sd` — add `&& any_re_level`
- Line 49: `log_lambda_gp_02_level_intercept_sd` — add `&& any_re_level`
- Line 78: `log_lambda_gp_12_s_level_intercept_sd` — add `&& any_re_level`
- Line 91: `log_lambda_gp_12_t_level_intercept_sd` — add `&& any_re_level`
- Line 112: `log_lambda_gp_03_level_intercept_sd` — add `&& any_re_level`
- Line 129: `log_lambda_gp_32_s_level_intercept_sd` — add `&& any_re_level`

**Do NOT gate `level_alpha` and `level_rho`** — leave their sizing unchanged (`enable_ms_XX ? n_levels : 0`). These GP kernel parameters have unconditional priors in `priors.stan` and are harmless when unused. Gating them would require matching prior gates, adding complexity for no correctness benefit.

- [ ] **Step 2: Commit**

```bash
git add stan/modules/multistate/parameters.stan
git commit -m "feat: gate RE level_intercept_sd and GP kernel params on any_re_level"
```

---

### Task 5: Conditional Scaling in `transformed_parameters.stan`

**Files:**
- Modify: `stan/modules/multistate/transformed_parameters.stan` — 6 transition blocks

- [ ] **Step 1: Update the 0→1 block (lines ~35-71)**

Change the scaling logic. Before (line 35-42):
```stan
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      // Scale intercepts (shared by both intercept-only and GP modes)
      log_lambda_gp_01_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
        log_lambda_gp_01_level_intercept_sd[lv];
```

After:
```stan
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      // Scale intercepts: FE uses fixed data SD, RE uses estimated SD
      if (enable_ms_level_baseline_hazard[lv] == 1) {
        // Fixed effects: no pooling
        log_lambda_gp_01_level_intercept[lv_start:lv_end] =
          raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
          fe_log_lambda_gp_01_level_intercept_sd[lv];
      } else {
        // Random effects (2 or 3): hierarchical pooling
        log_lambda_gp_01_level_intercept[lv_start:lv_end] =
          raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
          log_lambda_gp_01_level_intercept_sd[lv];
      }
```

Then change the GP branch from `== 2` to `== 3` (line ~44):
```stan
      if (enable_ms_level_baseline_hazard[lv] == 3) {
```

- [ ] **Step 2: Repeat for the 0→2 block (lines ~173-205)**

Same pattern: add FE/RE conditional for scaling, change GP check from `== 2` to `== 3`.

- [ ] **Step 3: Repeat for 1→2 sojourn block (lines ~287-330)**

Same pattern. GP check at ~295 changes from `== 2` to `== 3`.

- [ ] **Step 4: Repeat for 1→2 clock-forward block (lines ~376-420)**

Same pattern. GP check at ~384 changes from `== 2` to `== 3`.

- [ ] **Step 5: Repeat for 0→3 block (lines ~472-510)**

Same pattern. GP check at ~480 changes from `== 2` to `== 3`.

- [ ] **Step 6: Repeat for 3→2 block (lines ~538-578)**

Same pattern. GP check at ~546 changes from `== 2` to `== 3`.

- [ ] **Step 7: Commit**

```bash
git add stan/modules/multistate/transformed_parameters.stan
git commit -m "feat: add FE/RE conditional scaling, update GP threshold to == 3"
```

---

### Task 6: Update Priors in `priors.stan`

**Files:**
- Modify: `stan/modules/multistate/priors.stan`

- [ ] **Step 1: Gate RE SD priors on flag >= 2**

For each transition's prior block, add a guard on `level_intercept_sd` priors. The existing code applies `level_intercept_sd` priors **unconditionally** (no flag check). Add a `>= 2` guard so FE mode doesn't try to access the zero-sized parameter.

Example for 0→1 (lines ~29-33). The existing code has NO condition — add one:

Before:
```stan
  if (enable_ms_01) {
    // ...
    log_lambda_gp_01_level_intercept_sd ~ normal(
      0, log_lambda_gp_01_level_intercept_sd_sd
    );
```

After — wrap the `level_intercept_sd` prior in an `any_re_level` guard:
```stan
  if (enable_ms_01) {
    // ...
    if (any_re_level) {
      log_lambda_gp_01_level_intercept_sd ~ normal(
        0, log_lambda_gp_01_level_intercept_sd_sd
      );
    }
```

Apply to all 6 transitions (01, 02, 12_s, 12_t, 03, 32_s). The `raw_log_lambda_gp_XX_level_intercept ~ std_normal()` lines remain unchanged (used by both FE and RE).

**Do NOT change `level_alpha`/`level_rho` priors** — these parameters are not conditionally sized (see Task 4), so their priors remain unconditional.

- [ ] **Step 2: Commit**

```bash
git add stan/modules/multistate/priors.stan
git commit -m "feat: gate RE SD priors on flag >= 2, GP kernel priors on flag == 3"
```

---

### Task 7: Stan Syntax Check

- [ ] **Step 1: Verify all three Stan models compile**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/ms-standalone.stan
```

Expected: all three compile without errors.

- [ ] **Step 2: Fix any compilation errors and re-check**

- [ ] **Step 3: Commit any fixes**

---

### Task 8: R Priors — Add FE Defaults

**Files:**
- Modify: `r/priors.R` — `get_multistate_priors()` (~line 43)

- [ ] **Step 1: Add FE defaults to `get_multistate_priors()`**

After the existing `log_lambda_gp_XX_level_intercept_sd_sd` entries, add:

```r
# Fixed-effect prior SD for level intercepts (used when flag == 1)
fe_log_lambda_gp_01_level_intercept_sd = rep(1.0, n_levels),
fe_log_lambda_gp_02_level_intercept_sd = rep(1.0, n_levels),
fe_log_lambda_gp_12_s_level_intercept_sd = rep(1.0, n_levels),
fe_log_lambda_gp_12_t_level_intercept_sd = rep(1.0, n_levels),
fe_log_lambda_gp_03_level_intercept_sd = rep(1.0, n_levels),
fe_log_lambda_gp_32_s_level_intercept_sd = rep(1.0, n_levels),
```

Note: `get_tumor_priors()` already delegates to `get_multistate_priors()` via
`list_assign(!!!get_multistate_priors(...))` at line 423, so no separate addition needed there.

- [ ] **Step 2: Commit**

```bash
git add r/priors.R
git commit -m "feat: add fe_level_intercept_sd defaults (1.0) to priors"
```

---

### Task 9: R Initializers — Conditional RE/GP Sizing

**Files:**
- Modify: `r/initializers_ms.R:30,55,74,93,112,133,157`
- Modify: `r/initializers_fixed.R:27` (and corresponding lines)

- [ ] **Step 1: Update GP threshold in `initializers_ms.R`**

Line 30 — change `== 2L` to `== 3L`:
```r
n_gp_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard == 3L])
```

- [ ] **Step 2: Conditionally skip `level_intercept_sd` init when all-FE**

Add `any_re_level` computation and gate each `level_intercept_sd` init:

After line 20, add:
```r
any_re_level <- any(enable_ms_level_baseline_hazard >= 2L)
```

Then wrap each `level_intercept_sd` init. Example for 0→1 (line 55):
Before:
```r
log_lambda_gp_01_level_intercept_sd = if (enable_ms_01) abs(rnorm(n_levels, sd = log_lambda_gp_01_level_intercept_sd_sd)) else numeric(0),
```
After:
```r
log_lambda_gp_01_level_intercept_sd = if (enable_ms_01 && any_re_level) abs(rnorm(n_levels, sd = log_lambda_gp_01_level_intercept_sd_sd)) else numeric(0),
```

Apply to all 6 transitions: lines 55, 74, 93, 112, 133, 157.

**Do NOT gate `level_alpha`/`level_rho` init** — leave unchanged (matches Task 4 decision to keep GP kernel parameter sizing unconditional).

- [ ] **Step 3: Apply same changes to `initializers_fixed.R`**

Line 27: change `== 2L` to `== 3L`. Gate `level_intercept_sd` on `any_re_level`.

- [ ] **Step 4: Commit**

```bash
git add r/initializers_ms.R r/initializers_fixed.R
git commit -m "feat: update initializers for FE mode — GP threshold 3L, conditional RE sizing"
```

---

### Task 10: Encoding Migration in Targets Files

**Files:**
- Modify: `targets/sclc_targets.R:~977,~2320`
- Modify: `targets/pioneer_targets.R:~412`

- [ ] **Step 1: Migrate sclc GP config (line ~977)**

Before:
```r
enable_ms_level_baseline_hazard = c(trial = 2L, patient = 0L),
```
After:
```r
enable_ms_level_baseline_hazard = c(trial = 3L, patient = 0L),
```

- [ ] **Step 2: Migrate sclc RE intercept config (line ~2320)**

Before:
```r
enable_ms_level_baseline_hazard = c(trial = 1L, patient = 0L),
```
After:
```r
enable_ms_level_baseline_hazard = c(trial = 2L, patient = 0L),
```

- [ ] **Step 3: Set pioneer to FE mode (line ~412)**

Before:
```r
enable_ms_level_baseline_hazard = as.array(as.integer(c(trial = TRUE, arm = FALSE, patient = FALSE))),
```
After:
```r
enable_ms_level_baseline_hazard = as.array(c(trial = 1L, arm = 0L, patient = 0L)),
```

This sets the trial level to FE intercept (flag=1) — the primary motivation for this feature.

- [ ] **Step 4: Commit**

```bash
git add targets/sclc_targets.R targets/pioneer_targets.R
git commit -m "feat: migrate level baseline hazard encoding (old 1→2, old 2→3), set pioneer to FE"
```

---

### Task 11: Final Syntax Check and Verification

- [ ] **Step 1: Re-run Stan syntax check after all changes**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan && echo "OK: sf-ssm-log-space"

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan && echo "OK: pioneer"

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/ms-standalone.stan && echo "OK: ms-standalone"
```

- [ ] **Step 2: Verify R initializer loads without error**

```bash
Rscript -e "
Sys.setenv(TAR_RUN = 'burden-endpoints')
store <- file.path('/mnt/data/analysis-results', Sys.getenv('DOMINO_STARTING_USERNAME'), 'pioneer/burden-endpoints/_targets')
sd <- targets::tar_read(ms_standalone_stan_data_posterior_combined, store = store)
# Override flag to FE for testing
sd\$enable_ms_level_baseline_hazard <- as.array(c(1L, 0L, 0L))
# Add FE hyperparameters
sd\$fe_log_lambda_gp_01_level_intercept_sd <- rep(1.0, 3)
sd\$fe_log_lambda_gp_02_level_intercept_sd <- rep(1.0, 3)
sd\$fe_log_lambda_gp_12_s_level_intercept_sd <- rep(1.0, 3)
sd\$fe_log_lambda_gp_12_t_level_intercept_sd <- rep(1.0, 3)
sd\$fe_log_lambda_gp_03_level_intercept_sd <- rep(1.0, 3)
sd\$fe_log_lambda_gp_32_s_level_intercept_sd <- rep(1.0, 3)
source('r/initializers_ms.R')
init <- ms_init_values(sd)
cat('Initializer keys:', paste(names(init), collapse = ', '), '\n')
cat('SUCCESS\n')
"
```

- [ ] **Step 3: Fix any issues and commit**
