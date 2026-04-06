# Burden Endpoints — HMC Geometry Ablation Study

**Date**: 2026-04-06  
**Branch**: `karim/burden-endpoints`  
**Question**: Why does combining trial data (FPI-2265-202, 5 arms) with RWD (Flatiron Pluvicto)
cause HMC geometry degradation — tiny stepsizes (0.001–0.003), max treedepth hits — when
the Flatiron-only model runs cleanly (stepsizes 0.015–0.035)?

---

## Background

The baseline `flatiron_markov_level_only` model (Flatiron patients only, Markov 0→1 time
scale, patient-level RE only) sampled with healthy stepsizes ~0.015–0.016.

When trial patients are added (`combined_markov_level_only`), stepsizes collapse to 0.001–0.003.
This ablation series isolates the cause.

---

## Experiments

### Exp 1 — filtered-markov-level-only
**Job**: #489 · Running · 6h 12m  
**Store**: `filtered-markov-level-only`  
**Target**: `pioneer_fit_posterior_filtered_markov_level_only`  
**Purpose**: Establish Flatiron-only healthy baseline with the **filtered** Flatiron dataset
(Pluvicto patients only, stricter eligibility criteria). Flat hierarchy — patient RE only,
no trial or arm intercepts.

**Config**: 1,666 patients (filtered Flatiron) · 2 arms · patient RE · Markov · latent PSA · 16 covariates

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 48% warmup | -9815 | 0.0123 | 12571 | 0 | ✓ healthy |
| 2 | 42% warmup | -9862 | 0.0179 | 12586 | 0 | ✓ healthy |
| 3 | 45% warmup | -9764 | 0.0315 | 12562 | 0 | ✓ healthy |
| 4 | 43% warmup | -9647 | 0.0256 | 12432 | 0 | ✓ healthy |

**Verdict**: ✅ Healthy. Stepsizes 0.012–0.032, zero divergences. Reference baseline.

---

### Exp 2 — unfiltered-markov-level-only
**Job**: #490 · Running · 5h 28m  
**Store**: `unfiltered-markov-level-only`  
**Target**: `pioneer_fit_posterior_unfiltered_markov_level_only`  
**Purpose**: Test whether **Flatiron patient filtering** explains the healthy baseline.
Uses unfiltered Flatiron dataset (+436 patients with broader case mix). Same flat hierarchy.

**Config**: 2,102 patients (unfiltered Flatiron) · 2 arms · patient RE only · Markov · latent PSA

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 17% warmup | -11820 | 0.0255 | 15331 | 0 | ✓ OK |
| 2 | 14% warmup | -11920 | **0.0037** | 15325 | 1 | ⚠️ slow |
| 3 | 16% warmup | -11739 | 0.0238 | 15202 | 0 | ✓ OK |
| 4 | 8% warmup | -12203 | **0.0016** | 15780 | 0 | ⚠️ struggling |

**Verdict**: ⚠️ Mixed. Two chains healthy, two struggling. Filtering helps but isn't the whole story.

---

### Exp 3 — combined-markov-level-observed
**Job**: #492 · Running · 4h 20m  
**Store**: `combined-markov-level-observed`  
**Target**: `pioneer_fit_posterior_combined_markov_level_observed`  
**Purpose**: Test whether **latent vs observed PSA** causes the geometry problem in the
combined model. Uses **observed PSA values** directly as the 0→1 transition covariate
(no latent state space for PSA), with trial + Flatiron combined data.

**Config**: ~1,773 patients (combined, filtered Flatiron) · 7 arms · patient RE only · Markov · **observed PSA**

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 14% warmup | -10228 | 0.0202 | 13154 | 0 | ✓ OK |
| 2 | 13% warmup | -10171 | 0.0111 | 13122 | 0 | ✓ OK |
| 3 | 16% warmup | -10205 | 0.0228 | 13107 | 0 | ✓ OK |
| 4 | 13% warmup | -10159 | **0.0026** | 13053 | 1 | ⚠️ slow |

**Verdict**: ⚠️ Mostly OK. Three chains healthy (0.011–0.023), one struggling. Observed PSA
helps but chain 4 still slow — combined data geometry remains partially problematic.

---

### Exp 4 — combined-markov-level-only (full joint model)
**Job**: #492 (shared store) · Running · 4h 17m  
**Store**: `combined-markov-level-only`  
**Target**: `pioneer_fit_posterior_combined_markov_level_only`  
**Purpose**: **Primary problem case.** Full joint model (PSA + multistate) on combined
trial + Flatiron data, flat hierarchy (patient RE only), latent PSA. Reproduces the
geometry degradation seen in production `combined`.

**Config**: ~2,102 patients (trial + unfiltered Flatiron) · 7 arms · patient RE only · Markov · latent PSA

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 9% warmup | -10374 | **0.0022** | 13305 | 0 | ⚠️ stuck |
| 2 | 7% warmup | -10135 | **0.0012** | 13060 | 1 | ⚠️ stuck |
| 3 | 9% warmup | -10460 | **0.0008** | 13474 | 2 | ⚠️ stuck |
| 4 | 9% warmup | -10127 | **0.0017** | 13085 | 1 | ⚠️ stuck |

**Verdict**: ❌ Stuck. All chains at ~8% warmup after 4h 17m. Stepsizes 0.0008–0.0022.
This is the problem we are trying to fix.

---

### Exp 5 — PSA standalone (combined, flat hierarchy)
**Job**: #493 · Running · 4h 17m  
**Store**: `combined-markov-level-only`  
**Target**: `pioneer_psa_standalone_combined_markov_level_only`  
**Purpose**: **Isolate the PSA module.** Run PSA dynamics alone (no multistate hazard)
on combined data with flat hierarchy. Tests whether the geometry problem is in the PSA
state-space module or the joint model.

**Config**: Same data as Exp 4 · patient RE only · PSA dynamics only (no multistate fit)

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 15% warmup | -5668 | 0.0178 | 8642 | 0 | ✓ OK |
| 2 | 14% warmup | -5618 | 0.0221 | 8448 | 0 | ✓ OK |
| 3 | 14% warmup | -5663 | 0.0314 | 8605 | 0 | ✓ OK |
| 4 | 14% warmup | -5667 | **0.0064** | 8585 | 1 | ⚠️ slow |

**Verdict**: ⚠️ Mostly OK. Three chains healthy (0.018–0.031), chain 4 slow. PSA-only
combined runs better than the full joint model — the multistate module interaction with
combined data is partly responsible.

---

### Exp 6 — PSA standalone (arm-level RE)
**Job**: #488 · Running · ~6h+  
**Store**: `student-t-standalone`  
**Target**: `pioneer_psa_standalone_arm_markov_level_only`  
**Purpose**: Test whether adding **arm-level RE** to the PSA standalone model (combined
data) helps or hurts. Arm intercepts allow each of the 7 treatment arms to have its own
PSA baseline.

**Config**: ~1,773 patients · 7 arms · **arm RE + patient RE** · Markov · latent PSA

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 29% warmup | -8230 | **0.0016** | 11178 | 0 | ⚠️ slow |
| 2 | 30% warmup | -9060 | **0.0077** | 12040 | 0 | ⚠️ slow |
| 3 | 31% warmup | -8386 | **0.0090** | 11331 | 0 | ⚠️ slow |
| 4 | 29% warmup | -8855 | **0.0075** | 11751 | 0 | ⚠️ slow |

**Verdict**: ⚠️ Slow. RE arm-level does NOT improve geometry — stepsizes 0.002–0.009.
Adding a random-effect SD parameter for arm creates funnel geometry.

---

### Exp 7 — Full joint model (arm-level RE, observed PSA)
**Job**: #494 · Running · 4h 14m  
**Store**: `student-t`  
**Target**: `pioneer_fit_posterior_arm_markov_level_only_obs_psa`  
**Purpose**: Full joint PSA + multistate model, combined data, **arm-level RE** hierarchy,
observed PSA. Tests whether arm RE hierarchy + observed PSA together help.

**Config**: ~1,773 patients · **arm RE + patient RE** · Markov · **observed PSA** · 16 covariates

**Diagnosis** _(last updated: 2026-04-06 ~21:00 UTC)_

| Chain | Progress | lp | stepsize | energy | divs | Status |
|-------|----------|----|----------|--------|------|--------|
| 1 | 8% warmup | -14310 | **0.0026** | 17560 | 0 | ⚠️ stuck |
| 2 | 7% warmup | -15297 | **0.0023** | 18430 | 1 | ⚠️ stuck |
| 3 | 7% warmup | -13531 | **0.0018** | 16570 | 0 | ⚠️ stuck |
| 4 | 7% warmup | -14347 | **0.0023** | 17303 | 0 | ⚠️ stuck |

**Verdict**: ❌ Stuck. Arm RE + observed PSA makes things worse. The funnel from the arm-level
SD parameter dominates. Stepsizes ~0.002 — comparable to the worst combined models.

---

### Exp 8 — combined-markov-fe-arm-level-only (KEY EXPERIMENT)
**Job**: #504 · Running · <1h  
**Store**: `combined-markov-fe-arm-level-only`  
**Target**: `pioneer_fit_posterior_combined_markov_fe_arm_level_only`  
**Purpose**: **Primary fix candidate.** Full joint model (PSA + multistate), combined
trial + Flatiron data, **arm-level FIXED EFFECTS** (SD is a fixed hyperparameter, not
estimated). FE avoids the SD funnel geometry while still allowing each arm its own
PSA intercept.

Theory: The population PSA mean cannot simultaneously describe RWD Pluvicto (slow responders)
and 5 trial dose arms (dose-dependent fast responses). Arm-level intercepts let each arm
offset from the population mean. FE avoids the RE SD funnel that degraded geometry in Exp 6/7.

**Config**: ~2,102 patients (combined) · **arm FE + patient RE** · Markov · latent PSA  
**FE SD hyperpriors**: tr=0.50, frac=0.35, init=0.50

**Diagnosis**: No CSVs yet — just launched. Check back in ~30–60 min.

---

## Failed Launch Attempts (infrastructure only)

Jobs #495–503 were failed attempts to launch Exp 8 due to:
- **#495**: missing `-s` flag
- **#496**: `NA` in `tr_raw_level_intercept` — `tr_sd_level` had 2 elements for 3-level hierarchy (fixed: `855cee2`)
- **#497**: Same NA error — fix was applied to wrong file (`initializers_fixed.R` instead of `r/pioneer/initializers.R`)
- **#498–499**: Wrong API field: used `commitId` (silently ignored) instead of `mainRepoGitRef`
- **#500–501**: Job ran against `main` branch (mirror lag / no SHA pinning)
- **#502**: `tr_raw_level_intercept` NA still present — three bugs in `r/pioneer/initializers.R`:
  mode-as-multiplier in `cumsum()`, wrong parameter name (`tr_sd_level_intercept` → `tr_sd_level_intercept_raw`), `==1` instead of `!=0` (fixed: `162aa30`)
- **#503**: Wrong branch (Domino used `main` mirror, not `karim/burden-endpoints`)
- **#504 attempt #1 / job #502**: Ragged list error — `tr_fe_sd_level_intercept` missing from `get_pioneer_priors()`, so `[[2]] <- 0.5` on NULL created `list(NULL, 0.5)` (fixed: `f891c197`)

---

## Summary Table

| Exp | Job | Store | Patients | arm hier | PSA | Status | Stepsize range |
|-----|-----|-------|----------|----------|-----|--------|---------------|
| 1 | #489 | filtered-markov-level-only | 1,666 (filtered RWD) | none | latent | ✅ Running | 0.012–0.032 |
| 2 | #490 | unfiltered-markov-level-only | 2,102 (RWD only) | none | latent | ⚠️ Running | 0.002–0.026 |
| 3 | #492 | combined-markov-level-observed | ~1,773 (combined) | none | observed | ⚠️ Running | 0.003–0.023 |
| 4 | #492 | combined-markov-level-only | ~2,102 (combined) | none | latent | ❌ Running | 0.0008–0.002 |
| 5 | #493 | combined-markov-level-only | ~2,102 (combined) | none | latent (PSA standalone) | ⚠️ Running | 0.006–0.031 |
| 6 | #488 | student-t-standalone | ~1,773 (combined) | RE | latent (PSA standalone) | ⚠️ Running | 0.002–0.009 |
| 7 | #494 | student-t | ~1,773 (combined) | RE | observed | ❌ Running | 0.002–0.003 |
| 8 | #504 | combined-markov-fe-arm-level-only | ~2,102 (combined) | **FE** | latent | 🔄 Starting | TBD |

**Key finding so far**: The geometry problem is caused by the combination of (a) trial patients
with dose-dependent PSA kinetics and (b) a single population mean forced to span all 7 treatment
arms. Arm-level RE makes it worse (SD funnel). Arm-level FE (Exp 8) is the key test.
