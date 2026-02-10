# Inventory & Refactor Spec: Total Rate (tr), Fraction Mix (frac), Initial Proportions (init)

Scope: Parameters and related hyperparameters/priors involved in (a) total tumor rate magnitude, (b) decrease vs growth rate fraction (mix), and (c) initial state proportions. Excludes GP/process noise, measurement noise, growth lag / transition, and other endpoints. Based on files:

- `stan/tumor/sf-ssls-parameters.stan`
- `stan/tumor/sf-ssls-transformed_parameters.stan`
- `stan/tumor/sf-ssls-priors.stan`
- `stan/tumor/sf-ssls-hyperparam.stan`
- `stan/base_data.stan` (flags)
- Targets pipeline flag assignments in `targets/sclc_targets.R`

Legend columns:
- role: total_rate | frac | init_prop
- level: pop | trial | patient
- type: loc (intercept), coef_qr (QR-space coefficient), coef (original-scale coefficient), sd (hierarchical SD parameter), raw (standard normal latent), effect (realized hierarchical effect), derived (deterministic algebraic including all linpred objects)
- gating flags: data flags required for parameter/effect to influence model. Some objects always allocated but may be zeroed if flags disable flows.
- source: where *declared* (P=parameters, H=hyperparam/data, T=transformed parameters, G=generated quantities, R=prior statement location)

Linpred classification note: Previously listed as its own type; now considered a derived deterministic object (matrix * coefficient + optional hierarchical slope effects) and marked type = derived with a parenthetical (linpred) in mapping tables for clarity.

## 1. Current Parameter Inventory (Unmodified Model)

| # | Name | Role | Level | Type | Decl Src | Prior Src | Hyperparams | Gating Flags (must be TRUE & not negated) | Usage / Notes |
|---|------|------|-------|------|----------|-----------|-------------|-------------------------------------------|---------------|
| 1 | pop_log_total_rate | total_rate | pop | loc | P | R (`sf-ssls-priors`) | pop_log_total_rate_mean, pop_log_total_rate_sd | (always) | Base log scale for total rate |
| 2 | trial_log_total_rate_sd | total_rate | trial | sd | P | R | trial_log_total_rate_sd_sd | add_trial_level_total_rate && !pop_rates_param_only | SD for trial RE |
| 3 | raw_trial_log_total_rate | total_rate | trial | raw | P | R (std_normal) | (scales via #2) | add_trial_level_total_rate && !pop_rates_param_only | Latent trial RE draws |
| 4 | patient_log_total_rate_sd | total_rate | patient | sd | P | R | patient_log_total_rate_sd_sd | !pop_rates_param_only | SD for patient RE |
| 5 | raw_patient_log_total_rate | total_rate | patient | raw | P | R (std_normal) | (scales via #4) | !pop_rates_param_only | Latent patient RE draws |
| 6 | trial_log_total_rate_effect | total_rate | trial | effect | T | (implicit) | (#2,#3) | add_trial_level_total_rate && !pop_rates_param_only | Added to patient total rate |
| 7 | patient_log_total_rate_effect | total_rate | patient | effect | T | (implicit) | (#4,#5) | !pop_rates_param_only | Patient total rate deviation |
| 8 | patient_log_total_rate | total_rate | patient | derived loc | T | — | uses 1,6,7 | (see above) | Input to per-process rates |
| 9 | pop_decrease_frac_logit | frac | pop | loc (logit) | P | R | pop_decrease_frac_logit_mean, pop_decrease_frac_logit_sd | (always) | Governs allocation decrease vs growth |
| 10 | pop_decrease_frac_logit_coef_qr | frac | pop | coef_qr | P | (prior on #11) | pop_decrease_frac_logit_coef_mean, pop_decrease_frac_logit_coef_sd | (always) | QR-space coeff vector |
| 11 | pop_decrease_frac_logit_coef | frac | pop | coef | T | R | (derived via backsolve) | (always) | Original-scale betas (prior applied here) |
| 12 | patient_decrease_frac_logit_linpred | frac | patient | linpred | T | — | (#10) | (always) | Q * coef_qr; added to patient logit |
| 13 | raw_trial_decrease_frac_logit_coef | frac | trial | raw | P | R (std_normal) | trial_decrease_frac_logit_coef_sd (scales via #14) | (intended: add_trial_level_total_rate analog? none defined) | Currently NOT applied (transformation commented) -> orphan risk |
| 14 | trial_decrease_frac_logit_coef_sd | frac | trial | sd (per-cov) | P | R | trial_decrease_frac_logit_coef_sd_sd | (same as #13) | Prior sampled though effect unused |
| 15 | raw_patient_decrease_frac_logit_coef | frac | patient | raw (matrix) | P | R (std_normal) | patient_decrease_frac_logit_coef_sd (scales via #16) | add_patient_level_frac | Slopes currently not added to linpred (needs confirm) |
| 16 | patient_decrease_frac_logit_coef_sd | frac | patient | sd (per-cov) | P | R | patient_decrease_frac_logit_coef_sd_sd | add_patient_level_frac | Prior sampled; effect may be unused |
| 17 | patient_decrease_frac_logit_sd | frac | patient | sd (intercept RE) | P | R | patient_decrease_frac_logit_sd_sd | add_patient_level_frac | Active intercept variability |
| 18 | raw_patient_decrease_frac_logit | frac | patient | raw (intercept) | P | R (std_normal) | (scaled by #17) | add_patient_level_frac | Added to patient logit |
| 19 | patient_decrease_frac_logit | frac | patient | loc (post-RE) | T | — | uses 9,12,(17,18) | add_patient_level_frac (for RE part) | Core patient fraction logit |
| 20 | patient_log_decrease_frac | frac | patient | derived | T | — | from #19 | (always) | log(f_decrease) |
| 21 | patient_log_growth_frac | frac | patient | derived | T | — | from #19 | (always) | log(1 - f_decrease) |
| 22 | pop_log_decrease_frac | frac | pop | derived | G | — | from #9 | (always) | Generated qty convenience |
| 23 | pop_log_growth_frac | frac | pop | derived | G | — | from #9 | (always) | Generated qty convenience |
| 24 | pop_decrease_prop_logis | init_prop | pop | loc (logit) | P | R | pop_decrease_prop_logis_mean, pop_decrease_prop_logis_sd | (always) | Initial decrease-state proportion |
| 25 | pop_decrease_prop_logis_coef_qr | init_prop | pop | coef_qr | P | (prior on #26) | pop_decrease_prop_logis_coef_mean, pop_decrease_prop_logis_coef_sd | (always) | QR coefficients |
| 26 | pop_decrease_prop_logis_coef | init_prop | pop | coef | T | R | derived via backsolve | (always) | Original-scale betas |
| 27 | patient_decrease_prop_logis_linpred | init_prop | patient | linpred | T | — | (#25) | (always) | Q * coef_qr contribution |
| 28 | trial_decrease_prop_logis_sd | init_prop | trial | sd | P | R | trial_decrease_prop_logis_sd_sd | add_trial_level_prop && !pop_initial_states_param_only | SD for trial proportion loc |
| 29 | raw_trial_decrease_prop_logis | init_prop | trial | raw | P | R (std_normal) | (scaled by #28) | add_trial_level_prop && !pop_initial_states_param_only | Trial RE draws |
| 30 | patient_decrease_prop_logis_sd | init_prop | patient | sd | P | R | patient_decrease_prop_logis_sd_sd | !pop_initial_states_param_only | SD for patient proportion loc |
| 31 | raw_patient_decrease_prop_logis | init_prop | patient | raw | P | R (std_normal) | (scaled by #30) | !pop_initial_states_param_only | Patient RE draws |
| 32 | trial_decrease_prop_logis_effect | init_prop | trial | effect | T | — | (#28,#29) | add_trial_level_prop && !pop_initial_states_param_only | Added to patient values |
| 33 | patient_decrease_prop_logis_effect | init_prop | patient | effect | T | — | (#30,#31) | !pop_initial_states_param_only | Added to patient values |
| 34 | patient_decrease_prop_logis | init_prop | patient | loc (post-RE) | T | — | 24 + linpred + effects | see gating above | Core patient initial proportion logit |
| 35 | patient_log_decrease_prop | init_prop | patient | derived | T | — | from #34 | (always) | log(p_decrease_init) |
| 36 | patient_log_growth_prop | init_prop | patient | derived | T | — | from #34 | (always) | log(1 - p_decrease_init) |
| 37 | patient_log_decrease_rate | total_rate+frac | patient | derived | T | — | patient_log_total_rate + patient_log_decrease_frac | (always) | Decrease process log-rate |
| 38 | patient_log_growth_rate | total_rate+frac | patient | derived | T | — | patient_log_total_rate + patient_log_growth_frac | (always) | Growth process log-rate |
| 39 | pop_log_decrease_rate | total_rate+frac | pop | derived | G | — | #1 + #22 | (always) | Convenience output |
| 40 | pop_log_growth_rate | total_rate+frac | pop | derived | G | — | #1 + #23 | (always) | Convenience output |
| 41 | trial_log_decrease_rate | total_rate+frac | trial | derived | G | — | trial_log_total_rate + pop_log_decrease_frac | add_trial_level_total_rate | Uses implicit trial_total_rate_effect |
| 42 | trial_log_growth_rate | total_rate+frac | trial | derived | G | — | trial_log_total_rate + pop_log_growth_frac | add_trial_level_total_rate |  |
| 43 | trial_log_decrease_rate_residual | total_rate+frac | trial | derived | G | — | trial_log_decrease_rate - pop_log_decrease_rate | add_trial_level_total_rate | Centered residual |
| 44 | trial_log_growth_rate_residual | total_rate+frac | trial | derived | G | — | trial_log_growth_rate - pop_log_growth_rate | add_trial_level_total_rate | Centered residual |
| 45 | patient_log_decrease_rate_residual | total_rate+frac | patient | derived | G | — | patient_log_decrease_rate - pop_log_decrease_rate | !pop_rates_param_only | Residual wrt population |
| 46 | patient_log_growth_rate_residual | total_rate+frac | patient | derived | G | — | patient_log_growth_rate - pop_log_growth_rate | !pop_rates_param_only | Residual wrt population |

## Flags (Relevant Subset Definitions)
- `add_trial_level_total_rate`: Enable trial hierarchy for total rate (trial sd + raw + effect).
- `pop_rates_param_only`: When TRUE, suppresses trial & patient total rate random effects (but population param remains).
- `add_patient_level_frac`: Enables patient-level fraction intercept (and infrastructure for slopes, though slopes currently not applied).
- `pop_covar_coef_only`: In transformation logic, if TRUE, skips hierarchical covariate effects for fraction (trial & patient slopes). Name suggests inverse semantics; contributes to confusion.
- `add_trial_level_prop`: Enable trial hierarchy for initial proportions.
- `pop_initial_states_param_only`: Disables trial & patient initial proportion random effects & related draws when TRUE.

## Orphan / Possibly Unused Parameters
- Fraction trial-level random slopes (#13, #14) sampled but not applied (transformed code lines handling them are commented out). Increases dimension without affecting likelihood.
- Fraction patient-level random slopes (#15, #16) appear similarly unused unless omitted lines (not shown in context) add them; needs confirmation by inspecting omitted lines 28–30 in `sf-ssls-transformed_parameters.stan`.

## Asymmetries & Inconsistencies
- Total rate has both trial and patient REs with clear gating; fraction has only patient intercept RE (slopes infrastructure dormant; no trial intercept RE) producing conceptual imbalance.
- Different gating naming styles (`pop_*_param_only` vs `add_patient_level_frac`) increase cognitive load.
- Mixed naming for logit vs log-derived transformations (e.g., `patient_log_decrease_frac` is log(f) while original variable is on logit scale) can mislead.

## Candidate Cleanup Themes (For Later Steps)
1. Remove or activate dormant fraction hierarchical slopes to eliminate orphan parameters.
2. Harmonize gating: a uniform schema like `enable_trial_fraction`, `enable_patient_fraction_slopes`, etc.
3. Explicitly separate fraction intercept vs slopes in naming.
4. Provide consistent naming for log(fraction) vs logit(fraction).
5. Co-locate each module's (total / fraction / init_prop) hyperparams, priors, parameters, transformed components, and derived outputs.

## Open Questions (Require Confirmation Before Refactor)
- Should fraction module support trial-level variation (intercept and/or slopes)? If yes, revive code; if no, remove parameters (#13,#14).
- Are patient-level fraction covariate slopes intended? If not, remove (#15,#16) and simplify initializer/prior generation.
- Is `pop_covar_coef_only` meant to mean "only population (no hierarchical)"? If so, rename to `disable_hier_fraction` or similar.

## Next Step (Not Executed Yet)
Proceed to classification & mapping (Step 2) once orphan usage is confirmed.

---
Generated: Inventory Step (Total Rate, Fraction, Initial Proportions). No code behavior altered.

---

## 2. Conceptual Framework (Refactor Target)

This section proposes a unified conceptual decomposition and a future naming / structural scheme for the three focused modules.

### 1. Conceptual Modules

We explicitly separate three orthogonal layers that combine to produce process rates:

1. TotalRate: Magnitude baseline (log total rate) with hierarchical variation.
2. FractionMix: Allocation of total rate between decrease and growth processes (on logit scale -> probabilities -> log fractions).
3. InitialStateProp: Allocation of initial latent state proportions (distinct from ongoing process fractions).

Each module optionally supports hierarchical random effects at trial and patient levels and covariate linear predictors (population + optional hierarchical slopes).

### 2. Design Principles

* Uniform naming pattern (see below) reduces cognitive switching between modules.
* Zero-length parameters only when a hierarchy is disabled; no sampling of unused parameters (eliminate orphan priors).
* QR decomposition for covariate coefficients always reflected with `_coef_qr` (QR-space) and `_coef` (original space) naming (replaces earlier placeholder `theta`).
* Clear semantic distinction between logit-scale parameters (`*_logit_loc_*`) and log of a probability (`*_log_*`).
* Hierarchical intercept RE vs hierarchical slope RE separated explicitly.

### 2.1 Naming Pattern (Final)

General token order:

`<module>_<scale?>_<component>_<level>` with optional sub-component tokens.

Where:
* module: `tr` (total rate), `frac` (fraction mix), `init` (initial state proportion)
* scale?: optional `logit`, `log`, `log_decrease`, `log_growth`, etc., to make transformation explicit
* component: one of:
	* `loc` (population or realized intercept / baseline on the given scale)
	* `coef_qr` (QR-space coefficients). Replaces previous term `theta` for clarity.
	* `coef` (original-scale coefficients after backsolve)
	* `linpred` (pure fixed-effects linear predictor = design matrix × QR-space coefficients + hierarchical slope contributions ONLY; never includes intercept or intercept random effects)
	* `sd` (hierarchical SD hyperparameter)
	* `raw` (standard normal latent draws)
	* `effect` (realized random effect = sd * raw)
	* `prob` (on probability scale when derived from logit)
	* residual tokens: `rate_resid`, etc.
* level: `pop`, `trial`, `patient` (omitted when not level-specific, e.g., derived mixture pieces at population level)

Definition of `linpred`:
`<module>_linpred_<level>` represents ONLY the contribution from covariates (fixed + optional hierarchical slopes) before adding any intercept (population) or intercept random effects. We retain log or logit space throughout; conversion to probability or non-log scale happens only when required for downstream quantities or reporting, preserving numerical stability and clarity.

Examples:
* `tr_loc_pop` (current `pop_log_total_rate`)
* `tr_sd_trial`, `tr_raw_trial`, `tr_effect_trial`
* `frac_logit_loc_pop` (current `pop_decrease_frac_logit`)
* `frac_coef_qr_pop` / `frac_coef_pop` (current `pop_decrease_frac_logit_coef_qr` / recovered original scale)
* `frac_linpred_patient` (current `patient_decrease_frac_logit_linpred` minus intercept portions)
* `frac_logit_loc_patient` (current `patient_decrease_frac_logit`)
* `frac_log_decrease_patient` (current `patient_log_decrease_frac`)
* `init_logit_loc_pop` (current `pop_decrease_prop_logis`)
* `init_coef_qr_pop`, `init_coef_pop`
* `init_linpred_patient`

### 2.2 Module & Level Mapping (Existing + New)

Tables below enumerate (a) existing parameters slated for renaming and (b) newly proposed parameters that do NOT yet exist (status = new). Blank Current indicates a new scaffolded element.

#### Legend
Type: loc | coef_qr | coef | linpred | sd | raw | effect | derived
Level: pop | trial | patient
Status: existing (present today) | new (to be added) | remove (to be deleted) | activate (present but unused; will be wired) | prune (will be dropped to reduce dimension)

#### 2.2.1 Total Rate (tr)

Shape notation: scalar, vector[n_patient], vector[n_trial], vector[K_tr], matrix[n_trial,K_tr], matrix[n_patient,K_tr]. (K_tr = number of total rate covariates.)

| Module | Level | Current | Proposed | Type | Shape | Status | Notes |
|--------|-------|---------|----------|------|-------|--------|-------|
| tr | pop | pop_log_total_rate | tr_loc_pop | loc | scalar | existing | Core log intercept |
| tr | pop | — | tr_coef_qr_pop | coef_qr | vector[K_tr] | new | Population covariate QR coef (currently absent) |
| tr | pop | — | tr_coef_pop | coef | vector[K_tr] | new | Original-scale betas |
| tr | pop | — | tr_linpred_pop | derived (linpred) | vector[n_patient] | new | Population covariate contribution broadcast to patients |
| tr | patient | patient_log_total_rate_sd | tr_sd_patient_intercept | sd | scalar | existing | Rename clarifies intercept |
| tr | patient | raw_patient_log_total_rate | tr_raw_patient_intercept | raw | vector[n_patient] | existing | Standard normal draws |
| tr | patient | patient_log_total_rate_effect | tr_effect_patient_intercept | effect | vector[n_patient] | existing | sd * raw |
| tr | patient | patient_log_total_rate | tr_loc_patient | loc | vector[n_patient] | existing | Realized patient intercept |
| tr | trial | trial_log_total_rate_sd | tr_sd_trial_intercept | sd | scalar | existing | Trial intercept SD |
| tr | trial | raw_trial_log_total_rate | tr_raw_trial_intercept | raw | vector[n_trial] | existing | Trial intercept draws |
| tr | trial | trial_log_total_rate_effect | tr_effect_trial_intercept | effect | vector[n_trial] | existing | Trial effect |
| tr | trial | — | tr_loc_trial | loc | vector[n_trial] | new | Realized trial intercept (pop + effect) |
| tr | patient | — | tr_linpred_patient | derived (linpred) | vector[n_patient] | new | Covariate contribution only |
| tr | trial | — | tr_linpred_trial | derived (linpred) | vector[n_trial] | new | For trial-scope summaries if needed |
| tr | trial | — | tr_sd_trial_slope | sd | vector[K_tr] | new | Trial slope SD |
| tr | trial | — | tr_raw_trial_slope | raw | matrix[n_trial,K_tr] | new | Trial slope draws |
| tr | trial | — | tr_effect_trial_slope | effect | matrix[n_trial,K_tr] | new | Applied to linpred (if enabled) |
| tr | patient | — | tr_sd_patient_slope | sd | vector[K_tr] | new | Patient slope SD |
| tr | patient | — | tr_raw_patient_slope | raw | matrix[n_patient,K_tr] | new | Patient slope draws |
| tr | patient | — | tr_effect_patient_slope | effect | matrix[n_patient,K_tr] | new | Patient slope deviation |
| tr | patient | patient_log_decrease_rate | tr_log_decrease_rate_patient | derived | vector[n_patient] | existing | Derived with fraction |
| tr | patient | patient_log_growth_rate | tr_log_growth_rate_patient | derived | vector[n_patient] | existing | Derived with fraction |
| tr | pop | pop_log_decrease_rate | tr_log_decrease_rate_pop | derived | scalar | existing | |
| tr | pop | pop_log_growth_rate | tr_log_growth_rate_pop | derived | scalar | existing | |
| tr | trial | trial_log_decrease_rate | tr_log_decrease_rate_trial | derived | vector[n_trial] | existing | |
| tr | trial | trial_log_growth_rate | tr_log_growth_rate_trial | derived | vector[n_trial] | existing | |
| tr | trial | trial_log_decrease_rate_residual | tr_log_decrease_rate_resid_trial | derived | vector[n_trial] | existing | Centered residual |
| tr | trial | trial_log_growth_rate_residual | tr_log_growth_rate_resid_trial | derived | vector[n_trial] | existing | Centered residual |
| tr | patient | patient_log_decrease_rate_residual | tr_log_decrease_rate_resid_patient | derived | vector[n_patient] | existing | Centered residual |
| tr | patient | patient_log_growth_rate_residual | tr_log_growth_rate_resid_patient | derived | vector[n_patient] | existing | Centered residual |

#### 2.2.2 Fraction Mix (frac)

Shape notation: scalar, vector[n_patient], vector[n_trial], vector[K_frac], matrix[n_trial,K_frac], matrix[n_patient,K_frac]. (K_frac = number of fraction covariates.)

| Module | Level | Current | Proposed | Type | Shape | Status | Notes |
|--------|-------|---------|----------|------|-------|--------|-------|
| frac | pop | pop_decrease_frac_logit | frac_logit_loc_pop | loc | scalar | existing | Population logit intercept |
| frac | pop | pop_decrease_frac_logit_coef_qr | frac_coef_qr_pop | coef_qr | vector[K_frac] | existing | QR coefficients |
| frac | pop | pop_decrease_frac_logit_coef | frac_coef_pop | coef | vector[K_frac] | existing | Original-scale betas (prior target) |
| frac | patient | patient_decrease_frac_logit_linpred | frac_linpred_patient | derived (linpred) | vector[n_patient] | existing | Will exclude intercept once refactored |
| frac | patient | patient_decrease_frac_logit_sd | frac_sd_patient_intercept | sd | scalar | existing | Intercept SD |
| frac | patient | raw_patient_decrease_frac_logit | frac_raw_patient_intercept | raw | vector[n_patient] | existing | Intercept raw |
| frac | patient | patient_decrease_frac_logit | frac_logit_loc_patient | loc | vector[n_patient] | existing | Realized patient logit |
| frac | pop | pop_log_decrease_frac | frac_log_decrease_pop | derived | scalar | existing | log(f_decrease_pop) |
| frac | pop | pop_log_growth_frac | frac_log_growth_pop | derived | scalar | existing | log(f_growth_pop) |
| frac | patient | patient_log_decrease_frac | frac_log_decrease_patient | derived | vector[n_patient] | existing | log(f_decrease_patient) |
| frac | patient | patient_log_growth_frac | frac_log_growth_patient | derived | vector[n_patient] | existing | log(f_growth_patient) |
| frac | trial | — | frac_sd_trial_intercept | sd | scalar | new | Trial intercept SD |
| frac | trial | — | frac_raw_trial_intercept | raw | vector[n_trial] | new | Trial intercept draws |
| frac | trial | — | frac_effect_trial_intercept | effect | vector[n_trial] | new | Trial intercept effect |
| frac | trial | — | frac_logit_loc_trial | loc | vector[n_trial] | new | Realized trial logit intercept |
| frac | trial | raw_trial_decrease_frac_logit_coef | frac_raw_trial_slope | raw | matrix[n_trial,K_frac] | activate | Will become slope raw draws |
| frac | trial | trial_decrease_frac_logit_coef_sd | frac_sd_trial_slope | sd | vector[K_frac] | activate | Will become slope SD |
| frac | trial | — | frac_effect_trial_slope | effect | matrix[n_trial,K_frac] | new | Derived slope deviations |
| frac | patient | raw_patient_decrease_frac_logit_coef | frac_raw_patient_slope | raw | matrix[n_patient,K_frac] | activate | To be wired into linpred |
| frac | patient | patient_decrease_frac_logit_coef_sd | frac_sd_patient_slope | sd | vector[K_frac] | activate |  |
| frac | patient | — | frac_effect_patient_slope | effect | matrix[n_patient,K_frac] | new | Derived slope deviations |
| frac | patient | — | frac_linpred_trial | derived (linpred) | vector[n_trial] | new | Trial-level aggregate linpred (optional) |

#### 2.2.3 Initial Proportions (init)

Shape notation: scalar, vector[n_patient], vector[n_trial], vector[K_init], matrix[n_trial,K_init], matrix[n_patient,K_init]. (K_init = number of initial proportion covariates.)

| Module | Level | Current | Proposed | Type | Shape | Status | Notes |
|--------|-------|---------|----------|------|-------|--------|-------|
| init | pop | pop_decrease_prop_logis | init_logit_loc_pop | loc | scalar | existing | Logit intercept |
| init | pop | pop_decrease_prop_logis_coef_qr | init_coef_qr_pop | coef_qr | vector[K_init] | existing | QR coefficients |
| init | pop | pop_decrease_prop_logis_coef | init_coef_pop | coef | vector[K_init] | existing | Original-scale betas |
| init | pop | — | init_linpred_pop | derived (linpred) | vector[n_patient] | new | Population covariate contribution |
| init | patient | patient_decrease_prop_logis_linpred | init_linpred_patient | derived (linpred) | vector[n_patient] | existing | Linpred (will exclude intercept) |
| init | trial | trial_decrease_prop_logis_sd | init_sd_trial_intercept | sd | scalar | existing | Trial intercept SD |
| init | trial | raw_trial_decrease_prop_logis | init_raw_trial_intercept | raw | vector[n_trial] | existing | Trial raw |
| init | trial | trial_decrease_prop_logis_effect | init_effect_trial_intercept | effect | vector[n_trial] | existing | Trial effect |
| init | trial | — | init_logit_loc_trial | loc | vector[n_trial] | new | Realized trial logit |
| init | patient | patient_decrease_prop_logis_sd | init_sd_patient_intercept | sd | scalar | existing | Patient intercept SD |
| init | patient | raw_patient_decrease_prop_logis | init_raw_patient_intercept | raw | vector[n_patient] | existing | Patient raw |
| init | patient | patient_decrease_prop_logis_effect | init_effect_patient_intercept | effect | vector[n_patient] | existing | Patient effect |
| init | patient | patient_decrease_prop_logis | init_logit_loc_patient | loc | vector[n_patient] | existing | Realized logit patient |
| init | patient | patient_log_decrease_prop | init_log_decrease_patient | derived | vector[n_patient] | existing | log(p_decrease) |
| init | patient | patient_log_growth_prop | init_log_growth_patient | derived | vector[n_patient] | existing | log(p_growth) |
| init | trial | — | init_sd_trial_slope | sd | vector[K_init] | new | Trial slope SD |
| init | trial | — | init_raw_trial_slope | raw | matrix[n_trial,K_init] | new | Trial slope raw |
| init | trial | — | init_effect_trial_slope | effect | matrix[n_trial,K_init] | new | Trial slope effect |
| init | patient | — | init_sd_patient_slope | sd | vector[K_init] | new | Patient slope SD |
| init | patient | — | init_raw_patient_slope | raw | matrix[n_patient,K_init] | new | Patient slope raw |
| init | patient | — | init_effect_patient_slope | effect | matrix[n_patient,K_init] | new | Patient slope effect |
| init | patient | — | init_linpred_trial | derived (linpred) | vector[n_trial] | new | Optional trial linpred summary |

#### 2.2.4 Derived Naming (Standardized)
For any logit intercept `<m>_logit_loc_<level>` produce:
- `<m>_prob_<level>` = inv_logit(...)
- `<m>_log_<level>` = log(prob)
- `<m>_log_comp_<level>` = log1m(prob)

These may be emitted only if required by downstream computations to avoid clutter.

### 2.3 Flag Specification (Explicit, No Arrays)

All modules have their own explicit flags; no index-based arrays (improves greppability & reduces cognitive load). Intercept vs slope (covariate) gating separated. Population intercepts always on and therefore have no flag.

#### Total Rate (tr) Flags
| Flag | Meaning | Default |
|------|---------|---------|
| enable_pop_cov_tr | Enable population covariate linear model (QR) | 0 |
| enable_trial_intercept_tr | Trial random intercept hierarchy | 0 |
| enable_trial_cov_tr | Trial slope deviations (QR) | 0 |
| enable_patient_intercept_tr | Patient random intercept hierarchy | 0 |
| enable_patient_cov_tr | Patient slope deviations (QR) | 0 |

#### Fraction Mix (frac) Flags
| Flag | Meaning | Default |
|------|---------|---------|
| enable_pop_cov_frac | Population covariate model | 1 |
| enable_trial_intercept_frac | Trial random intercept | 0 |
| enable_trial_cov_frac | Trial slope deviations | 0 |
| enable_patient_intercept_frac | Patient random intercept | 1 (current behavior) |
| enable_patient_cov_frac | Patient slope deviations | 0 (scaffolded) |

#### Initial Proportions (init) Flags
| Flag | Meaning | Default |
|------|---------|---------|
| enable_pop_cov_init | Population covariate model | 1 |
| enable_trial_intercept_init | Trial random intercept | 1 (retains current) |
| enable_trial_cov_init | Trial slope deviations | 0 |
| enable_patient_intercept_init | Patient random intercept | 1 (retains current) |
| enable_patient_cov_init | Patient slope deviations | 0 |

#### Flag Semantics
- If an intercept hierarchy flag is 0: omit (do not declare) its sd/raw/effect parameters.
- If a covariate hierarchy flag is 0: omit slope sd/raw/effect objects for that level; corresponding contribution to `linpred` is skipped.
- `linpred` objects exist only at levels where at least one covariate source (population or hierarchical) is enabled; otherwise they can be elided (implementation choice—may still allocate for simplicity with size 0).

#### Removed Legacy Concepts
- Array-based `enable_*[3]` flags removed (never implemented in code).
- No backward compatibility alias layer (user decision) — names will switch in-place during refactor.
- Legacy flags (`pop_rates_param_only`, `pop_initial_states_param_only`, etc.) will be deleted instead of mapped.

### 2.4 Parameter Allocation Rules

For each module M in {tr, frac, init} and each level L in {trial, patient}:

If `enable_<L>[M] == 1`:
* Include `M_sd_<L>` (scalar or row_vector per-covariate if slopes) with prior.
* Include `M_raw_<L>` of matching dimension with std_normal() prior.
* Define `M_effect_<L> = M_sd_<L> * M_raw_<L>`.
Else: allocate zero length to minimize dimension (Stan optimization) OR simply skip declaration via conditional include.

Slopes: If `enable_<L>_slopes[M] == 1`, declare `M_raw_<L>_slope` and `M_sd_<L>_slope` as row/column vectors per covariate, and add their contribution to `M_linpred_patient` (or trial-level aggregated onto patients by trial index).

### 2.5 Linear Predictor Assembly (Pattern)

```
// Population base
frac_logit_loc_patient = frac_logit_loc_pop;

// Population covariates
frac_linpred_patient = Q * frac_theta_pop;  // always

// Trial slopes (if enabled)
if (enable_trial_slopes[2]) frac_linpred_patient += rows_dot_product(Q, frac_theta_trial_qr[patient_trial]);

// Patient slopes
if (enable_patient_slopes[2]) frac_linpred_patient += rows_dot_product(Q, frac_theta_patient_qr);

// Intercept RE
if (enable_patient_intercept[2]) frac_logit_loc_patient += frac_sd_patient_intercept * frac_raw_patient_intercept;

// Final logit
frac_logit_loc_patient += frac_linpred_patient;
```

Analogous pattern for `init_*` module.

### 2.6 Derived Quantities Standardization

For any logit intercept `x_logit_loc_*`, define:
* Probability: `x_prob_* = inv_logit(x_logit_loc_*)`
* Log probability: `x_log_* = log(x_prob_*)`
* Complement log probability: `x_log_comp_* = log1m_inv_logit(x_logit_loc_*)`

Apply this schema to both fraction and initial state modules, removing hand-written `-log1p_exp()` calls for clarity (can wrap in inline function if performance sensitive).

### 2.7 Orphan Parameter Resolution

Option A (prune): Remove current trial & patient slope latent structures for fraction until a use case emerges.
Option B (activate): Implement hierarchical slope contributions now using the standardized pattern, updating priors and initializers accordingly.

Given current design complexity, recommend Option A initially; document activation steps.

### 2.8 Backward Compatibility
None required. Downstream code must update concurrently with refactor PR. No alias block will be added.

### 2.9 Refactor Execution Phases (Updated)

1. Introduce new names side-by-side (add transformed block computing new variables from legacy ones + aliases) – no sampler dimension change.
2. Remove orphan / unused fraction slope parameters (if pruning) and adjust initializers.
3. Introduce structured enable arrays & map legacy flags in transformed data.
4. Rewrite parameter blocks using new naming; keep legacy alias layer.
5. Drop legacy names after downstream update.

### 2.10 Open Decisions
1. Activate or postpone wiring of slope deviations (trial/patient) for each module? (Default: scaffold but off.)
2. Keep defaults where fraction patient intercept is on but trial intercept off? (Currently specified.)
3. Standardize on omission vs zero-length declarations for disabled pieces (recommended: omit to reduce parameter count). Confirm.

---
Step 2 draft appended (no code changes yet). Awaiting confirmation / adjustments before moving to naming implementation.

## 3. Final Naming Tokens (Summary)
`<module>_<scale?>_<component>_<level?>` with modules `tr|frac|init`; components: `loc|coef_qr|coef|linpred|sd|raw|effect` plus derived forms `log_decrease|log_growth|rate_resid_*`. `linpred` excludes intercept terms.

## 4. Defaults Snapshot
- Total Rate: only `tr_loc_pop` active (all hierarchies & covariates off).  
- Fraction: pop intercept + pop covariates + patient intercept; trial & slope deviations off.  
- Initial: pop intercept + pop covariates + trial & patient intercepts; slope deviations off.

## 5. Next Step
Implement explicit flags in data blocks and begin parameter block re-write using new naming (Task 4 / 5). This document now supersedes earlier duplicated flag/naming sections (removed above).

## 6. Orphan Elimination Guarantee
Scaffolding ensures: no `sd`/`raw` pair or slope deviation is declared unless its module & level flag is enabled. This prevents prior-only (orphan) sampling.

## 7. Open Decisions Recap
See Section 2.10 — please confirm omission vs zero-length approach (recommended: omit) before code changes.

---
Document cleaned: superseded sections removed; no backward compatibility layer; explicit per-module flags adopted.
