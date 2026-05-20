# Multi-level dispersion (SD) hierarchy — design

**Issue:** [#110](https://github.com/azu-oncology-rd/cds-dsi-pioneer-pioneer-2026/issues/110)
**Related:** [#109](https://github.com/azu-oncology-rd/cds-dsi-pioneer-pioneer-2026/issues/109) — posterior geometry pathology
**Branch:** `karim/multilevel-sd` (based on `karim/re-cp`)
**Date:** 2026-04-24
**Scope:** Phase 1 (infrastructure, bit-exact when inactive) + Phase 2 (activate on `psa_standalone`)

---

## 1. Motivation

The current multilevel parameter framework supports a **location hierarchy** with level-indexed intercepts and slopes (modes: FE / RE / RE_CP / RE_GP). At each random-effect level, the raw deviations are scaled by a single scalar SD per level:

```stan
// current (tr/parameters.stan)
array[n_re_levels_tr_intercept] real<lower=0> tr_sd_level_intercept_raw;
vector[n_raw_groups_tr_intercept] tr_raw_level_intercept;

// assembly (tr/transformed_parameters.stan)
tr_scaled_level_intercept[e_lo:e_hi] =
  tr_sd_level_intercept[lv] * tr_raw_level_intercept[r_lo:r_hi];
```

This forces every group at a given level to share a single dispersion. In the pioneer joint trial+RWD fit, 2 trial arms and 5 RWD arms of patients share one `sd_patient_intercept`. Diagnostic work in issue #109 identified this pooling as a major driver of the posterior geometry pathology: sampling-phase E-BFMI 0.08–0.17, maxtree hits 100+ per chain. Even switching patient intercepts to CP (#109's partial fix) leaves the dispersion pooled — one chain of four remained stuck at tiny stepsize.

**Goal.** Extend the framework so the dispersion scaling raw draws at a location level can itself vary across the levels above it, rather than being a single scalar. Each location level L with an active SD sub-hierarchy gets:

- a **population log-SD** `log_sd_L_pop`
- a **sub-hierarchy of deviations** across sub-levels 1..L−1
- a per-group `sd_L[g]` used in the existing location-side assembly

Under the pioneer use case, trial arms learn a tight patient-level SD and RWD arms learn a loose one — no single compromise.

---

## 2. Non-goals

- Crossed (non-nested) level structures. The implementation assumes nested levels (trial ⊇ arm ⊇ patient). Vocabulary is neutral (no "parent"/"child") so a crossed extension later is non-breaking.
- SD sub-hierarchies for slopes (covariate coefficients). Deferred to Phase 4.
- Rollout across tr / frac / init / multistate modules. Deferred to Phase 3. Phase 2 activates on `psa_standalone` only.
- RE_GP mode on the sub-hierarchy side. Sub-levels are restricted to `{NONE, FE, RE, RE_CP}`.
- Changes to behavior when no sub-hierarchy is active. When the mode tibble is empty, draws must be bit-identical to the current branch.

---

## 3. Vocabulary

The implementation stays nested, but the design vocabulary is neutral so a future crossed-levels extension is non-breaking.

| Do NOT use            | Use instead                                        |
|-----------------------|----------------------------------------------------|
| "parent level"        | "sub-level"  or  "level ℓ of the sub-hierarchy"    |
| "parent of"           | (no equivalent — sub-levels aren't owned)          |
| "nests in"            | "partitions"  or  "groups"                         |
| "child level"         | (no equivalent)                                    |
| "finer / coarser"     | "level L has sub-levels 1..L−1" (positional)       |
| "hierarchy" (in prose)| "multilevel structure" / "level stack"             |

Structural indexing words (`level`, `lv`, `n_levels`, `patient_level_groups`) are retained — they are positional, not relational.

---

## 4. Mathematical model

**Notation convention.** This section uses mathematical shorthand (`log_sd_L`, `δ_L_ℓ`, `sd_L[g]`) for the model. Actual Stan variable names follow the existing `<module>_<kind>_level_<parameter>` template — see §7.2 for declarations and the table below:

| Math shorthand            | Stan variable (for `tr` module, intercept parameter)       |
|---------------------------|------------------------------------------------------------|
| `sd_L[g]` (per-group SD)  | `tr_sd_intercept_pergroup[lv][g]`                          |
| `log_sd_L_pop`            | `tr_log_sd_level_intercept_pop[lv]`                        |
| `hyperscale_L_ℓ`          | `tr_sd_hyperscale_level_intercept[lv, sub_lv]`             |
| raw bucket for δ          | `tr_raw_log_sd_level_intercept` (flat; indexed by `raw_pos`) |
| cp bucket for δ           | `tr_cp_log_sd_level_intercept`  (flat; indexed by `cp_pos`)  |

Other modules replace the `tr_` prefix with `frac_`, `init_`, `ms_`, `psa_`. Other SD-scaled parameters replace `intercept` with the parameter name (e.g., `slope` in Phase 4).

For each location level L with an active SD sub-hierarchy, the log-SD applied to group `g` at level L for a unit `i` in that group is:

```
log_sd_L[g]  =  log_sd_L_pop  +  Σ_{ℓ=1}^{L-1}  δ_L_ℓ[ g_ℓ(unit-at-g) ]

sd_L[g]      =  exp( log_sd_L[g] )
```

Each sub-level deviation `δ_L_ℓ` has a **mode** ∈ `{NONE, FE, RE, RE_CP}`:

| Mode  | Sub-level ℓ's contribution                                      | Parameter storage          |
|-------|-----------------------------------------------------------------|----------------------------|
| NONE  | `δ_L_ℓ[g] ≡ 0`. No allocation.                                  | Size-0 raw & CP buckets    |
| FE    | `δ_L_ℓ[g] = fixed_hyperscale_L_ℓ × raw_L_ℓ[g]`, `raw ~ std_normal()` | `raw` bucket, fixed scale |
| RE    | `δ_L_ℓ[g] = free_hyperscale_L_ℓ  × raw_L_ℓ[g]`, `raw ~ std_normal()` | `raw` bucket, free scale  |
| RE_CP | `δ_L_ℓ[g] ~ normal(0, free_hyperscale_L_ℓ)` directly (centered) | `cp` bucket, free scale    |

Mirrors the existing location-hierarchy mode vocabulary verbatim. FE/RE/RE_GP go into the `raw` bucket, RE_CP goes into the `cp` bucket — same split as `split_cp_ncp_pos` on the location side.

**Location-side assembly is unchanged by activation.** The existing mode dispatch (FE / RE / RE_CP) at the location level continues to operate. The only change: scalar `tr_sd_level_intercept[lv]` is replaced by per-group `tr_sd_intercept_pergroup[lv]`.

```stan
// Location mode at L = RE (NCP):  form UNCHANGED, sd lookup changes
tr_scaled_level_intercept[e_lo:e_hi] =
  tr_sd_intercept_pergroup[lv] .* tr_raw_level_intercept[r_lo:r_hi];

// Location mode at L = RE_CP:  form UNCHANGED, prior scale changes
// In priors.stan:
tr_cp_level_intercept[c_lo:c_hi] ~ normal(0, tr_sd_intercept_pergroup[lv]);

// Location mode at L = FE:  tr_sd_intercept_pergroup[lv] is data-filled; otherwise identical
```

When the sub-hierarchy is inactive for level L, `tr_sd_intercept_pergroup[lv][g]` is a constant vector filled with the existing scalar `tr_sd_level_intercept[lv]` — so the scaled assembly is bit-identical to the current code.

---

## 5. Semantic validity (nesting constraint)

Because the implementation is nested, only sub-levels positionally *before* the location level are valid — sub-level ℓ must satisfy `ℓ < L` in the level stack order. In SCLC (trial=1, arm=2, patient=3):

| Location level | Valid sub-levels        |
|----------------|-------------------------|
| trial (L=1)    | (none — no sub-levels)  |
| arm (L=2)      | {trial}                 |
| patient (L=3)  | {trial, arm}            |

This is a **structural constraint imposed by the nested implementation**. It is not expressed as "coarser than" in the design vocabulary; it is simply that sub-levels are positionally 1..L−1 in the level stack. A future crossed-levels extension would replace this positional rule with a per-pair partition-constancy check, but the config shape and Stan internals would stay the same.

Invalid configurations (e.g., `(location_level='trial', sub_level='arm', mode='re')`) are rejected by R-side validation pre-launch.

---

## 6. R-side configuration

### 6.1 Per-module mode tibble

Each module with an SD-scaled parameter gets a new list-column in the variant tribble. The column value is a flat long-form tibble:

```r
# targets/pioneer_targets.R
psa_sd_intercept_modes = tibble::tribble(
  ~location_level, ~sub_level, ~mode,
  "patient",       "arm",      "re_cp"
)
# Unmentioned (location_level, sub_level) pairs default to "none".
# An empty tibble ≡ all NONE ≡ current scalar-SD behavior (bit-exact).
```

Columns:
- `location_level`: name of the location level (must be in the active level stack)
- `sub_level`: name of a sub-level (must be positionally before `location_level`)
- `mode`: one of `"none"`, `"fe"`, `"re"`, `"re_cp"` (the mode vocabulary from the location hierarchy; `"gp"` is rejected)

One such column per module per SD-scaled parameter:
- `tr_sd_intercept_modes`
- `frac_sd_intercept_modes`
- `init_sd_intercept_modes`
- `ms_sd_intercept_modes` (one row per transition via `location_level` namespacing, or a nested list — TBD in Phase 3)
- `psa_sd_intercept_modes` (Phase 2 target)

### 6.2 Validation helper

Added to the existing `r/multi_level_hierarchy.R` (alongside `create_hierarchy_structure()`, which is the canonical home for level-stack helpers):

```r
validate_sd_modes <- function(sd_modes, level_stack, allowed_modes = c("none", "fe", "re", "re_cp"))
```

Checks, for each row:
1. `location_level` ∈ `level_stack`
2. `sub_level` ∈ `level_stack`
3. `position(sub_level) < position(location_level)` (positional rule)
4. `mode` ∈ `allowed_modes` (rejects `"gp"`)

Called in the module-specific `prepare_*_stan_data()` step. Fails fast via `stop()` with a message naming the offending row. Empty tibble passes trivially.

### 6.3 Prior columns

Per-module, per-SD-parameter columns for sub-hierarchy priors (used only when the mode tibble is non-empty):

- `<module>_log_sd_<param>_pop_mean`
- `<module>_log_sd_<param>_pop_sd`
- `<module>_sd_<param>_hyperscale_sd` (half-normal prior scale for RE/RE_CP hyperscales)

Prior defaults (in `r/priors.R`): chosen so that, under an active sub-hierarchy with small hyperscales, the marginal prior on `tr_sd_intercept_pergroup[lv]` approximates the current `half_normal(0, σ_current)` prior on `tr_sd_level_intercept[lv]`. (When the mode tibble is empty these parameters are not allocated at all — bit-exactness does not depend on calibration.) Specifically:
- `<module>_log_sd_<param>_pop_mean = log(median of current half-normal prior)`
- `<module>_log_sd_<param>_pop_sd ≈ 0.5` (roughly matches the tail width of the current half-normal on log-scale)
- `<module>_sd_<param>_hyperscale_sd ≈ 0.3`; hyperscales are only active when the user opts in, so the prior can be mildly regularizing

These defaults are not exercised when the mode tibble is empty (the new parameters aren't allocated), so calibration choices are non-blocking for Phase 1.

### 6.4 Initializers

`r/initializers.R` extended to route initial draws for the new parameters:
- `<module>_log_sd_level_<param>_pop` → `log(current scalar sd default)`
- `<module>_raw_log_sd_level_<param>` → zero
- `<module>_cp_log_sd_level_<param>` → zero
- `<module>_sd_hyperscale_level_<param>` → small positive (e.g., 0.1)

Zero-sized when all modes are NONE.

---

## 7. Stan implementation

### 7.1 New helper in `stan/hierarchy.stanfunctions`

Parallel to `split_cp_ncp_pos`, operating on a 2D mode matrix:

```stan
/**
 * Split sub-hierarchy (L, ℓ) pairs into NCP (raw) and CP (cp) buckets by mode.
 * For each location level L and sub-level ℓ ∈ 1..L-1:
 *   - modes in {FE=1, RE=2} → raw bucket
 *   - mode == RE_CP (4)      → cp bucket
 *   - mode == NONE (0)       → skipped
 * RE_GP (3) is rejected via fatal_error.
 *
 * @param n_levels          number of levels in the stack
 * @param n_groups_per_level  group counts per level
 * @param sd_mode           [n_levels, n_levels] matrix; entries with ℓ >= L are 0
 *
 * @return tuple(raw_total, raw_pos[L][n_levels+1], cp_total, cp_pos[L][n_levels+1])
 */
tuple(int, array[,] int, int, array[,] int) split_sd_cp_ncp_pos(
  int n_levels,
  array[] int n_groups_per_level,
  array[,] int sd_mode
);
```

Emits per-location-level `raw_pos[L][ℓ]` and `cp_pos[L][ℓ]` position arrays: `raw_pos[L][ℓ+1] - raw_pos[L][ℓ]` is the number of groups at sub-level ℓ for location level L's sub-hierarchy (in the raw bucket).

### 7.2 New declarations in each module

**`<module>/flags.stan`** additions:

```stan
// 2D mode matrix: [location_level L, sub-level ℓ]
array[n_levels, n_levels] int enable_sd_level_intercept_mode_<module>;
// Derived boolean: has_sd_subhierarchy_<module>_intercept[lv]
// computed in transformed_data as (sum(enable_sd_level_intercept_mode_<module>[lv, :]) > 0)
```

**`<module>/parameters.stan`** additions (sizes zero when all modes NONE):

```stan
// Population log-SD per location level (allocated when that level's row in the mode matrix has any non-NONE entry):
vector[n_levels_with_subhier_<module>_intercept] <module>_log_sd_level_intercept_pop;

// Hyperscales for RE / RE_CP sub-level deviations:
vector<lower=0>[n_re_hyperscales_<module>_intercept] <module>_sd_hyperscale_level_intercept_raw;

// Flat raw bucket across all (L, ℓ) with mode ∈ {FE, RE}:
vector[n_raw_groups_<module>_log_sd_intercept] <module>_raw_log_sd_level_intercept;

// Flat cp bucket across all (L, ℓ) with mode == RE_CP:
vector[n_cp_groups_<module>_log_sd_intercept]  <module>_cp_log_sd_level_intercept;
```

Size-of scalars (`n_levels_with_subhier_*`, `n_re_hyperscales_*`, `n_raw_groups_*_log_sd_*`, `n_cp_groups_*_log_sd_*`) are computed in `transformed_data` from the mode matrix and emitted by `split_sd_cp_ncp_pos`.

### 7.3 `transformed_parameters.stan` structure

Changes to the existing location-side code are minimal. For each SD-scaled location parameter:

**Step A — build per-group SD vector `tr_sd_intercept_pergroup[lv]`.**

```stan
// For EACH location level lv, allocate a per-group SD vector.
// When sub-hierarchy inactive for lv: fill with scalar sd (current behavior).
// When active: sum log-deviations, then exp.
array[n_levels] vector[n_forecast_groups_per_level[lv]] tr_sd_intercept_pergroup;

for (lv in 1:n_levels) {
  int ngroups_lv = n_forecast_groups_per_level[lv];
  if (!has_sd_subhierarchy_tr_intercept[lv]) {
    // Inactive: constant fill with existing scalar.
    tr_sd_intercept_pergroup[lv] = rep_vector(tr_sd_level_intercept[lv], ngroups_lv);
  } else {
    // Active: build log_sd[g] additively, then exp.
    vector[ngroups_lv] log_sd_g = rep_vector(tr_log_sd_level_intercept_pop[<idx>], ngroups_lv);
    for (sub_lv in 1:(lv - 1)) {
      int mode = enable_sd_level_intercept_mode_tr[lv, sub_lv];
      if (mode == LEVEL_MODE_NONE) continue;
      // deviation[g_ℓ] gathered from raw/cp bucket, scaled/sampled per mode:
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_pos_tr_log_sd_intercept[lv][sub_lv];
        int c_hi = cp_pos_tr_log_sd_intercept[lv][sub_lv + 1] - 1;
        // deviation per sub-level group → map to location-level groups via sub_level_groups_at_level:
        log_sd_g += tr_cp_log_sd_level_intercept[c_lo:c_hi][ subgroup_idx_L_to_ell[lv, sub_lv] ];
      } else {
        // FE or RE
        int r_lo = raw_pos_tr_log_sd_intercept[lv][sub_lv];
        int r_hi = raw_pos_tr_log_sd_intercept[lv][sub_lv + 1] - 1;
        real hyperscale_lv_ell = (mode == LEVEL_MODE_FE)
          ? tr_fe_sd_hyperscale_intercept[lv, sub_lv]
          : tr_sd_hyperscale_level_intercept_raw[<idx>];
        log_sd_g += hyperscale_lv_ell * tr_raw_log_sd_level_intercept[r_lo:r_hi][ subgroup_idx_L_to_ell[lv, sub_lv] ];
      }
    }
    tr_sd_intercept_pergroup[lv] = exp(log_sd_g);
  }
}
```

`subgroup_idx_L_to_ell[lv, sub_lv]` is a flat index, computed in `transformed_data`, mapping each location-level-lv group to its sub-level-ℓ group — derived from `patient_level_groups` (one representative patient per location-lv group gives the ℓ-group it belongs to; well-defined because implementation is nested).

**Step B — use `tr_sd_intercept_pergroup[lv]` in the existing location assembly.**

Existing location assembly, touched minimally:

```stan
// BEFORE:
tr_scaled_level_intercept[e_lo:e_hi] =
  tr_sd_level_intercept[lv] * tr_raw_level_intercept[r_lo:r_hi];

// AFTER — scalar replaced by per-group lookup (still raw × sd, just vector × vector):
tr_scaled_level_intercept[e_lo:e_hi] =
  tr_sd_intercept_pergroup[lv] .* tr_raw_level_intercept[r_lo:r_hi];

// RE_CP branch prior (in priors.stan):
// BEFORE:  tr_cp_level_intercept[c_lo:c_hi] ~ normal(0, tr_sd_level_intercept[lv]);
// AFTER:   tr_cp_level_intercept[c_lo:c_hi] ~ normal(0, tr_sd_intercept_pergroup[lv]);
```

Note: `tr_raw_level_intercept[r_lo:r_hi]` has length `n_forecast_groups_per_level[lv]` (for RE levels), matching `tr_sd_intercept_pergroup[lv]` for the `.*` product to type-check.

**Blast radius on existing files:** exactly two substitutions per SD-scaled parameter in `transformed_parameters.stan`, plus one in the corresponding `priors.stan` RE_CP branch. No other edits to the existing location-side code.

### 7.4 New priors

**`<module>/priors.stan`** additions (fire only when sub-hierarchy is active):

```stan
// Population log-SD:
tr_log_sd_level_intercept_pop ~ normal(log_sd_pop_mean, log_sd_pop_sd);

// Hyperscales (half-normal via <lower=0> constraint):
tr_sd_hyperscale_level_intercept_raw ~ normal(0, hyperscale_prior_sd);

// NCP raw deviations (std_normal or student_t if enable_student_t_hierarchy):
tr_raw_log_sd_level_intercept ~ std_normal();

// CP deviations — sampled at their hyperscale directly:
// Handled within the per-level loop in priors.stan using the cp_pos array:
for (lv in 1:n_levels) {
  if (!has_sd_subhierarchy_tr_intercept[lv]) continue;
  for (sub_lv in 1:(lv - 1)) {
    if (enable_sd_level_intercept_mode_tr[lv, sub_lv] == LEVEL_MODE_RE_CP) {
      int c_lo = cp_pos_tr_log_sd_intercept[lv][sub_lv];
      int c_hi = cp_pos_tr_log_sd_intercept[lv][sub_lv + 1] - 1;
      real hs = tr_sd_hyperscale_level_intercept_raw[<idx>];
      tr_cp_log_sd_level_intercept[c_lo:c_hi] ~ normal(0, hs);
    }
  }
}
```

---

## 8. Bit-exactness guarantee

When every variant's mode tibble is empty (all-NONE):

1. `enable_sd_level_intercept_mode_tr[L, ℓ] = 0` for all `(L, ℓ)` → `has_sd_subhierarchy_tr_intercept[lv] = 0` for all `lv`.
2. In `transformed_parameters`, the `if (!has_sd_subhierarchy_tr_intercept[lv])` branch fires for every level → `tr_sd_intercept_pergroup[lv]` is a constant fill of `tr_sd_level_intercept[lv]`.
3. `tr_sd_intercept_pergroup[lv] .* tr_raw_level_intercept[r_lo:r_hi]` is mathematically identical to `tr_sd_level_intercept[lv] * tr_raw_level_intercept[r_lo:r_hi]` (scalar broadcast × vector ≡ constant-fill vector ⊙ vector).
4. New parameter blocks (`_raw_log_sd_`, `_cp_log_sd_`, `_log_sd_pop`, `_sd_hyperscale_raw`) have size zero → not sampled, no RNG consumption shift, no prior contributions.
5. No structural change to existing `transformed_data` pre-computation (the new position arrays are emitted alongside existing ones; when all modes NONE, they are all zeros/empty).

Draws from identical seeds are bit-identical to the pre-change branch.

---

## 9. Phase 2: `psa_standalone` activation

New pioneer variant: `psa_standalone_sd_arm_intercept` with:

```r
psa_sd_intercept_modes = tibble::tribble(
  ~location_level, ~sub_level, ~mode,
  "patient",       "arm",      "re_cp"
)
```

This gives each arm its own patient-level log-SD (sampled via RE_CP at the arm sub-level), on top of a population log-SD. All other sub-levels are NONE.

**Diagnostic targets** (posterior fit only):
- Sampling-phase E-BFMI > 0.3 across all chains (currently 0.08–0.17 in combined fit)
- Cumulative max-treedepth hits per chain < 20 (currently 100+)
- Post-warmup stepsize stabilizes in 0.01–0.05
- Posterior on per-arm `exp(log_sd_patient)` separates trial arms from RWD arms by ≥ 1 log-unit (~2.7× SD ratio)

Prior-only fit: passes as a sanity check — posterior should not be degenerate or push against priors.

A paired baseline run on `psa_standalone_combined` (the existing variant) uses the same random seed for convergence-delta comparison.

---

## 10. Testing

### 10.1 Bit-exact regression (the primary gate)

A single test script (likely `tests/testthat/test-bit-exact-sd-hierarchy.R`) runs a known pioneer fit target (e.g., `pioneer_fit_posterior_combined_markov_fe_psa_re_ms`) with empty mode tibbles everywhere, and diffs the Stan output CSVs against a frozen reference from `main`. Test fails on any byte-level difference.

### 10.2 Helper unit tests

`tests/testthat/test-split-sd-cp-ncp-pos.R`:
- All-zero mode matrix → all-zero position arrays, zero totals
- Single RE_CP entry → cp bucket sized correctly, raw bucket size zero
- Mixed FE / RE / RE_CP → raw and cp buckets sized correctly per (L, ℓ)
- RE_GP entry → `fatal_error`

`tests/testthat/test-validate-sd-modes.R`:
- Empty tibble → passes
- Valid rows → passes
- Unknown `location_level` / `sub_level` → errors with named row
- Positional rule violation (`sub_level` at or after `location_level`) → errors
- Disallowed mode (`"gp"` or typo) → errors

### 10.3 Functional test

`psa_standalone_sd_arm_intercept` fit completes without errors and meets the Phase 9 convergence targets. Posterior summary for per-arm `log_sd_patient` reviewed by eye.

### 10.4 Stan syntax

`stanc --include-paths=stan,stan/psa stan/psa/pioneer.stan` must parse the updated module cleanly.

---

## 11. Implementation plan summary

1. **Helper + type scaffolding** — `split_sd_cp_ncp_pos`, `LEVEL_MODE_*` constant reuse, position-array conventions.
2. **PSA module plumbing** — add flags, parameters, transformed_parameters, priors blocks for the sub-hierarchy in `stan/modules/psa/` (or wherever the psa module lives).
3. **R side** — `validate_sd_modes()`, add `psa_sd_intercept_modes` list-column to pioneer tribble, initializer routing.
4. **Bit-exact test passes** on empty mode tibbles for an existing pioneer variant.
5. **Phase 2 variant** — add `psa_standalone_sd_arm_intercept` to tribble with the arm sub-level configured; launch a Domino job; verify convergence.

Each step produces a commit that builds and passes the existing test suite.

---

## 12. Open questions

These are pragmatic choices that can be resolved during implementation without changing the spec:

1. **Prior calibration for `<module>_log_sd_level_<param>_pop`** — the lognormal-to-half-normal matching is approximate. Prior-predictive checks in Phase 2 will confirm the defaults are sensible; adjust as needed.
2. **Naming convention for multi-transition modules** (`ms_*`) — how to namespace per-transition SD sub-hierarchies in the mode tibble (e.g., `location_level = "patient_01"` vs. nested list column). Only becomes relevant in Phase 3.
3. **Initializer values for RE_CP cp buckets** — zero is safe; optionally small-random may help escape degenerate starts. Not a spec concern.

---

## 13. References

- Issue #110 — this spec's originating issue
- Issue #109 — posterior geometry diagnostic experiments (jobs #658, #660, #663, #666, #670, #673)
- `docs/multi_level_hierarchy_design.md` — current location-hierarchy design
- `stan/hierarchy.stanfunctions` — existing helpers (`compute_level_module_flags`, `split_cp_ncp_pos`)
- `karim/re-cp` branch — RE_CP mode for the location hierarchy; SD sub-hierarchy closely mirrors that work
