# Laplace Surrogate — Phase 2: Production Integration

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Integrate the *validated* quadratic surrogate (Phase 1 PASSED — see
`2026-06-08-laplace-surrogate-backgrounded-trials.md` §5a) into
`stan/tumor/sf-ssls-lfo.stan` as a **modular, runtime-toggleable** feature, and
**remove** the old dead hand-coded `stan/modules/laplace/`.

**Validated design (source of truth = `stan/experiments/laplace_surrogate_mixed.stan`):**
- Quadratic-in-time surrogate, **intercept pinned** at `β_pop[1]`, marginalize
  `(β₁,β₂)` only → `hessian_block_size = 2`.
- **GH-3 pushforward bridge** `surrogate_bridge(...)` returns BOTH the population
  mean `β_pop = E[β(θ)]` (kills the Jensen-gap location bias) AND the `(β₁,β₂)`
  covariance `Σ_β = Cov[β(θ)]`, from one 27-point quadrature.
- Gate on the **PFS endpoint**; PASS confirmed (median PFS identical, ≤1 wk
  quantile diff, ≤0.019 event-rate diff).

**Three requirements driving this phase (user, 2026-06-09):**
1. **Modular** — self-contained `stan/modules/laplace_surrogate/`; removing the
   module + its `#include`s + one guarded block fully removes the feature.
2. **Runtime switch** — `enable_background_surrogate` data flag (default 0); flip
   without recompiling. When 0, model behaves exactly as today.
3. **Replace old laplace completely** — delete `stan/modules/laplace/` (5 files)
   + dead `r/pioneer/prepare_laplace_data.R`. Confirmed unused (no `#include`,
   no targets reference).

**Production reality the toy abstracted (from code-explorer trace):**
- Population means in the real model: `tr_loc_pop`, **`frac_logit_loc_pop`**,
  **`init_logit_loc_pop`** (note `_loc_` infix for frac/init). All in `model{}`
  scope (parameters).
- The toy's scalar `tr_sd`/`frac_sd`/`init_sd` becomes the **patient-level total
  marginal SD** = `sqrt(Σ_lv *_sd_level_intercept[lv]²)` over enabled levels.
  `*_sd_level_intercept` is an `array[n_levels] real` transformed parameter, in
  `model{}` scope, in all three modules (`tr/frac/init/transformed_parameters.stan:13`).
- **TR-only SD sub-hierarchy** makes the TR per-patient SD vary by patient when
  active. Detect via existing scalar `n_subhier_active_tr_intercept`
  (`tr/transformed_data.stan:54`). **Decision: guard it out** — `fatal_error` if
  `enable_background_surrogate=1 && n_subhier_active_tr_intercept>0`. frac/init
  never have this.
- Per-patient **LOD offset**: the real likelihood uses
  `log_lod - log_baseline_sld[i]` (`sf-ssls-lfo.stan:98`) because burden is
  normalized to each patient's baseline. The surrogate must apply the same
  per-patient offset, NOT a single global `log_lod`.
- Measurement noise param is **`measure_sd_sld`**.
- Visit week index per row: **`t_patient_visit_idx`** (used at line 140); the
  surrogate's `t = t_patient_visit_idx[v] - 1`.

**Tech Stack:** Stan (CmdStan 2.39), CmdStanR, R (tidyverse), `stan/modules/` includes.

---

## Conventions (read before starting)

- Native pipe `|>`; tidyverse; no `install.packages()`; no hardcoded subject IDs;
  no `<<-`; explicit `store=` in targets.
- Stan: built-in zero constructors; `fatal_error` (not `reject`) for data-validation
  in transformed data; NCP for hierarchical params.
- CmdStan path for any standalone check: `~/.cmdstan/cmdstan-2.39.0`.
- Stan syntax check:
  `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor <file>`

---

## Task 1: Remove the old dead hand-coded laplace (separate cleanup commit)

**Files:**
- Delete: `stan/modules/laplace/{flags,data,transformed_data,model}.stan`,
  `stan/modules/laplace/laplace.stanfunctions`
- Delete: `r/pioneer/prepare_laplace_data.R`

- [ ] **Step 1: Re-confirm nothing references them** (safety before deletion)

Run:
```bash
grep -rn "modules/laplace/" stan/ --include=*.stan | grep -v "modules/laplace_surrogate"
grep -rln "prepare_laplace_data\|prepare_laplace" r/ targets/ --include=*.R
```
Expected: NO output from either (the surrogate module is `laplace_surrogate`, a
different path). If anything prints, STOP and report — do not delete.

- [ ] **Step 2: Delete the files**

Run:
```bash
git rm stan/modules/laplace/flags.stan stan/modules/laplace/data.stan \
       stan/modules/laplace/transformed_data.stan stan/modules/laplace/model.stan \
       stan/modules/laplace/laplace.stanfunctions \
       r/pioneer/prepare_laplace_data.R
rmdir stan/modules/laplace 2>/dev/null || true
```

- [ ] **Step 3: Verify the tumor model still compiles** (it never used these)

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: clean parse (no output).

- [ ] **Step 4: Commit**

```bash
git commit -m "chore(laplace): remove dead hand-coded laplace module

The hand-coded per-patient Newton solver (stan/modules/laplace/) and
r/pioneer/prepare_laplace_data.R are unused — no model #includes them and
no targets pipeline references the R script. Superseded by the validated
built-in-Laplace surrogate (stan/modules/laplace_surrogate/, added next).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Create the modular `laplace_surrogate` module

Follow the standard module file pattern so the feature is self-contained and
removable. The surrogate functions are promoted **verbatim** from the validated
`stan/experiments/laplace_surrogate_mixed.stan` (do NOT re-derive — that file is
the validated artifact).

**Files (create all under `stan/modules/laplace_surrogate/`):**
- `flags.stan` — the runtime toggle
- `data.stan` — anchor times
- `surrogate.stanfunctions` — the 5 validated functions
- `transformed_data.stan` — constant Vandermonde + GH nodes + Laplace knobs + guard
- `likelihood.stan` — the guarded backgrounded-patient `laplace_marginal_tol` block

- [ ] **Step 1: `flags.stan`**

```stan
// Runtime switch for marginalizing backgrounded-trial patients via the
// log-concave quadratic surrogate. 0 = off (backgrounded patients contribute
// nothing to the LFO likelihood, i.e. current behavior); 1 = on.
int<lower=0, upper=1> enable_background_surrogate;
```

- [ ] **Step 2: `data.stan`**

```stan
// Fixed calendar anchors (weeks) for the surrogate quadratic. The first anchor
// MUST be 0 (baseline normalization pins g(0)=0; the surrogate intercept relies
// on it). Tune the other two to the backgrounded trials' observation window.
vector[3] surrogate_anchor_times;
```

- [ ] **Step 3: `surrogate.stanfunctions`** — copy the five functions VERBATIM from
`stan/experiments/laplace_surrogate_mixed.stan` lines 20–126:
`bi_exp_log_burden`, `surrogate_anchor_betas`, `surrogate_bridge` (tuple
mean+cov), `surrogate_ll`, `surrogate_K_fn`. Prepend a module docstring
summarizing the validated design (block-size 2, pinned intercept, GH-3 bridge).
Do NOT alter the bodies.

- [ ] **Step 4: `transformed_data.stan`**

```stan
// --- Laplace surrogate: constants + validity guard --------------------------
// Constant inverse Vandermonde for the 3 fixed anchors (first anchor must be 0).
matrix[3, 3] surrogate_V;
for (k in 1:3) {
  surrogate_V[k, 1] = 1.0;
  surrogate_V[k, 2] = surrogate_anchor_times[k];
  surrogate_V[k, 3] = surrogate_anchor_times[k] * surrogate_anchor_times[k];
}
matrix[3, 3] surrogate_Vinv = inverse(surrogate_V);

// 3-point Gauss-Hermite nodes/weights for the standard normal (weights sum 1).
vector[3] surrogate_gh_x = [-sqrt(3.0), 0.0, sqrt(3.0)]';
vector[3] surrogate_gh_w = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]';

// Laplace inner-solver knobs (linear-Gaussian surrogate => solver 1; may fall
// back to solver 2 on the collinear (t,t^2) basis — harmless, validated).
real surrogate_tolerance = 1e-8;
int surrogate_max_num_steps = 100;
int surrogate_hessian_block_size = 2;   // marginalize (beta_1, beta_2)
int surrogate_solver = 1;
int surrogate_max_steps_line_search = 0;
int surrogate_allow_fallback = 1;
real surrogate_jitter = 1e-10;

// GUARD: the bridge SD assumes the simple summed-variance patient SD, which
// holds for frac/init always and for tr UNLESS the tr SD sub-hierarchy is
// active (then the per-patient SD varies and the single bridge is wrong).
if (enable_background_surrogate == 1 && n_subhier_active_tr_intercept > 0)
  fatal_error("enable_background_surrogate=1 is incompatible with an active tr ",
              "SD sub-hierarchy (n_subhier_active_tr_intercept=",
              n_subhier_active_tr_intercept, "): the surrogate bridge assumes a ",
              "single patient-level tr SD. Disable the surrogate or the tr SD ",
              "sub-hierarchy.");

// First anchor must be 0 (the pinned intercept depends on g(0)=0).
if (enable_background_surrogate == 1 && surrogate_anchor_times[1] != 0.0)
  fatal_error("surrogate_anchor_times[1] must be 0 (baseline); got ",
              surrogate_anchor_times[1]);

vector[enable_background_surrogate == 1 ? n_background_patients * 2 : 0]
  surrogate_theta_0 = rep_vector(0.0,
    enable_background_surrogate == 1 ? n_background_patients * 2 : 0);
```

- [ ] **Step 5: `likelihood.stan`** — the guarded backgrounded-patient term. Note
the per-patient LOD offset and the patient-level total SD via summed variance.

```stan
// --- Laplace surrogate: backgrounded-patient marginalized SLD ---------------
// Backgrounded patients inform population params only; their per-patient latents
// are integrated out via the validated quadratic surrogate. Forecast patients
// (and the current default behavior) are untouched.
if (enable_background_surrogate == 1 && n_background_patients > 0) {
  // Patient-level total marginal SD = sqrt(sum over enabled levels of sd^2).
  // (frac/init always; tr guaranteed simple here by the transformed_data guard.)
  real tr_patient_sd = 0;
  real frac_patient_sd = 0;
  real init_patient_sd = 0;
  for (lv in 1:n_levels) {
    if (enable_level_intercept_tr[lv]   != 0) tr_patient_sd   += square(tr_sd_level_intercept[lv]);
    if (enable_level_intercept_frac[lv] != 0) frac_patient_sd += square(frac_sd_level_intercept[lv]);
    if (enable_level_intercept_init[lv] != 0) init_patient_sd += square(init_sd_level_intercept[lv]);
  }
  tr_patient_sd   = sqrt(tr_patient_sd);
  frac_patient_sd = sqrt(frac_patient_sd);
  init_patient_sd = sqrt(init_patient_sd);

  // Compact background-only views; per-patient LOD offset (burden normalized to
  // each patient's baseline => LOD shifts by log_baseline_sld[p]).
  array[n_background_patients + 1] int bg_pos;
  bg_pos[1] = 1;
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    bg_pos[j + 1] = bg_pos[j] + (ve - vs + 1);
  }
  int n_bg_visits = bg_pos[n_background_patients + 1] - 1;
  vector[n_bg_visits] bg_obs;
  array[n_bg_visits] int bg_time;
  // Per-patient LOD offset replicated per visit (surrogate_ll uses one global
  // log_lod, so we pass each patient's already-offset value via a vector form).
  vector[n_bg_visits] bg_log_lod;
  {
    int w = 1;
    for (j in 1:n_background_patients) {
      int p = background_patient_idx[j];
      int vs, ve;
      (vs, ve) = get_pos(patient_visit_pos, p);
      for (v in vs:ve) {
        bg_obs[w]     = normalized_sld[v];
        bg_time[w]    = t_patient_visit_idx[v];
        bg_log_lod[w] = log_lod - log_baseline_sld[p];
        w += 1;
      }
    }
  }

  vector[3] bg_beta_pop;
  matrix[2, 2] bg_Sigma;
  (bg_beta_pop, bg_Sigma) = surrogate_bridge(
    tr_loc_pop, frac_logit_loc_pop, init_logit_loc_pop,
    tr_patient_sd, frac_patient_sd, init_patient_sd,
    surrogate_Vinv, surrogate_anchor_times,
    surrogate_gh_x, surrogate_gh_w, surrogate_jitter);

  target += laplace_marginal_tol(
    surrogate_ll,
    (bg_beta_pop, measure_sd_sld, bg_log_lod, n_background_patients,
     bg_obs, bg_pos, bg_time),
    surrogate_hessian_block_size,
    surrogate_K_fn,
    (bg_Sigma, n_background_patients),
    (surrogate_theta_0, surrogate_tolerance, surrogate_max_num_steps,
     surrogate_solver, surrogate_max_steps_line_search, surrogate_allow_fallback)
  );
}
```

**IMPORTANT — `surrogate_ll` signature change:** the experiment passed a scalar
`log_lod`; here we need a **per-visit** `bg_log_lod` vector. Update
`surrogate.stanfunctions`'s `surrogate_ll` to take `data vector log_lod_per_visit`
instead of `real log_lod`, and use `log_lod_per_visit[v]` in the `normal_lcdf`
branch. Make the SAME change in `laplace_surrogate_mixed.stan` so the experiment
stays a faithful mirror: update `surrogate_ll`'s signature AND its `model{}`
call site there — replace the scalar `log_lod` argument with a constant vector
`rep_vector(log_lod, n_bg_visits)` built in that model's transformed data. The
toy used a constant LOD, so a constant vector reproduces the validated PASS
exactly (no re-run needed; signature-only change).

- [ ] **Step 6: Syntax-check the module via a throwaway includer**

```bash
cat > $CLAUDE_JOB_DIR/tmp/_inc.stan <<'EOF'
functions {
#include "modules/laplace_surrogate/surrogate.stanfunctions"
}
EOF
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor $CLAUDE_JOB_DIR/tmp/_inc.stan
```
Expected: clean parse.

- [ ] **Step 7: Commit**

```bash
git add stan/modules/laplace_surrogate/ stan/experiments/laplace_surrogate_mixed.stan
git commit -m "feat(laplace): modular laplace_surrogate module (flag-gated)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Wire the module into `sf-ssls-lfo.stan`

**Files:**
- Modify: `stan/tumor/sf-ssls-lfo.stan`

- [ ] **Step 1: Add the functions include** (after line 12,
`#include "modules/tumor/tumor.stanfunctions"`):

```stan
  #include "modules/laplace_surrogate/surrogate.stanfunctions"
```

- [ ] **Step 2: Add the data includes** in the `data {}` block (after line 35,
`#include "modules/state_space/lfo_data.stan"`):

```stan
  #include "modules/laplace_surrogate/flags.stan"
  #include "modules/laplace_surrogate/data.stan"
```

- [ ] **Step 3: Add the transformed-data include** in `transformed data {}` (after
line 57, `#include "_lfo_transformed_data.stan"`). It must come AFTER the tr
module's transformed_data (which defines `n_subhier_active_tr_intercept`) and the
forecast routing (`n_background_patients`) — both already included earlier in the
block, so appending at the end is safe:

```stan
  #include "modules/laplace_surrogate/transformed_data.stan"
```

- [ ] **Step 4: Add the likelihood include** inside `model{}`'s
`if (fit_tumor_data) { ... }`, AFTER the multistate block (after line 109's
closing brace for the `enable_ms_01` block, before the `fit_tumor_data` block's
closing brace at line 110). The surrogate functions and all referenced symbols
(`tr_loc_pop`, `frac_logit_loc_pop`, `init_logit_loc_pop`, `*_sd_level_intercept`,
`measure_sd_sld`, `log_lod`, `log_baseline_sld`, `normalized_sld`,
`t_patient_visit_idx`, `background_patient_idx`, `n_background_patients`,
`patient_visit_pos`, `n_levels`, `enable_level_intercept_*`) are in scope there:

```stan
    #include "modules/laplace_surrogate/likelihood.stan"
```

- [ ] **Step 5: Syntax-check the full model**

```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: clean parse. If a symbol-not-in-scope error appears, fix by confirming
the exact name with `grep` in the relevant module's `transformed_parameters.stan`
/ `parameters.stan` and correcting the reference in `likelihood.stan`.

- [ ] **Step 6: Commit**

```bash
git add stan/tumor/sf-ssls-lfo.stan
git commit -m "feat(laplace): wire flag-gated surrogate into sf-ssls-lfo

Backgrounded patients contribute a solver-1 Laplace-marginalized SLD term
when enable_background_surrogate=1; default 0 preserves current behavior.
Forecast path untouched.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: R pipeline defaults (flag off, anchors)

**Files:**
- Modify: `r/priors.R`
- Modify: the tumor Stan-data assembly (locate in Step 1)

- [ ] **Step 1: Locate the tumor stan-data assembly + prior constructor**

```bash
grep -rn "log_lod\|fit_tumor_data\|measure_sd_sld_alpha" r/ --include=*.R | grep -iv test | head
grep -n "measure_sd_sld_alpha\|tr_loc_pop_mean" r/priors.R
```

- [ ] **Step 2: Add defaults in `r/priors.R`** (near the tumor SLD entries):

```r
    # Backgrounded-trial Laplace surrogate (OFF by default). First anchor must be
    # 0 (baseline); tune the other two to the historical trials' visit window.
    enable_background_surrogate = 0L,
    surrogate_anchor_times = c(0, 12, 28),
```

- [ ] **Step 3: Ensure both fields reach the Stan data list.** If the data list
spreads the prior list, Step 2 suffices; else add the two fields explicitly.
Verify:
```bash
grep -rn "surrogate_anchor_times\|enable_background_surrogate" r/
```
Expected: present in prior defaults and (transitively or explicitly) the data list.

- [ ] **Step 4: Syntax-check**

```bash
Rscript -e 'invisible(parse("r/priors.R"))'
```
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add r/priors.R
git commit -m "feat(laplace): surrogate defaults (off; anchors 0/12/28) in priors

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: End-to-end compile (default OFF) + on/off smoke

**Files:** none (verification only)

- [ ] **Step 1: Compile the production model through CmdStanR**

```bash
Rscript -e 'library(cmdstanr); set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0"); cmdstan_model("stan/tumor/sf-ssls-lfo.stan", include_paths=c("stan","stan/tumor"), compile=TRUE)'
```
Expected: compiles without error.

- [ ] **Step 2: Confirm default-off is a no-op.** The feature is data-gated
(`enable_background_surrogate=0`), so an existing tumor fit/test must be
byte-unaffected. Run the existing tumor test suite:
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter="tumor")'
```
Expected: pass (no behavioral change when the flag is 0). If no tumor-filtered
tests exist, run the full suite and confirm no new failures vs `main`.

- [ ] **Step 3: Commit (if any fixes were needed); otherwise no-op.**

---

## Task 6: Update spec + memory to reflect Phase 2 done

**Files:**
- Modify: `docs/superpowers/specs/2026-06-08-laplace-surrogate-backgrounded-trials-design.md`
- Modify: `/home/ubuntu/.claude/projects/-mnt-code/memory/laplace-builtin-2-39.md`

- [ ] **Step 1:** In the spec §4, mark integration COMPLETE; record the module
file layout, the `enable_background_surrogate` flag, the per-patient LOD offset,
the summed-variance patient SD, and the tr-sub-hierarchy `fatal_error` guard.
- [ ] **Step 2:** In the memory file, append that Phase 2 integration landed
(modular, flag-gated, old laplace removed), with the production symbol mapping
(`frac_logit_loc_pop`/`init_logit_loc_pop`, summed-variance SD).
- [ ] **Step 3: Commit** (repo files only; memory file saved in place):

```bash
git add docs/superpowers/specs/2026-06-08-laplace-surrogate-backgrounded-trials-design.md
git commit -m "docs(laplace): mark Phase 2 integration complete

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Done criteria

- Old `stan/modules/laplace/` + `r/pioneer/prepare_laplace_data.R` removed; tumor
  model still compiles.
- `stan/modules/laplace_surrogate/` exists as a self-contained module (flags/data/
  functions/transformed_data/likelihood); removing it + its 4 `#include`s fully
  removes the feature.
- `sf-ssls-lfo.stan` compiles; `enable_background_surrogate=0` (default) is a
  verified no-op; `=1` adds the surrogate term, guarded against the tr SD
  sub-hierarchy and a nonzero first anchor.
- R defaults wired (flag off, anchors). Spec + memory updated.
