# Joint SLD+MS Laplace Surrogate — Full-Model Index Contract (2026-06-09)

Status: **REVIEWED (go-with-changes) — awaiting human sign-off.** This is a
*contract* document, not an implementation plan. It specifies the index spaces,
variable inventory, invariants, guards, and a numeric ground-truth example for
wiring the **validated joint SLD+multistate Laplace surrogate** into the **full
publication model** `stan/tumor/sf-ssm-log-space.stan`.

Adversarial review: workflow `wf_9c08ae27` (4 dimension reviewers + 2× skeptic
verification + synthesis, all Opus). Verdict **go-with-changes**: one
skeptic-confirmed BLOCKER (the existing `multistate(...)` call reads out-of-bounds
on forecast-sized matrices once a real split exists — see §1 "REQUIRED model-block
edit") plus several major/minor framing corrections, ALL APPLIED 2026-06-09. The
review also REFUTED one reviewer false-positive (a GQ-OOB claim citing dead code;
§5 item 5 corrected). §8 is now a fully-worked, code-verified example.

Prerequisites already satisfied (do not re-litigate):
- Joint SLD+MS(0→1) surrogate validated standalone (`laplace_joint_test.stan`,
  gate `wf_10cb1329` PASS; `coef_01` Laplace 1.212 vs HMC 1.221).
- SLD-only surrogate integrated as modular `stan/modules/laplace_surrogate/`,
  flag-gated behind `enable_background_surrogate` (default 0 = no-op).
- Forecast/background routing already exists and is shared:
  `_full_model_transformed_data.stan:100-120` builds `forecast_patient_idx`,
  `background_patient_idx`, `n_background_patients`.
- cmdstan 2.39 toolchain pinned (env + `.Rprofile` + hash + guard).

---

## 0. Design decisions locked in this session

1. **Background MS transition set is config-driven**, resolving to whichever
   transitions have burden coupling and/or patient-level frailty in the active
   config. For the *publication* config that is **{0→1, 0→3}** (0→2 coupling is
   `FALSE`).
2. **Background 0→1 hazard honors observed visit weeks** (visit-gated latent
   mode), NOT a dense weekly grid. 0→3 is continuous (every week).
3. **Marginalize the patient-level frailty too** (do not forbid it). The
   publication enables a correlated patient-level RE-NCP intercept on slots
   {0→1, 0→3} (`publication_targets.R:401-406`); background patients' frailty is
   integrated out alongside the burden latents.
4. **Deliverable this session = this contract only.** No implementation plan;
   hand to adversarial review first.

---

## 1. The central object: hazard decomposition for background patients

For a background patient, each coupled transition's log-hazard at week `w`
splits into a part that is **constant in the marginalized latents** and a part
that is **linear in them**:

```
loghaz_T(w) = [ baseline_T(group(p), w)            ]   ← θ-INDEPENDENT (static)
            + [ time_invariant_coef_T · covar(p)   ]   ← θ-INDEPENDENT (static)
            + [ u_T                                ]   ← θ-dependent: frailty latent (if slot T has frailty)
            + [ tv_coef_T[1]·std_level(g(w)) + tv_coef_T[2]·std_vel(g'(w)) ]  ← θ-dependent: burden coupling
```
where `T ∈ {0→1, 0→3}`, `group(p)` is the patient's trial group,
`g(w) = b0 + b1·w + b2·w²` is the surrogate quadratic log-burden, and
`g'(w) = b1 + 2·b2·w` its velocity (the resolved 2-feature basis — §6 O1). Both
features are linear in `(b1,b2)`.

The static parts are computed once in the model block (parameter-dependent but
`θ`-independent) and **passed into the functor as ordinary arguments** — they are
NOT recomputed inside the Laplace integral, and NOT gathered from any
forecast-sized matrix (background patients have no row there). Only the
`u_T` and `tv_coef` terms vary with `θ`, and both are **linear in log-hazard
space**, preserving log-concavity (the validated property).

### Why no shared-matrix-SHAPE change (the BLOCKER B2 resolution)
`log_cond_surv_01/03` are sized `[n_forecast_patients, max_all_t]` and assembled
in forecast-local row space; they are `#include`d by 8 models (pioneer indexes
them in forecast-local/compact space). **The SURROGATE WIRING touches only**
`stan/modules/laplace_surrogate/*` **and the model block of**
`sf-ssm-log-space.stan` **— zero changes to shared multistate array *shapes*.**
The static baseline is *gathered* from the already-computed, group-keyed
`log_level_lambda_01_residual` / `log_level_lambda_03_residual` (NOT
patient-keyed). The gather uses `patient_ms_baseline_flat_idx_slot`
(`transformed_data.stan:384-398`).

**Gather-safety, stated precisely (review-corrected).** Correctness does NOT rely
on a general "is-background" guard. It relies on gathering the **trial-level**
(`lv < n_levels`) residual, where `get_global_group_idx` returns the correct trial
group for *every* patient including background. The `→ 1` fallback at
`transformed_data.stan:391` fires only at the **patient level** (`lv == n_levels`,
keyed on `patient_level_groups[i,lv] > n_forecast_patients`) — a different
predicate than the routing key — and that patient-level column is the marginalized
frailty, which is never gathered (it's in `θ`). The gather therefore **silently
depends on a trial-level (non-patient) split**; see Guard 6 (§7).

**SEPARATELY**, the upstream O1 (level, velocity) coupling-basis change DOES edit a
shared file (`stan/_ms_burden_tv_covar.stan`, `#include`d by `pioneer.stan:99`
+ 4 tumor models): it alters `n_time_varying_covar` (3→2) and adds
velocity-standardization constants. That is a **gated prerequisite (§9 step 0),
NOT part of the surrogate wiring** — the blanket "zero shared changes" claim
applies to the surrogate wiring only.

### REQUIRED model-block edit: restrict the existing likelihood to forecast patients (BLOCKER, review-confirmed)
The surrogate wiring is NOT purely additive. `sf-ssm-log-space.stan:101-119` calls
`ms_final_state ~ multistate(ones_vector(n_patients), … ms_time_01, … log_cond_surv_01,
… log_cond_surv_03, …)` passing **full unified arrays** (length `n_patients`) and
never slicing by `ms_patient_idx`. Inside `multistate_lpmf` the loop is
`for (i in 1:n_patients)` (`multistate.stanfunctions:963,966`), reading
`log_cond_surv_01[i]` / `log_cond_surv_03[i]` — matrices sized
`[n_forecast_patients, max_all_t]` filled only for `j in 1:n_forecast_patients`.
Under the publication split (`forecast_split_level=1L`, background patients exist)
rows `i > n_forecast_patients` **read OUT OF BOUNDS → runtime crash before the
surrogate contributes.** The dense tumor SLD loop (`:91-95`) has the dual problem:
it slices forecast-local `states` with unified `patient_visit_pos`.

When `enable_background_surrogate == 1`, the model-block edit MUST:
1. **Slice the multistate call by `ms_patient_idx`** exactly as
   `modules/multistate/likelihood.stan:5-13` already does
   (`ms_final_state[ms_patient_idx]`, `ms_time_01[ms_patient_idx]`,
   `ones_vector(n_forecast_patients)` etc.) so it runs forecast-only and never
   indexes a background row.
2. **Gate the dense SLD obs loop over `forecast_patient_idx`** so background SLD is
   handled SOLELY by the surrogate (the surrogate **replaces, not supplements**,
   the dense likelihood for background patients — avoid double-counting).
3. **Fix the hardcoded `enable_ms_visit_gated_01 = 0` argument** at
   `sf-ssm-log-space.stan:118` to the actual flag value (publication: `TRUE`,
   `publication_targets.R:312`) — a pre-existing inconsistency that would otherwise
   assemble the 0→1 hazard differently in the dense path vs the surrogate.

---

## 2. Parametric latent layout (NOT hardcoded to 4)

Per-background-patient latent vector:

```
d = 2 + n_frailty_slots

θ_i = ( b1_i, b2_i,         ← burden RE (always present; drive SLD + burden coupling)
        u_{s1,i}, …, u_{sF,i} )   ← one frailty latent per slot in `frailty_slots`
```

- `n_frailty_slots = |frailty_slots|`, where `frailty_slots` = coupled slots with
  a patient-level RE-NCP baseline intercept (`ms_level_intercept_mode[slot,
  patient_lv] == 2`).
- **Publication**: `frailty_slots = {0→1, 0→3}` ⇒ `d = 4`. **SLD-only fallback**:
  `frailty_slots = ∅` ⇒ `d = 2`, collapsing to the already-validated surrogate.
- Two INDEPENDENT config-driven sets — the contract must not assume they coincide:
  - `coupled_tv_slots`: slots whose hazard reads the burden coupling
    (`tv_coef_T`). Drives **functor body only**.
  - `frailty_slots`: slots with a marginalized patient-level intercept. Drives
    **latent dim + the `K` covariance block**.
  - They happen to both be `{0→1, 0→3}` for the publication, but are specified
    and asserted separately.

### Stacking convention (CONTRACT-FIXED)
Within patient `i`'s block of length `d`: burden latents first
`[b1, b2]`, then frailty latents in **ascending MS slot order**
(`MS_SLOT_01 < MS_SLOT_02 < MS_SLOT_03 < …`). So for the publication:
`[b1, b2, u_01, u_03]`. The frailty sub-order MUST match the `Σ_u` block's member
order (§3) and the `ms_corr_block_member_slot` ordering — a single source of
truth; assert equality.

---

## 3. Prior covariance `K_i` (block-diagonal, dim `d×d`)

```
K_i =  ⎡ Σ_β    0  ⎤      Σ_β : 2×2 burden-RE covariance (EXISTING: surrogate_bridge output)
       ⎣  0    Σ_u ⎦      Σ_u : (n_frailty_slots)² correlated-frailty covariance
```

- **`Σ_β`** — unchanged from the SLD-only surrogate: `surrogate_bridge(...)`
  returns it as the GH-3 pushforward marginalized-coef covariance.
- **`Σ_u`** — the SAME construction the module already uses
  (`transformed_parameters.stan:37`):
  `Σ_u = diag(σ) · (L Lᵀ) · diag(σ)`, where
  - `σ[m] = log_lambda_gp_<slot_m>_level_intercept_sd[patient_lv]` for each
    `slot_m ∈ frailty_slots` (reuses the existing per-transition patient-level
    intercept SD — no new scale hyperparameters).
  - `L = L_ms_intercept_corr[b]` for the patient-level correlation block `b`
    (`ms_corr_block_level[b] == patient_lv`).
- **Zero cross-covariance** between burden RE and frailty RE: they are independent
  a priori (burden is a tumor-trajectory RE; frailty is a hazard-intercept RE).
  This is a modeling assumption to be stated in the spec and checked in the gate.
- `surrogate_K_fn` is extended to tile this `d×d` block per patient (currently
  tiles a 2×2). `hessian_block_size` becomes `d` (currently 2).
- **Single-source-of-truth for `d` (review-added implementation note):** the
  literal `int d = 2` is hardcoded in THREE places today — `surrogate_ll`
  (`surrogate.stanfunctions:88`), `surrogate_K_fn` (`:109`), and the
  `transformed_data` `theta_0` / `hessian_block_size` setup. ALL must be driven
  from one `d = 2 + n_frailty_slots`, with an assert that
  `surrogate_hessian_block_size == d`.

---

## 4. Extended `surrogate_ll` functor signature

Current (SLD-only): `surrogate_ll(theta, beta_pop, measure_sd,
log_lod_per_visit, n_bg, bg_obs, bg_pos, bg_time)`.

Extended (joint). New arguments are all `data`-qualified static inputs OR
parameter-derived statics passed positionally (Stan's `laplace_marginal_tol`
allows non-data args in the parameter tuple):

| Arg | Qual | Meaning | Index space |
|---|---|---|---|
| `theta` | — | stacked latents, length `n_bg·d` | per-patient blocks of `d` (§2) |
| `beta_pop` | — | burden bridge mean `[b0,b1,b2]` | length 3 (existing) |
| `measure_sd` | — | SLD obs SD | scalar |
| `bg_obs/bg_pos/bg_time/log_lod_per_visit` | data | SLD term inputs | **bg-compact** (existing) |
| `base01_static[bg, week]` | param-static | 0→1 baseline + TI-cov offset | bg-compact rows; **observed visit weeks** |
| `base03_static[bg, week]` | param-static | 0→3 baseline + TI-cov offset | bg-compact rows; **all weeks 1..te** |
| `bg_visit_wk_01` (+pos) | data | observed visit weeks per bg patient (0→1 gating) | bg-compact ragged |
| `bg_event_wk_01/03`, `bg_cens_01/03` | data | event week + censor flag per transition | bg-compact |
| `tv_coef_01[k], tv_coef_03[k]` | param | burden-coupling coefficients | length 2 per slot: [level, velocity] (§6 O1) |
| `median_vel, iqr_vel` | data | velocity standardization (NEW, §6 O1) | scalar |
| `median_lb, iqr_lb` | data | burden standardization | scalar |
| `frailty_offset[bg, f]` | (in θ) | — | frailty latents live in θ, NOT a separate arg |

**Functor body** per background patient `i` (log-concave; mirrors
`laplace_joint_test.stan` but vector baseline + frailty):
1. SLD term over `bg_pos[i]` visits (existing, with LOD `normal_lcdf`).
2. 0→1 term: for each observed visit week `w` in `bg_visit_wk_01[i]` up to
   event/censor: `loghaz = base01_static[i,w] + u_01_i + tv_coef_01[1]·std_level(g(w))
   + tv_coef_01[2]·std_vel(g'(w))`, where `g(w)=b0+b1·w+b2·w²`,
   `g'(w)=b1+2·b2·w` (both linear in the latents); accumulate `−exp(loghaz)`; add
   `loghaz` at the event week if uncensored.
3. 0→3 term: same, over **all** weeks `1..te_03(i)`, with `u_03_i`, `tv_coef_03`.

**Week-index convention (PIN before the standalone gate — review-flagged
off-by-one risk).** Let `w` be the single week index used consistently for BOTH
the SLD-μ term and the hazard burden `g(w)`. Production currently carries THREE
conventions: surrogate SLD `t = bg_time[v] − 1.0` (anchored idx, baseline→0,
`surrogate.stanfunctions:97`); visit-gated 0→1 `wk = t_patient_visits[v]` (RAW
calendar week used to index `log_cond_surv_01[j,wk]`, `transformed_parameters.stan:197`);
toy `w ∈ 1:te` with no offset (`laplace_joint_test.stan:50,59`). The joint functor
MUST feed `g(w)` with the SAME `w` at which `base01_static[i,w]` is gathered (the
residual column = calendar week), and the `−1` SLD offset MUST be applied
identically to the SLD and hazard burden origins, or dropped from both. Resolve
this explicitly before the gate; an unpinned offset makes `g(w)` inconsistent
across the SLD and hazard terms.

---

## 5. Index-space facts the contract pins down

1. **bg-compact space** (index `1..n_background_patients`) via
   `background_patient_idx[i] = p` (unified patient id). All functor inputs are
   in bg-compact rows; the existing `surrogate_bg_pos/obs/time`
   (`laplace_surrogate/transformed_data.stan:48-76`) establish the pattern.
2. **Static baseline gather**: for `p = background_patient_idx[i]`,
   `base01_static[i, ·] = log_level_lambda_01_residual[ patient_ms_baseline_flat_idx_slot[MS_SLOT_01, p, lv], · ]`
   summed over enabled baseline levels `lv`, plus `time_invariant_coef_01 ·
   covar(p)`. The gather table is already full-patient-sized and
   background-safe. **0→1 restricted to observed visit weeks; 0→3 over all weeks.**
3. **Visit weeks** come from `t_patient_visits[patient_visit_pos[p] ..]` (unified),
   compacted to `bg_visit_wk_01`. Mirror the existing `surrogate_bg_time` build.
4. **Event/censor weeks** from `ms_time_01/ms_time_03` + `ms_censored_*` (unified),
   compacted to bg space.
5. **GQ is background-safe** (review-scoped — do NOT assert a blanket "all GQ loops
   are `1:n_forecast_patients`"). The GQ actually compiled into
   `sf-ssm-log-space.stan` is `modules/state_space/generated_quantities.stan`
   (included at line 201), which loops `1:n_forecast_patients`. The one GQ loop
   over `1:n_patients` (`tumor/generated_quantities.stan:22`) indexes ONLY
   unified-sized arrays (`recist`/`rep_recist`, sized `n_total_visits`), so it is
   in-bounds. (`_tumor_endpoints_generated_quantities.stan` is dead code — not
   `#include`d — and irrelevant.) **Conclusion holds:** background patients produce
   no per-patient forecast curves; the two trials have separate KMs and the
   backgrounded trial simply has no simulated curve, as intended. No
   KM-denominator change.

---

## 6. OPEN modeling decisions (must resolve before the standalone gate)

These are NOT index issues; they are modeling choices the contract surfaces:

- **(O1) RESOLVED by a model change — burden coupling basis = (level, velocity).**
  Decision (2026-06-09): replace the current 3-feature coupling (standardized
  log-burden + log-decrease-rate + log-growth-rate, `_ms_burden_tv_covar.stan:64-80`)
  with a **2-feature trajectory-derivative basis, MODEL-WIDE** (forecast AND
  background):
  - **Feature 1 = level**: standardized log-burden `g(w)` (unchanged).
  - **Feature 2 = velocity**: `g'(w)`, the rate of change of log-burden.
    Forecast patients: numerical derivative across `states_full_grid` columns.
    Background patients: analytic derivative of the quadratic, `b1 + 2·b2·w`.
    Analytically consistent in **FORM** (level + velocity), but NOT bit-identical
    (review-corrected): the background path is the uncapped analytic quadratic
    `g(w)=b0+b1·w+b2·w²` / `g'(w)=b1+2·b2·w`, a 3-anchor quadratic **surrogate** of
    the forecast trajectory, whereas the forecast feature is the capped
    (`fmin(·,10)`), bi-exponential, process-noise-perturbed `log_sum_exp` burden
    (`_ms_burden_tv_covar.stan:59-62`) with a numerical-derivative velocity. The
    surrogate is a deliberate log-concave approximation. Because `tv_coef` is
    shared model-wide and calibrated on the FORECAST feature scale, the standalone
    gate (§9 step 2) MUST add a forecast-vs-surrogate **feature-scale** comparison.
  - Both are **linear in (b1,b2)** ⇒ survival term stays log-concave; the
    surrogate is exact-in-form, not an approximation of a different model.
  - The bi-exponential decrease/growth-rate covariates are **dropped**.
  - **NEW transformed-data constants required (BOTH paths):** velocity needs its
    own centering/scaling (`median_velocity_obs`/`iqr_velocity_obs` or analogous),
    computed from observed-trajectory numerical derivatives — log-burden's
    `median/iqr` do not apply to a velocity. Without this, `time_varying_coef`
    for velocity is on an uninterpretable scale and its `Normal(0,σ)` prior is
    mis-calibrated.
  - **This is a CORE MS model change**, not a surrogate-only detail. It alters
    `coef`'s clinical interpretation (now a tumor-kinetics "where + which way"
    coupling) and affects forecast patients. It needs its OWN spec + clinical
    re-justification + a refit comparison vs the current 3-feature model, SEPARATE
    from this index contract. This contract assumes the (level, velocity) basis as
    its input.
- **(O2) Burden–frailty prior independence is PRIOR-EXACT** (§3 zero cross-block) —
  not an approximation (review-corrected). Burden REs use independent `std_normal`
  NCP priors (`tr/priors.stan:36`); the frailty `ms_corr_u` is a separate
  LKJ+Gaussian block with no cross-covariance term. So the block-diagonal `K`
  introduces **zero prior approximation error**. The gate check is confirmation,
  not validation of an assumption.
- **(O3) 0→3 continuous vs visit-gated. CORRECTED 2026-06-11 (impl review
  `wf_e88e7d91`): production integrates 0→3 VISIT-GATED, not continuous.** The
  earlier "confirmed continuous" reading was WRONG — every `log_cond_surv_03`
  integration site in `multistate.stanfunctions` (state-0 censored :458, state-1
  :488, state-2 :659/:684, ...) uses `sum_at_visits_below`, i.e. assessment-visit
  weeks only, for ALL patient patterns regardless of `enable_ms_visit_gated_01`.
  `enable_ms_03_time_varying_cov = TRUE` controls whether the burden COVARIATE is
  read, NOT the integration grid. The production surrogate functor was corrected to
  visit-gate 0→3 (commit `53139a75` on branch karim/laplace). **RESOLVED 2026-06-11:
  the d=4 gate was RE-RUN with a visit-gated-0→3 DGP (sparse 12/36-week visit grid,
  67%/64% event rates) and PASSED** — HMC-vs-Laplace max |diff| = 0.10 (frailty SDs,
  a shared identifiability limit), all coupling/baseline diffs ≤ 0.04, sd_ratio
  0.96–1.09 (Laplace NOT overconfident — the key correctness signal), HMC ref Rhat
  ≤ 1.03. Materially cleaner than the original continuous-0→3 gate's borderline 5.59
  diff/MCSE. The visit-gated marginalization is validated for production.

---

## 7. Guards / asserts (in `laplace_surrogate/transformed_data.stan`)

Triggered only when `enable_background_surrogate == 1`:

1. **No patient-level GP on any coupled slot**:
   `enable_ms_level_gp[slot, patient_lv] == 0` for `slot ∈ coupled_tv_slots ∪
   frailty_slots`. A time-varying patient-level residual is not a marginalizable
   scalar → `fatal_error`. (Publication: these are 0, lines 403-404 ✓.)
2. **Patient-level baseline intercept must be marginalizable**: for every coupled
   slot, `ms_level_intercept_mode[slot, patient_lv] ∈ {0 (off), 2 (RE-NCP)}`. Any
   other mode (e.g. FE) on a coupled slot → `fatal_error` with message
   "patient-level baseline beyond marginalized Gaussian frailty is unsupported by
   the surrogate; revisit when enabling it."
3. **Frailty-set consistency**: `n_frailty_slots` (the latent-dim driver) MUST
   equal the count of coupled slots with `ms_level_intercept_mode[·,patient_lv]
   == 2`. The latent vector and `K` can never silently drift from the model's
   actual frailty structure → `fatal_error` on mismatch.
4. **Correlation block alignment**: the patient-level correlation block's member
   slots (`ms_corr_block_member_slot[b,·]` where `ms_corr_block_level[b] ==
   patient_lv`) MUST equal `frailty_slots` in the same order as the §2 stacking →
   `fatal_error` on mismatch (else `Σ_u` rows misalign with the latents).
5. **Existing guards retained**: tr SD sub-hierarchy incompatibility
   (`transformed_data.stan:28-33`), first anchor == 0 (lines 36-38).
6. **Split must be at a non-patient (trial) level** (review-added): assert
   `forecast_split_level > 0 && forecast_split_level < n_levels` when
   `enable_background_surrogate == 1`. The baseline gather silently depends on
   this — a patient-level split would (a) select at most ONE forecast patient
   (`forecast_group` is a single scalar, `_full_model_data.stan:23`; patient level
   is identity-enforced, `_hierarchy_transformed_data.stan:29`) and (b) push the
   gather into the `→1` patient-level fallback. → `fatal_error` otherwise.

**The "last value always zero" intent, generalized:** the patient level (last
hierarchy column) of the MSM *baseline* may carry ONLY a marginalized Gaussian
intercept (frailty) on coupled slots — never a GP, FE, or any structure the
surrogate does not integrate out. Guards 1–4 enforce this; if we later want a
richer patient-level baseline under the surrogate, these asserts force a
deliberate revisit rather than a silent wrong answer.

---

## 8. Numeric ground-truth example (worked, verified against code)

**Setup.** `n_patients = 5`, split at the **TRIAL level**:
`forecast_split_level = 1` (trial level; matches `publication_targets.R:258`
`forecast_split_level = 1L`), `forecast_group = 1` (= trial A's integer code).
Patients in id order: `p1,p2 ∈ trial A` (code 1), `p3,p4,p5 ∈ trial B` (code 2).
(NB a *patient*-level split is impossible here: `forecast_group` is a single
scalar — `_full_model_data.stan:23` — and the patient level is identity-enforced
— `_hierarchy_transformed_data.stan:29` — so it would pick exactly one forecast
patient. The split MUST be at the trial level for a 2/3 cohort split.)

`frailty_slots = {0→1, 0→3}` ⇒ `d = 2 + 2 = 4`; latent vector length `3·4 = 12`.

### (a) Routing (`_full_model_transformed_data.stan:104-119`)
Loop `p = 1..5` on key `patient_level_groups[p,1] != forecast_group`:
- p1 → A == 1 → forecast; `forecast_patient_idx[1] = 1`
- p2 → A == 1 → forecast; `forecast_patient_idx[2] = 2`
- p3 → B != 1 → background; `background_patient_idx[1] = 3`
- p4 → B → background; `background_patient_idx[2] = 4`
- p5 → B → background; `background_patient_idx[3] = 5`

⇒ `forecast_patient_idx = [1,2]` (`n_forecast_patients = 2`),
`background_patient_idx = [3,4,5]` (`n_background_patients = 5−2 = 3`).
Guard `j_fc == n_forecast_patients` (line 117) passes. bg-compact `i ∈ {1,2,3}`
maps to unified patient `{3,4,5}`.

### (b) Static baselines for bg-compact `i = 1` (⇒ `p = 3`, trial B)
Publication enables the MS baseline hazard at the **trial level only**
(`enable_ms_level_baseline_hazard` trial=3, patient=0; `publication_targets.R:292`),
so only the trial level contributes a residual row; the patient-level baseline is
mode 0 (no row). With `lv = trial` (`< n_levels`) the gather guard at
`transformed_data.stan:391` never fires, so:

- **0→1** (visit-gated; `enable_ms_visit_gated_01 = TRUE`,
  `publication_targets.R:312`):
  `base01_static[1, w] = log_pop_lambda_01[w]
   + log_level_lambda_01_residual[ g01, w ] + time_invariant_coef_01·covar(p3)`,
  where `g01 = patient_ms_baseline_flat_idx_slot[MS_SLOT_01, 3, trial_lv]`
  = `get_global_group_idx(...)` = trial B's **group row** (NOT a patient slot),
  the same group-keyed gather the forecast path uses at
  `transformed_parameters.stan:164`. The 0→1 hazard sum runs **only over p3's
  observed visit weeks** `w ∈ bg_visit_wk_01[1]`.

- **0→3** (continuous; `enable_ms_03_time_varying_cov = TRUE`,
  `publication_targets.R:329`; dense full-matrix add at
  `transformed_parameters.stan:780`, no visit-gating branch):
  `base03_static[1, w] = log_pop_lambda_03[w]
   + log_level_lambda_03_residual[ g03, w ] + time_invariant_coef_03·covar(p3)`,
  `g03 = patient_ms_baseline_flat_idx_slot[MS_SLOT_03, 3, trial_lv]`
  (same group-keyed gather, `transformed_parameters.stan:767`;
  `enable_ms_03_time_invariant_cov = TRUE` so the TI-cov-03 offset is folded in).
  The 0→3 hazard sum runs over **all** weeks `1..te_03(p3)`.

**Week-index convention (pin before gate).** `w` is the single week index used
identically for the SLD-μ term and the hazard burden `g(w)` (see §4 pin).

### (c) bg-compact ragged arrays (`laplace_surrogate/transformed_data.stan:48-76`)
Let p3 visit weeks `{0,1,2,3}` (4 visits), p4 `{0,1,2,3,4}` (5), p5 `{0,1,2}` (3),
baseline week 0 anchored to `t_patient_visit_idx = 1`. Pos build (lines 51-57):
`surrogate_bg_pos[1]=1`, `[2]=1+4=5`, `[3]=5+5=10`, `[4]=10+3=13` ⇒
`surrogate_bg_pos = [1,5,10,13]`, `surrogate_n_bg_visits = 12`.
`surrogate_bg_time`: p3→`[1,2,3,4]`, p4→`[1,2,3,4,5]`, p5→`[1,2,3]`.
The new `bg_visit_wk_01` (+ parallel pos array) carries the SAME per-patient
post-baseline candidate weeks the 0→1 hazard is gated on; the contract must state
whether the baseline/screening visit (anchored idx 1) is excluded from the 0→1
hazard sum.

### (d) Prior covariance `K` (`12×12`, block-diagonal)
`n_bg = 3`, `d = 4` ⇒ `K` is `12×12` = `blockdiag(B, B, B)` with three IDENTICAL
`4×4` blocks (Σ_β, Σ_u are population-level, tiled per patient exactly as
`surrogate_K_fn` tiles one Σ_β today, `surrogate.stanfunctions:110-113`):

```
B = ⎡ Σ_β   0  ⎤     Σ_β = ⎡ Var(b1)      Cov(b1,b2) ⎤   (surrogate_bridge GH-3 pushforward)
    ⎣  0   Σ_u ⎦           ⎣ Cov(b1,b2)   Var(b2)    ⎦

                          Σ_u = ⎡ s01²          ρ·s01·s03 ⎤   (transformed_parameters.stan:37)
                                ⎣ ρ·s01·s03     s03²      ⎦
```
`s01 = log_lambda_gp_01_level_intercept_sd[patient_lv]`,
`s03 = log_lambda_gp_03_level_intercept_sd[patient_lv]`
(`transformed_parameters.stan:29,31`); `ρ = (L Lᵀ)[1,2]`,
`L = L_ms_intercept_corr[b]` for the patient-level block. Off-diagonal blocks are
all 0: cross-patient by prior independence; burden↔frailty by §3 zero cross-block
(PRIOR-EXACT — see O2).

### (e) θ stacking and frailty sub-order
```
θ = [ b1_p3, b2_p3, u01_p3, u03_p3,   ← i=1
      b1_p4, b2_p4, u01_p4, u03_p4,   ← i=2
      b1_p5, b2_p5, u01_p5, u03_p5 ]  ← i=3
```
Burden-first matches the validated functor
(`laplace_joint_test.stan:43-44`, `surrogate.stanfunctions:92-93`). Frailty
sub-order `[u01, u03]`: `ms_corr_block_member_slot[b,·]` is filled by looping
`k in 1:N_MS_INTERCEPT_SLOTS` in ASCENDING slot order
(`transformed_data.stan:207-210`), with `MS_SLOT_01=1 < MS_SLOT_03=3` both in
patient-level corr group 1 (`publication_targets.R:405-406`) ⇒ member m=1→slot 01,
m=2→slot 03. `Σ_u` row/col order therefore equals the θ stacking `[u01,u03]`,
satisfying the §7 guard-4 assert.

---

## 9. Prerequisites before ANY full-model code (discipline + ordering)

Mirroring what has worked on this branch. **Ordered dependency chain:**

0. **(UPSTREAM, separate spec) Core model change: (level, velocity) coupling
   basis.** O1 is resolved in principle but the model change itself must land and
   be validated FIRST — it is not surrogate-specific. Needs: its own spec,
   clinical re-justification of `coef` as a kinetics coupling, the new velocity
   standardization constants on the forecast path, and a refit comparison vs the
   current 3-feature model. **This contract's functor assumes the (level,
   velocity) basis as input**, so it cannot be validated until step 0 exists.
1. **R premise-check**: confirm the `d`-latent joint log-density is log-concave
   with the frailty intercepts included (extend the existing Hessian-eig check;
   the validated joint check was for `d=2` burden-only — frailty adds linear
   shifts, expected to remain concave, but VERIFY numerically). Use the (level,
   velocity) basis from step 0.
2. **Standalone gate** (`laplace_joint_frailty_test.stan` +
   `test_laplace_joint_frailty.R`): joint Laplace vs full HMC on synthetic data
   with correlated frailty on {0→1,0→3} and the (level, velocity) coupling.
   **STATUS (2026-06-11): conditional PASS — marginalization accepted as sound.**
   Final run (n=45, 82%/64% event rates, d=4 = b1,b2,u01,u03; ~5h):
   - Burden params + baselines: excellent agreement (diff/MCSE 0.1–0.6).
   - Coupling + frailty (cf_lvl/vel_01/03, s01, s03, rho): agree directionally;
     max diff/MCSE **5.59** (just over the 5 threshold).
   - **min SD ratio 0.96** (Laplace posteriors appropriately wide, NOT
     overconfident — the key correctness signal); 0 divergences both modes.
   - The borderline 5.59 is a **noisy-reference artifact**, NOT a Laplace error:
     the HMC reference under-mixed (`min_ess_bulk=142`, `rhat 1.04`), inflating the
     MCSE denominator; raw param diffs are tiny and Laplace tracks HMC in the same
     direction on every parameter. The 01↔03 correlation `rho` remains
     under-identified even at these event rates (HMC −0.19 vs true −0.5), a data
     property both methods share, not a marginalization defect.
   - **Decision (Karim, 2026-06-11):** accept as effectively-PASS on the weight of
     evidence rather than burn another ~5h chasing a cleaner HMC reference. The
     correlated-frailty marginalization is validated for production wiring.
   - NOTE: the gate uses the toy (level, velocity) hazard directly; the
     forecast-vs-surrogate feature-SCALE check (per O1) still belongs in the
     production wiring validation, where the real capped/bi-exponential forecast
     feature meets the uncapped quadratic surrogate feature.
3. **This contract** adversarially reviewed (workflow `wf_9c08ae27`, verdict
   go-with-changes; all skeptic-confirmed required changes applied 2026-06-09) —
   awaiting human sign-off.

Only then: implement the split in `sf-ssm-log-space.stan` (a separate spec +
plan). **The split edit MUST include the model-block likelihood restriction
specified in §1 (slice multistate by `ms_patient_idx`, gate dense SLD over
`forecast_patient_idx`, fix the hardcoded `enable_ms_visit_gated_01`).**
