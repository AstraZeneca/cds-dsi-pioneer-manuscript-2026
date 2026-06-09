# Joint SLD+MS Laplace Surrogate — Full-Model Index Contract (2026-06-09)

Status: **DRAFT for adversarial review.** This is a *contract* document, not an
implementation plan. It specifies the index spaces, variable inventory,
invariants, guards, and a numeric ground-truth example for wiring the **validated
joint SLD+multistate Laplace surrogate** into the **full publication model**
`stan/tumor/sf-ssm-log-space.stan`. Per the 2026-06-09 findings doc, an
adversarial index-contract review must pass *before* any full-model Stan code.

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

### Why no shared-matrix change (the BLOCKER B2 resolution)
`log_cond_surv_01/03` are sized `[n_forecast_patients, max_all_t]` and assembled
in forecast-local row space; they are `#include`d by 8 models (pioneer indexes
them in forecast-local/compact space). **The contract edits only**
`stan/modules/laplace_surrogate/*` **and the model block of**
`sf-ssm-log-space.stan` **— zero changes to shared multistate array shapes.**
The static baseline is *gathered* from the already-computed, group-keyed
`log_level_lambda_01_residual` / `log_level_lambda_03_residual` (NOT
patient-keyed), using the **already-background-safe** gather table
`patient_ms_baseline_flat_idx_slot` (`transformed_data.stan:384-398`, which loops
`1:n_patients` and guards background patients at the patient level → index 1).

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
5. **GQ is untouched** — background patients already produce no generated
   quantities (GQ loops are `1:n_forecast_patients`). No KM-denominator change:
   the two trials have separate KMs; the backgrounded trial simply has no
   simulated curve, as intended.

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
    Same quantity, two computations — **no forecast/background asymmetry**.
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
- **(O2) Burden–frailty prior independence** (§3 zero cross-block). State
  explicitly and check it does not materially bias `coef_01`/`sigma_01` in the
  gate.
- **(O3) 0→3 continuous vs visit-gated.** Confirmed continuous in publication
  (`enable_ms_03_time_varying_cov = TRUE`, not visit-gated). Contract assumes
  continuous; assert if a future config visit-gates 0→3.

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

**The "last value always zero" intent, generalized:** the patient level (last
hierarchy column) of the MSM *baseline* may carry ONLY a marginalized Gaussian
intercept (frailty) on coupled slots — never a GP, FE, or any structure the
surrogate does not integrate out. Guards 1–4 enforce this; if we later want a
richer patient-level baseline under the surrogate, these asserts force a
deliberate revisit rather than a silent wrong answer.

---

## 8. Numeric ground-truth example (DELIVERABLE OF the review workflow)

NOTE: this section is an intentional handoff, not an unfinished spec. The
adversarial review constructs the worked example and verifies hand-computed
indices against the contract. Suggested setup:
- `n_patients = 5`, `n_forecast_patients = 2`, `n_background_patients = 3`;
  `forecast_split_level = patient_lv`, one trial group per cohort.
- `frailty_slots = {0→1, 0→3}` ⇒ `d = 4`; latent vector length `3·4 = 12`.
- Verify: (a) `background_patient_idx` picks the right 3 patients; (b) each bg
  patient's `base01_static` row gathers the correct **trial-group** residual (not
  a patient slot); (c) the 0→1 term sums only over that patient's observed visit
  weeks while 0→3 sums over all weeks; (d) `K` is `12×12` block-diagonal with
  three identical `4×4` blocks `diag(Σ_β, Σ_u)`; (e) the frailty sub-order in θ
  matches `ms_corr_block_member_slot`.

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
2. **Standalone gate** (`laplace_joint_frailty_test.stan` + R driver): joint
   Laplace vs full HMC on synthetic data with frailty on {0→1,0→3} and the
   (level, velocity) coupling; PASS = `coef_01`, `coef_03`, `sigma_01`,
   `sigma_03`, and the 01–03 correlation all agree HMC-vs-Laplace within MCSE,
   clean diagnostics.
3. **This contract** adversarially reviewed and signed off.

Only then: implement the split in `sf-ssm-log-space.stan` (a separate spec +
plan).
