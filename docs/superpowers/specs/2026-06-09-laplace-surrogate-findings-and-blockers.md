# Laplace Surrogate — Consolidated Findings & Blockers (2026-06-09)

Status snapshot of the `karim/laplace` work. The tumor (SLD-only) surrogate is
validated; the **JOINT SLD+multistate surrogate is now also validated** (the
chosen design); enabling either on real data requires correctly implementing the
forecast/background split, which surfaced pre-existing model defects.

## DESIGN DECISION (2026-06-09): JOINT SLD+MS surrogate, not SLD-only

The publication model uses the **full model** `stan/tumor/sf-ssm-log-space.stan`
(NOT the LFO model). Its MS hazards for 0→1 (progression, visit-gated) and 0→3
(dropout) are DRIVEN by each patient's tumor trajectory (a time-varying
covariate). Marginalizing background patients' tumor latents from the SLD term
ALONE (the original surrogate) would sever that burden→hazard coupling for
background patients — and `coef_01` (the coupling coefficient) is the model's
headline scientific quantity, with ~67% of patients in the background trial.

**Resolution (validated):** marginalize background patients' (b1,b2) from the
JOINT likelihood — SLD + the burden-coupled 0→1/0→3 piecewise-exponential hazard
— feeding the SURROGATE QUADRATIC burden (linear in (b1,b2)) into the hazard.
This keeps the survival term log-concave (event term linear; −∫exp(linear)
concave), so Laplace stays valid, AND preserves the coupling inside the
marginalization (no info loss, no imputation).

Validation evidence:
- **R log-concavity check**: joint log-density max Hessian eig = −6188 (concave
  everywhere); bi-exponential control = +188 (the original NO-GO, reproduced).
- **Option-A loss quantified**: dropping the coupling term for background
  discards ~67% of the patients informing `coef_01`.
- **Standalone joint gate (`laplace_joint_test.stan` + `test_laplace_joint.R`),
  workflow `wf_10cb1329`: PASS.** Joint-Laplace vs full HMC, 30 synthetic
  patients: max diff/MCSE 1.76 (<5), min SD ratio 0.98, all clean (0 div, R-hat
  1.00). `coef_01`: HMC 1.221 vs joint-Laplace 1.212 (diff/MCSE 1.53) — coupling
  recovered identically. (Both differ from true 0.8 by a shared small-data/8-event
  artifact, NOT a Laplace error — the gate compares Laplace-vs-HMC, not vs truth.)
  Caveat: solver 1 fell back to solver 2 (numerical conditioning of the (t,t²)
  basis, not loss of concavity); harmless, fit matched HMC. A centered/scaled time
  basis would likely let solver 1 hold.

## TL;DR dependency chain (revised)

```
[DONE] SLD-only surrogate (validated PFS gate; integrated into LFO, flag-gated)
[DONE] JOINT SLD+MS surrogate design (validated standalone gate — the CHOSEN design)
[DONE] cmdstan 2.39 toolchain management (env + pin + hash + guard)
   │
   └─ to run on publication (FULL model sf-ssm-log-space.stan) needs:
        (P) implement forecast/background split correctly in the FULL model
            - tumor likelihood: forecast patients via states; background via
              JOINT surrogate (SLD+MS), no states rows
            - MS: full-model-LOCAL likelihood (do NOT widen shared matrices —
              breaks pioneer); background patients learn hazards via the
              joint-marginalized z
            - index contract from wf_fca8a194 (Option-A-scoped) must be REVISED
              for the joint design before implementation
        (B1, issue #92) LFO model MS likelihood is 0→1-only [separate bug; the
            full model is unaffected — it uses multistate_lpmf over all transitions]
```

The surrogate math (both SLD-only and joint) works. What remains is implementing
the forecast/background split in the full model around the joint surrogate.

---

## DONE (validated, committed, pushed on `karim/laplace`)

### Tumor surrogate (Phase 1 + 2)
- Quadratic-in-time log-concave surrogate marginalizes backgrounded patients'
  tumor latents via `laplace_marginal_tol`. Design: intercept pinned, marginalize
  `(β₁,β₂)`, GH-3 pushforward bridge for both mean and covariance.
- Validated by mixed-cohort + **PFS-endpoint** gate: median PFS identical, max
  quantile diff 1 wk, event-rate diff 0.019. See
  `2026-06-08-laplace-surrogate-backgrounded-trials-design.md` §5a.
- Integrated as modular `stan/modules/laplace_surrogate/` (flags/data/
  surrogate.stanfunctions/transformed_data/likelihood), gated behind
  `enable_background_surrogate` (default 0 = no-op). Old hand-coded
  `stan/modules/laplace/` removed.

### cmdstan 2.39 toolchain management ("renv for cmdstan")
Root cause of jobs #1885-1893 chain crashes: model compiled against stanc 2.38
(no `laplace_marginal_tol`). Three layers now protect it:
- **env**: 2.39 is the latest/default cmdstan in the Domino image (the compile
  subprocess reads the default).
- **`.Rprofile` `init_project()`**: pins `set_cmdstan_path(2.39)`, env-overridable
  via `PIONEER_CMDSTAN_PATH`.
- **`build_model` guard**: fails loud if a `laplace_marginal`-using model meets
  cmdstan < 2.39.
- **`compute_stan_source_hash`**: folds cmdstan version into the binary-cache hash
  so a version switch auto-invalidates the compiled exe.
- Also: `renv.lock` snapshotted (renv 1.1.5→1.2.3, R 4.3.2→4.5.0, +httr2).
- CLAUDE.md note: never use `!!!` (rlang splice) inside a `tar_target()` command
  (evaluates at manifest-parse time → object-not-found during job setup).

---

## BLOCKER B1: LFO multistate likelihood is 0->1-only (PRE-EXISTING BUG)

**Confirmed defect, independent of the surrogate.** Decision (2026-06-09): this is
a bug — the LFO model should learn from all enabled transitions like the full model.

- **Full model** `stan/tumor/sf-ssm-log-space.stan:101`:
  `ms_final_state ~ multistate(... enable_ms_01, enable_ms_02, enable_ms_12,
  enable_ms_03, enable_ms_32 ...)` — full `multistate_lpmf`, learns from ALL
  enabled transitions.
- **LFO model** `stan/tumor/sf-ssls-lfo.stan:108`: only
  `if (enable_ms_01) calc_ms_single_transition_loglik(...)` — **0->1 only.**
- **LFO never builds** `cutoff_ms_time_02/12/03/32` or their censored arrays
  (`stan/tumor/_lfo_transformed_data.stan` builds only `cutoff_ms_time_01`).
  So the cutoff-censored data plumbing for the other four transitions does not
  exist — this is structural, not a one-line omission.

**Consequence:** the LFO posterior is informed by strictly less data than the full
model; they are not measuring the same thing. Publication config enables all five
transitions (`enable_ms_01..32 = TRUE`), so the intent is clearly all-transition.

**Scope to fix:** build cutoff-censored data + likelihood terms for 0->2, 1->2,
0->3, 3->2 in the LFO model, mirroring the full model's `multistate_lpmf` path but
cutoff-aware. Must reconcile with the detection-week censoring convention used for
0->1 (`_lfo_transformed_data.stan:156-164`).

**Note:** `ms_km_est` is documented (`.claude/rules/sclc.md:50`) as "Multistate
hazard (0→1 transition)" — that's an *endpoint reporting* definition, NOT a
justification for the likelihood ignoring other transitions. Verify this distinction
when fixing.

---

## BLOCKER B2: forecast/background split has an MS index-space inconsistency

The surrogate needs the forecast/background split (`forecast_split_level`,
`forecast_group`) to designate the backgrounded trial. That split has **never run
with split>0** before, and enabling it exposes:

### The live bug
`stan/tumor/sf-ssls-lfo.stan:111`:
`log_cond_surv_01[cutoff_observed_patients]` — `cutoff_observed_patients` holds
FULL patient IDs (1..n_patients, e.g. 497), but `log_cond_surv_01` is sized
`[n_forecast_patients, ...]` (e.g. 78). Out-of-bounds gather when forecast < full.
Worked before only because every prior run had `n_forecast_patients == n_patients`.

### Why the naive fix is wrong (adversarial workflow `wf_8e1d2666`, all-Opus)
- "Widen `log_cond_surv_*` to `n_patients` in the shared module" is **INFEASIBLE**:
  `stan/modules/multistate/transformed_parameters.stan` is included by **8 models**.
  pioneer uses `ms_split_level` with a non-identity `ms_patient_idx` and indexes
  these matrices in forecast-local/compact row space; widening them misaligns its
  likelihood. ms-standalone hardcodes `n_forecast_patients = n_patients`.
- 14 of 22 proposed changes were adversarially **rejected**; the 5 core ones
  (widen shared matrices) were "premise rejected."

### Validated direction (needs human sign-off before code)
- Keep the shared module forecast-sized (untouched).
- Build the full-patient 0->1 conditional survival **tumor-LFO-locally** for the
  likelihood; forecast-restrict the GQ.
- Validated foundation pieces (survived adversarial review): new
  `n_ms_groups_per_level` (MS-specific full group counts, decoupled from the shared
  `n_forecast_groups_per_level`); GQ-only mapping arrays `forecast_patient_mask`,
  `forecast_cutoff_observed_patients`, `patient_to_forecast_cutoff_idx`; initializer
  alignment in BOTH `r/sclc/initializers_fixed.R` AND `r/initializers_ms.R`.
- Full contract + numeric ground-truth test in workflow run `wf_8e1d2666-2a6`
  output (one tracer `trace:trans-03-32` died on an API error → 0->3/3->2 coverage
  thin; re-run before implementing).

### B1 ↔ B2 entanglement
B1 (fix LFO to learn all transitions) and B2 (make MS likelihood background-inclusive
+ forecast-restricted GQ) touch the same code paths. **Fixing B1 first changes the
B2 surface** (B2 would then need full-patient survival for all five transitions, not
just 0->1). They should be designed together, B1 before B2.

---

## Open human decisions
1. **B1 fix scope/approach** — build cutoff-censored all-transition LFO likelihood;
   confirm censoring convention per transition.
2. **B2 KM-denominator intent** — when GQ is forecast-restricted, does the reported
   trial-KM denominator shrink to forecast patients, or keep all cutoff-observed and
   only suppress per-patient curves? (Affects published numbers.)
3. **Sequencing** — B1 then B2, both before the publication surrogate run can produce
   a trustworthy fit. The surrogate validation (synthetic) already stands on its own.

## Clean-up state
- `laplace` store + `/mnt/artifacts/.../laplace/models` wiped; no jobs running.
- Uncommitted: 4 untracked Phase-1 experiment files (`r/experiments/test_laplace_*`,
  `stan/experiments/laplace_tumor_*`) — evidence, keep or remove.
