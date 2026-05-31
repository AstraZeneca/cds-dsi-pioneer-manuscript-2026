# Spop Dropout PFS Convention: PFS-from-OS for 0→3→2 Patients

**Date:** 2026-05-30
**Status:** Approved — ready for implementation planning
**Branch:** `karim/msm-frailty-impl`

## Problem

The superpopulation (spop) PFS Kaplan–Meier curve is pessimistic for the
`died_off_trial` cohort (n=57). After the correlated-frailty fit (run
`202605292133`), their dropout probability rose correctly (17% → 39%), but their
spop median PFS got *worse* (gap −20.5 → −25 wk).

Root cause is **not** the 0→3 hazard or the frailty — it is a definitional
mismatch in how PFS is derived for off-trial deaths:

- **Observed data convention** (`r/data_preparation_pipeline/.../03_data_processing_funcs.R:420`,
  `main_pipeline_data_sclc.R:375`): a patient who dies without recorded
  progression gets `pfs = death_week − 1` — i.e. **death without progression is
  a PFS event dated at death**. For the `died_off_trial` cohort this yields a
  median observed PFS of 42 (their median death week ≈ 41), even though their
  last on-trial assessment was week ≈ 21.
- **Model spop convention** (`stan/pfs.stanfunctions` `derive_spop_pfs`, cause 3,
  lines 1393–1401): a dropout (0→3) has its combined PFS **right-censored at the
  dropout time `t_03`** — the 3→2 sojourn never enters PFS. So the model censors
  these patients at ~week 17 while the observed KM credits a PFS *event* at ~week 41.

The −25 wk gap is therefore largely the model faithfully encoding a *different*
PFS definition than the data uses. The fix is to align the spop PFS derivation
with the observed convention.

## Goal

Make the spop PFS for a dropout patient follow the observed-data convention:
- If the patient goes on to die via 3→2 within the simulation window → **PFS
  event at the 3→2 death week**.
- If the patient does not die within the window → **PFS right-censored at the
  forecast horizon** (`max_all_t`).

Apply to **both** the combined PFS (`spop_pfs`) and the multistate-only variant
(`spop_ms_pfs`).

## Design

### Section 1 — Core semantic change (approach A: reorder + graft)

For a spop patient whose first state-0 exit is dropout (`spop_cause == 3`), the
combined PFS and ms-PFS both become the OS outcome:

- 3→2 death observed within window → PFS event at death week (`spop_os`, censored 0)
- 3→2 not reached by `max_all_t` → PFS censored at horizon (`spop_os`, censored 1)

**Mechanism — approach A (reorder + graft):**

Today the GQ per-patient loop (`stan/pfs.stanfunctions` ~1958–1980) calls
`derive_spop_pfs` *before* `derive_spop_os_rng`, so the PFS derivation cannot see
the 3→2 outcome. Reorder so the OS draw runs first, then for `spop_cause == 3`
overwrite:

```
spop_pfs[j]              = spop_os[j];
spop_right_censored[j]   = spop_os_censored[j];
spop_ms_pfs[j]           = spop_os[j];
spop_ms_right_censored[j]= spop_os_censored[j];
```

Non-dropout causes keep `derive_spop_pfs` unchanged.

**Why approach A over passing the sojourn into `derive_spop_pfs` (approach B):**
the OS outcome *is* exactly the (time, censored) pair the new PFS rule needs.
Reusing it guarantees **PFS == OS for dropouts by construction**, preserving the
PFS ≤ OS invariant. Approach B would re-derive the death time and risk drift
between the OS death week and the PFS death week. The reorder is safe because
`derive_spop_os_rng` depends only on `spop_cause`, `spop_exit`, and the survival
rows — all available before the PFS call.

**Censoring horizon decision:** a non-dying dropout is censored at the OS
forecast horizon (`max_all_t`), inherited verbatim from `derive_spop_os_rng` /
`sample_dropout_death_rng`. This is consistent with the spop's purpose — a clean
re-roll of the fitted process from week 1 — rather than importing the observed
cohort's administrative follow-up into the counterfactual.

### Section 2 — CIF decomposition consequence

The CIF cause-attribution (Stan GQ
`_endpoints_generated_quantities.stan:336–347`, mirrored in `r/util.R:452–459`)
currently infers cause from the censoring pattern:

```r
pfs_e <- sms_rc == 0L
dd    <- pfs_e & s_osc == 0L & sms_pfs == s_os   # 0→2
prog  <- pfs_e & !dd                             # 0→1
drop_ <- !pfs_e & sms_pfs <= max_all_t           # 0→3 (inferred from censoring)
```

Under the new rule a dropout-death has `sms_rc == 0` and `sms_pfs == s_os`, so it
would satisfy `dd` and be **misattributed from CIF_03 into CIF_02**.

**Fix (option 2 — chosen):** key the CIF cause-attribution off the explicit
`spop_is_dropout` flag rather than the censoring pattern. `spop_is_dropout[j]` is
already set from `spop_cause == 3` (pfs.stanfunctions:1956), independent of
PFS event/censoring status, and `util.R` already documents that `is_dropout` must
be set explicitly from the upstream cause. Re-express the predicates so:

- `drop_` = `is_dropout` (all 0→3 exits, dead or alive)
- `dd` = death-without-progression that is **not** a dropout (on-trial 0→2)
- `prog` = 0→1 progression (unchanged)

This keeps the three CIF curves meaning what they did (01 = progression,
02 = on-trial death, 03 = dropout regardless of eventual death) while PFS treats
the dropout-death as an event. Rejected option 1 (let dropout-deaths fold into
CIF_02) because it would empty out the CIF_03 curve on the website.

### Section 3 — Touch-points & consistency

Four coupled sites, changed in lockstep:

1. **`stan/pfs.stanfunctions` GQ loop (~1958–1980):** reorder OS-before-PFS;
   graft `spop_os`/`spop_os_censored` into PFS + ms-PFS for `spop_cause == 3`.
2. **`stan/.../_endpoints_generated_quantities.stan` CIF block (336–347):**
   `drop_`/`dd` predicates key off `spop_is_dropout`.
3. **`r/util.R` CIF mirror (433–459):** identical predicate change; must stay
   bit-identical to the Stan GQ (`compute_cif_from_draws` exists to mirror it).
4. **`derive_spop_pfs` (pfs.stanfunctions:1377):** the cause-3 branch becomes
   superseded for the spop path once the graft overwrites it. Leave it returning
   the current value with a comment noting it is overridden downstream; do not
   alter the sample-path symmetry.

**Sample (conditional) path is untouched.** `derive_sample_os_rng` (lines
1614–1618) already handles observed off-trial death exactly via real
`ms_time_03`/`ms_time_32`. The new rule applies only to the spop re-roll, which
is the curve compared against observed PFS.

### Section 4 — Validation & cost

- **No re-fit.** This is a generated-quantities + downstream-R change. The
  posterior draws are unchanged; only the endpoint *derivation* changes.
  Regenerate downstream targets (GQ → KM rvars → CIF), not the MCMC.
- **Build sequence:** Stan syntax check (`stanc`) → regenerate
  `tumor_ssls_draws_endpoints` / `km_rvar` / CIF targets → verify
  `died_off_trial` spop median PFS moves toward 42 → confirm CIF_03 still
  populated.
- **Stan↔R bit-identity** is the primary test anchor: sites 2 and 3 must produce
  identical CIF curves.
- **Success criterion:** `died_off_trial` spop median PFS rises from 17 toward
  the observed 42; overall pessimism gap shrinks; CIF_03 retains the dropout
  cohort.

## Out of scope

- 3→2 sojourn-hazard re-calibration. The sojourn timing now drives PFS-event
  timing for dropouts; if the gap only partially closes, that is the *next*
  separable lever — not part of this change.
- Re-fitting the MCMC.
- Any change to the observed-data PFS derivation (`death_week − 1` convention is
  the target we align to, not change).
- Sample/conditional-path PFS or OS.

## Success criterion

`died_off_trial` spop median PFS rises materially from 17 toward observed 42,
the overall spop-vs-observed PFS pessimism shrinks, the CIF_03 (dropout) curve
remains populated, and Stan and R CIF recomputations remain bit-identical.
