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
rows — all available before the PFS call. **RNG-stream stability:**
`derive_spop_pfs` is RNG-free (documented at `pfs.stanfunctions:1364`), so moving
it after the OS draw does not perturb the RNG sequence for any patient — only
cause-3 patients' final `spop_pfs` values change, and only via the graft.

**Retain the dropout time `t_03`.** The graft overwrites the PFS slot, but the
3→2 sojourn KM (Section 3, site 4) needs the dropout exit time. Capture
`spop_time_03` (the cause-3 exit time, already in scope at the graft site) into a
per-patient `spop_dropout_week[j]` array before/at the graft so downstream code
can recover the sojourn as `spop_os − spop_dropout_week`.

**Censoring horizon decision:** a non-dying dropout is censored at the OS
forecast horizon (`max_all_t`), inherited verbatim from `derive_spop_os_rng` /
`sample_dropout_death_rng`. This is consistent with the spop's purpose — a clean
re-roll of the fitted process from week 1 — rather than importing the observed
cohort's administrative follow-up into the counterfactual. **Caveat (validated in
Section 4):** observed non-dying dropouts are censored at their last contact
(~wk 21) while the spop censors at `max_all_t`, so the spop dropout arm stays in
the risk set longer. This can lift the PFS *tail* optimistically even when the
median moves for the right reason — validation must check the censoring-time
distribution, not just the median.

**Degenerate case.** When `enable_ms_32 == 0`, `derive_spop_os_rng` for cause 3
returns `(t_03, censored = 1)`, so the graft sets `spop_pfs = t_03` censored —
identical to today's behavior. The PFS-from-OS rule only changes anything when
3→2 is modeled (it is, in the publication config).

### Section 2 — CIF decomposition consequence

The live Stan CIF is the function **`compute_trial_cif`**
(`stan/pfs.stanfunctions:1536`), invoked at
`stan/tumor/_tumor_endpoints_generated_quantities.stan:472` (spop) and `:479`
(sample). It already keys the dropout count off the explicit `spop_is_dropout`
array and operates on the **combined** `spop_pfs`/`spop_os` (not `spop_ms_pfs`).
Its branch structure is:

```stan
if (!right_censored[j]) {
  if (!os_censored[j] && pfs[j] == os[j]) cnt_02 += 1;   // 0→2 direct death
  else                                    cnt_01 += 1;   // 0→1 progression
} else if (is_dropout[j]) {
  cnt_03 += 1;                                            // 0→3 dropout
} // else admin-censored
```

**The problem the graft introduces:** after the graft a dropout-death has
`right_censored == 0` and `pfs == os`, so it enters the **first** branch and is
counted as `cnt_02` — it never reaches the `is_dropout` branch. CIF_03 would lose
every dropout that dies in-window (the exact n=57 cohort). Merely re-keying a
predicate is insufficient; the **branch must be structurally reordered** so the
`is_dropout` test runs first:

```stan
if (is_dropout[j]) {
  cnt_03 += 1;                                            // 0→3 dropout (dead or alive)
} else if (!right_censored[j]) {
  if (!os_censored[j] && pfs[j] == os[j]) cnt_02 += 1;    // on-trial 0→2 death
  else                                    cnt_01 += 1;    // 0→1 progression
} // else admin-censored
```

Update the `compute_trial_cif` docstring (`pfs.stanfunctions:1515–1523`) to match.
This keeps the three curves meaning what they did (01 = progression, 02 =
on-trial death, 03 = dropout regardless of eventual death) while PFS treats the
dropout-death as an event. (Rejected: letting dropout-deaths fold into CIF_02,
which would empty the CIF_03 curve.)

**No R mirror.** The former `r/util.R::compute_cif_from_draws` (a dead "TEMPORARY
workaround" with two identical copies and zero callers) was deleted in commit
`629d245f`. The pipeline reads `spop_cif_*`/`sample_cif_*` directly from the Stan
fit, so `compute_trial_cif` is now the **single source of truth** for the CIF —
no Stan↔R bit-identity requirement.

### Section 3 — Touch-points & consistency

Three coupled sites, changed in lockstep:

1. **`stan/pfs.stanfunctions` GQ loop (~1958–1980):** reorder so
   `derive_spop_os_rng` runs before the PFS assignment; for `spop_cause == 3`,
   graft `spop_os`/`spop_os_censored` into `spop_pfs`/`spop_right_censored` and
   `spop_ms_pfs`/`spop_ms_right_censored`. Also capture `spop_time_03` into
   `spop_dropout_week[j]` for the sojourn KM.
2. **`compute_trial_cif` (`stan/pfs.stanfunctions:1536`):** structural branch
   reorder — `is_dropout` test first (see Section 2). Update docstring.
3. **`spop_km_32` sojourn KM block
   (`stan/tumor/_tumor_endpoints_generated_quantities.stan:436–450`):** currently
   computes the sojourn as `soj = max(1, spop_os[i] − spop_pfs[i])` (line 443).
   After the graft `spop_os == spop_pfs` for dropouts, collapsing this to
   `max(1, 0) = 1` and destroying the sojourn KM. Recompute from the dropout
   week: `soj = max(1, spop_os[i] − spop_dropout_week[i])`.
4. **`spop_km_12` sojourn KM block
   (`stan/tumor/_tumor_endpoints_generated_quantities.stan:404–418`):** the
   1→2 (post-progression) sojourn filters on `!spop_ms_right_censored[i]`
   (line 410). After the graft flips `spop_ms_right_censored` to 0 for
   dropout-deaths, those patients would leak into the 1→2 progression-sojourn
   curve (with a spurious `soj = max(1, spop_os − spop_ms_pfs) = 1`). Exclude
   dropouts: add `&& !spop_is_dropout[i]` to the filter so only genuine 1→2
   progressors contribute. (Surfaced during planning; not in the original review.)

**`derive_spop_pfs` cause-3 branch is unreachable, not "superseded".** The inner
SLD-PD branch (`pfs.stanfunctions:1395–1397`) can never fire for a cause-3
patient: `classify_spop_exit` (`1344–1346`) only returns cause 3 when no
uncensored target PD precedes `t_03` (otherwise progression wins the priority
order at `1348`). So `t_target ≤ t_03` with `c_target == 0` cannot co-occur with
`cause == 3`. The graft therefore loses no live earlier-progression event. Either
delete the unreachable inner if/else (return the censored-at-dropout tuple) or
comment it accurately as "unreachable for spop". Do **not** alter the sample-path
symmetry.

**Sample (conditional) path is untouched.** `derive_sample_os_rng`
(`multistate.stanfunctions:1614–1624`) already handles observed off-trial death
exactly via real `ms_time_03`/`ms_time_32`, and `derive_sample_pfs`
(`pfs.stanfunctions:1491–1497`) keeps the current logic. The new rule applies
only to the spop re-roll — the curve compared against observed PFS.

### Section 4 — Validation & cost

- **No re-fit.** Generated-quantities + downstream change only. The posterior
  draws are unchanged; only the endpoint *derivation* changes. Regenerate
  downstream targets (GQ → KM rvars → CIF), not the MCMC.
- **Build sequence:** Stan syntax check (`stanc`) → regenerate
  `tumor_ssls_draws_endpoints` / `km_rvar` / CIF / sojourn-KM targets → checks below.
- **Validation (beyond the median):**
  1. `died_off_trial` spop median PFS rises from 17 toward observed 42.
  2. **CIF_03 count == n_dropout per trial**, and `cif_01 + cif_02 + cif_03 +
     admin-censored == n` — directly catches the dropout-death leak from Section 2.
  3. **PFS ≤ OS pointwise** across the whole spop curve (not just dropouts).
  4. **`spop_km_32` still meaningful** (not collapsed to a point mass at 1) after
     the site-3 fix.
  5. **Tail check:** report the spop dropout censoring-time distribution vs the
     observed last-contact distribution, and compare PFS at horizon / RMST against
     the observed `btype == "ub"` curve — to surface horizon-vs-last-contact tail
     optimism.

### Section 5 — Documentation of the PFS estimand

The model's PFS for off-trial deaths follows **SCLC's operational PFS
definition** — death-without-progression is a PFS *event* at the death week (via
the 0→3→2 path) — **not** canonical RECIST/FDA PFS, which would censor a death
that follows a long lost-to-follow-up gap at last assessment. These differ
precisely for off-trial deaths (death at ~wk 41 after dropout at ~wk 21). This is
a deliberate estimand choice to keep the model and the observed comparator on the
same PFS definition. Record it explicitly in the clinical-endpoints model spec
(`quarto/.../documentation/clinical-endpoints-specification.qmd`) as part of this
change, so readers know what the model's PFS means.

## Out of scope

- 3→2 sojourn-hazard re-calibration. The sojourn timing now drives PFS-event
  timing for dropouts; if the gap only partially closes, that is the *next*
  separable lever — not part of this change.
- Re-fitting the MCMC.
- Any change to the observed-data PFS derivation (`death_week − 1` convention is
  the target we align to, not change).
- Sample/conditional-path PFS or OS.
- Adopting canonical censor-after-gap PFS (the rejected estimand; would invert the fix).

## Success criterion

`died_off_trial` spop median PFS rises materially from 17 toward observed 42; the
overall spop-vs-observed PFS pessimism shrinks; CIF_03 retains the dropout cohort
(count == n_dropout); `spop_km_32` remains a genuine sojourn curve; PFS ≤ OS holds
pointwise; and the PFS estimand is documented in the model spec.
