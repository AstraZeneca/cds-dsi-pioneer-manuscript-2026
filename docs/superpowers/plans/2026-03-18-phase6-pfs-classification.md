# Plan: Phase 6 — PFS Classification and Derivation Functions

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/pfs.stanfunctions` (lines 1122–1316)
**Target assertions:** ~65 new assertions

## Background

The classification-first competing-risks pipeline (Issue 79) consists of four pure functions
(`classify_spop_exit`, `derive_spop_pfs`, `classify_sample_exit`, `derive_sample_pfs`) plus
`assessment_gated_survival_time_rng`. The first four are pure derivations with no RNG — they
encode the tie-breaking logic that determines PFS/OS endpoints from competing-risks samples.

The RNG function `assessment_gated_survival_time_rng` accumulates per-visit cumulative hazard
instead of sampling per-week. Its determinism can be verified by passing a survival probability
of exactly 1.0 (always survive) or 0.0 (always fail) — these collapse to deterministic outcomes.

## Functions to cover

| Function | Type |
|----------|------|
| `classify_spop_exit` | Pure: `(t_01,c_01,t_02,c_02,t_03,c_03,enable_02,enable_03) → (cause, exit_time)` |
| `derive_spop_pfs` | Pure: `(cause, exit_time, ...) → (pfs, rc, ms_pfs, ms_rc)` |
| `classify_sample_exit` | Pure: same structure with `ms_final_state` routing |
| `derive_sample_pfs` | Pure: same 4-tuple output |
| `assessment_gated_survival_time_rng` | RNG (2 overloads): tested with boundary survival probabilities |

---

## Task 1 — R helpers `helper-pfs-classify.R`

**File:** `tests/testthat/helper-pfs-classify.R`

```r
r_classify_spop_exit <- function(t_01, c_01, t_02, c_02, t_03, c_03,
                                  enable_ms_02, enable_ms_03) {
  progressed_first <- !c_01 && (c_02 || t_01 <= t_02) && (c_03 || t_01 <= t_03)
  died_directly    <- enable_ms_02 && !c_02 && (c_01 || t_02 < t_01) && (c_03 || t_02 <= t_03)
  dropped_out      <- enable_ms_03 && !c_03 && (c_01 || t_03 < t_01) && (c_02 || t_03 < t_02)
  if (progressed_first) list(cause = 1L, exit_time = t_01)
  else if (died_directly) list(cause = 2L, exit_time = t_02)
  else if (dropped_out)   list(cause = 3L, exit_time = t_03)
  else                    list(cause = 0L, exit_time = max(t_01, t_02, t_03))
}

r_derive_spop_pfs <- function(cause, exit_time, t_target, c_target, t_ms01, c_ms01,
                               t_03, enable_ms_02, enable_ms_03) {
  if (cause == 1L) return(c(exit_time, 0L, exit_time, 0L))
  if (cause == 2L) return(c(exit_time, 0L, exit_time, 0L))
  if (cause == 3L) {
    if (!c_target && t_target <= t_03) return(c(t_target, 0L, t_03, 1L))
    else                               return(c(t_03, 1L, t_03, 1L))
  }
  # cause == 0
  c(min(t_target, t_ms01), 1L, t_ms01, c_ms01)
}

r_classify_sample_exit <- function(st_01, sc_01, st_02, sc_02,
                                    ms_final_state, ms_censored_02, ms_time_03,
                                    enable_ms_02, enable_ms_03) {
  dropped_out    <- enable_ms_03 && ms_final_state == 3L
  died_directly  <- !dropped_out && enable_ms_02 && !sc_02 &&
                    (!ms_censored_02 || sc_01 || st_02 < st_01)
  progressed_first <- !dropped_out && !died_directly && !sc_01 &&
                      (sc_02 || st_01 <= st_02)
  if (dropped_out)       list(cause = 3L, exit_time = ms_time_03)
  else if (died_directly)  list(cause = 2L, exit_time = st_02)
  else if (progressed_first) list(cause = 1L, exit_time = st_01)
  else                       list(cause = 0L, exit_time = max(st_01, st_02))
}

r_derive_sample_pfs <- function(cause, exit_time, t_target, c_target,
                                 t_ms_raw, c_ms_raw, ms_time_03,
                                 enable_ms_02, enable_ms_03) {
  if (cause == 2L) return(c(exit_time, 0L, exit_time, 0L))
  if (cause == 3L) {
    if (!c_target && t_target <= ms_time_03) return(c(t_target, 0L, ms_time_03, 1L))
    else                                      return(c(ms_time_03, 1L, ms_time_03, 1L))
  }
  rc <- if (cause == 0L) 1L else 0L
  c(exit_time, rc, t_ms_raw, c_ms_raw)
}
```

**Test file:** `tests/testthat/test-helper-pfs-classify.R`
**Assertions (~18):**
- `classify_spop_exit`: progression wins; direct death wins; tie → progression priority
- `classify_spop_exit`: both censored → cause=0, exit_time=max(t_01,t_02,t_03)
- `derive_spop_pfs`: cause=1 → both PFS uncensored; cause=3 with early SLD → separate ms/pfs
- `classify_sample_exit`: ms_final_state=3 always routes to dropout regardless of forecasts
- `derive_sample_pfs`: cause=2 → both PFS uncensored at t_02

---

## Task 2 — Stan harness `test_pfs_classify_all.stan`

**File:** `tests/testthat/stan/test_pfs_classify_all.stan`

```stan
functions {
  #include pfs.stanfunctions
  #include util.stanfunctions
  #include pos.stanfunctions
}
data {
  int<lower=1> N_CASES;
  // classify_spop_exit inputs
  array[N_CASES] int t_01;  array[N_CASES] int c_01;
  array[N_CASES] int t_02;  array[N_CASES] int c_02;
  array[N_CASES] int t_03;  array[N_CASES] int c_03;
  array[N_CASES] int enable_02;
  array[N_CASES] int enable_03;
  // derive_spop_pfs additional inputs
  array[N_CASES] int t_target;  array[N_CASES] int c_target;
  array[N_CASES] int t_ms01;    array[N_CASES] int c_ms01;
  // classify_sample_exit inputs
  array[N_CASES] int ms_final_state;
  array[N_CASES] int ms_censored_02;
  array[N_CASES] int ms_time_03;
  // derive_sample_pfs additional inputs
  array[N_CASES] int t_ms_raw;  array[N_CASES] int c_ms_raw;
  // assessment_gated_survival_time_rng
  int<lower=1> MAX_T_RNG;
  int N_VISITS_RNG;
  array[N_VISITS_RNG] int visit_weeks;
  int max_all_t_rng;
}
generated quantities {
  array[N_CASES] int cause_spop;
  array[N_CASES] int exit_time_spop;
  array[N_CASES] int pfs_spop;       array[N_CASES] int rc_spop;
  array[N_CASES] int ms_pfs_spop;    array[N_CASES] int ms_rc_spop;
  array[N_CASES] int cause_sample;
  array[N_CASES] int exit_time_sample;
  array[N_CASES] int pfs_sample;     array[N_CASES] int rc_sample;
  array[N_CASES] int ms_pfs_sample;  array[N_CASES] int ms_rc_sample;

  for (i in 1:N_CASES) {
    (cause_spop[i], exit_time_spop[i]) = classify_spop_exit(
      t_01[i], c_01[i], t_02[i], c_02[i], t_03[i], c_03[i], enable_02[i], enable_03[i]
    );
    (pfs_spop[i], rc_spop[i], ms_pfs_spop[i], ms_rc_spop[i]) = derive_spop_pfs(
      cause_spop[i], exit_time_spop[i],
      t_target[i], c_target[i], t_ms01[i], c_ms01[i], t_03[i],
      enable_02[i], enable_03[i]
    );
    (cause_sample[i], exit_time_sample[i]) = classify_sample_exit(
      t_01[i], c_01[i], t_02[i], c_02[i],
      ms_final_state[i], ms_censored_02[i], ms_time_03[i],
      enable_02[i], enable_03[i]
    );
    (pfs_sample[i], rc_sample[i], ms_pfs_sample[i], ms_rc_sample[i]) = derive_sample_pfs(
      cause_sample[i], exit_time_sample[i],
      t_target[i], c_target[i], t_ms_raw[i], c_ms_raw[i], ms_time_03[i],
      enable_02[i], enable_03[i]
    );
  }

  // assessment_gated_survival_time_rng: all-zero log_surv → always event at first visit
  row_vector[MAX_T_RNG] log_surv_always_event = rep_row_vector(log(0.0), MAX_T_RNG);
  int ags_event; int ags_censored;
  (ags_event, ags_censored) = assessment_gated_survival_time_rng(
    log_surv_always_event, visit_weeks, N_VISITS_RNG, max_all_t_rng
  );
  // always survive → censored
  row_vector[MAX_T_RNG] log_surv_always_surv = rep_row_vector(0.0, MAX_T_RNG);
  int ags_cens_event; int ags_cens_flag;
  (ags_cens_event, ags_cens_flag) = assessment_gated_survival_time_rng(
    log_surv_always_surv, visit_weeks, N_VISITS_RNG, max_all_t_rng
  );
}
```

**Test file:** `tests/testthat/test-stan-pfs-classify.R`

Key test cases:

```r
# Case 1: progression wins (t_01=5 < t_02=8, both uncensored)
# Case 2: direct death wins (t_02=3 < t_01=7, 02 uncensored)
# Case 3: both censored → cause=0
# Case 4: dropout (ms_final_state=3) always wins in sample pathway
# Case 5: cause=2 → both PFS=t_02, uncensored in derive_spop_pfs

# assessment_gated_survival_time_rng:
# log_surv = log(0) → P(survive interval) = 0 → event at first visit
# log_surv = 0 → P(survive) = 1 → censored at max_all_t
```

**Assertion count:** ~27 Stan assertions

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-pfs-classify.R` + `test-helper-pfs-classify.R` | ~18 R |
| 2 | `test_pfs_classify_all.stan` + `test-stan-pfs-classify.R` | ~27 Stan |
| **Total** | | **~45 new assertions** |
