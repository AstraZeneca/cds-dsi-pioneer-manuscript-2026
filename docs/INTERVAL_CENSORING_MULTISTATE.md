# Interval Censoring Marginalization: Multistate Model

> **Status:** Implemented (semi-Markov only). Markov and Extended time scales
> fall back to no-IC (`gap` is forced to 0 when `ms_time_scale_12 != 1`).
> See GitHub issue #82 for tracking.

## Purpose: Why Interval Censoring Marginalization?

In oncology trials, disease progression is typically detected through periodic tumor assessments (e.g., imaging every 6–12 weeks). We rarely know the exact week a patient progressed; instead, we observe:

- **T_c** = week of last clean assessment (RECIST stable or response)
- **T_d** = week of first positive assessment (RECIST progressive disease)

The true progression week **s** lies somewhere in the open interval **(T_c, T_d]**. A naive approach would assume progression occurred at T_d (the first detection), but this underestimates the true progression risk and introduces bias, especially when assessment intervals are long.

**Interval censoring marginalization** solves this by treating s as a latent variable and marginalizing over all plausible values:

$$\mathcal{L} = \sum_{s=T_c+1}^{T_d} P(\text{survive to s, progress at s, observe next state-1/2/3 progression})$$

This accounts for all possible true progression weeks consistent with the observed assessment pattern, giving a more accurate likelihood and posterior inference.

## Mathematical Formulation

For a patient with progression detected in the interval (T_c, T_d], the interval-censored 0→1 (progression) log-likelihood is:

$$\log L = \log \Phi + \log \text{surv}_{01}(1 \ldots T_c) + \log \sum_{s=T_c+1}^{T_d} \alpha(s) + \log \left( \sum_{s=T_c+1}^{T_d} \beta(s) \right)$$

where:

- **Φ** = base survival (0→2 and 0→3 transitions to non-progression endpoints)
- **surv_{01}(1…T_c)** = survival from week 1 to T_c *not* progressing via 0→1
- **α(s)** = progression component: $\text{surv}_{01}(T_c+1 \ldots s-1) + \log h_{01}(s)$
  - Hazard h_{01}(s) is the 0→1 transition hazard at week s
- **β(s)** = sojourn component: survival in state 1 (post-progression) until the next event
  - For patients progressing to state 2 (death post-progression): $\text{surv}_{12}(\tau_s)$ where $\tau_s = T_{\text{end}} - s$
  - For patients progressing to state 3 (other event): observed at visit weeks only

### Incremental Sojourn Update Trick

Computing the full marginalization naively requires iterating over every possible s ∈ {T_c+1, …, T_d} and recalculating the sojourn survival log-likelihood each time, which is expensive.

The **incremental sojourn update** optimizes this by computing the sojourn term once and updating it incrementally as s steps backward from T_d:

1. **Seed** (at s = T_d):
   - Start with the full sojourn survival from s=T_d to T_end (kk=0 in the algorithm)

2. **Step backward** (s = T_d−1, T_d−2, …, T_c+1):
   - Remove one week of sojourn survival at the *beginning* of the state-1 period
   - For state-2 progressors: update the state-1→2 hazard and survival terms
   - For state-3 progressors: only update at observed visit weeks

This reduces computation from O(T_d − T_c) evaluations of the full sojourn survival to O(1) amortized cost per s.

## Implementation Details

### ms_ic_gap_01 Derivation (transformed_data.stan)

In `stan/modules/multistate/transformed_data.stan`, the IC gap is computed for each patient:

```stan
for (i in 1:n_patients) {
  ms_ic_gap_01[i] = enable_ms_01 && ms_0_1_event[i]
                    ? interval_censored[i] + 1
                    : 0;
}
```

- **If progression detected** (`ms_0_1_event[i] = 1`):
  - `ms_ic_gap_01[i] = interval_censored[i] + 1`
  - `interval_censored[i]` is the input gap = T_d − T_c in weeks
  - `+1` accounts for the inclusive loop range: s ∈ {T_c+1, …, T_d}

- **If censored or deterministic progression** (no interval censoring):
  - `ms_ic_gap_01[i] = 0`
  - This flag disables IC marginalization, recovering the no-IC likelihood

### Time Scale Scope: Semi-Markov Only

**Fully implemented:**
- **Semi-Markov** (ms_time_scale_12 = 1): State-1 sojourn time is measured from s (weeks since progression), independent of the 0→1 hazard history. IC marginalization integrates cleanly.

**Falls back to gap=0 (no IC):**
- **Markov** (ms_time_scale_12 = 0): Sojourn time measured from calendar week; requires re-indexing the sojourn hazard for each s. Complex bookkeeping.
- **Extended** (ms_time_scale_12 = 2): Shared time scale with baseline hazard; requires careful handling of age/time/epoch interactions.

Both Markov and Extended are flagged with a fallback in `transformed_data.stan`:

```stan
if (!is_semi_markov) {
  // Markov or Extended time scale: disable IC for now
  for (i in 1:n_patients) {
    if (ms_ic_gap_01[i] > 0) {
      ms_ic_gap_01[i] = 0;  // Fall back to no-IC
    }
  }
}
```

See GitHub issue #82 for tracking these extensions.

## Testing

Three unit tests in `tests/testthat/test-stan-multistate-ic.R` validate the IC marginalization:

1. **test_multistate_ic_gap_zero_identity()**
   - Verifies that gap=0 recovers the no-IC likelihood
   - Ensures no regression when IC is disabled (e.g., censored patients)

2. **test_multistate_ic_hand_computed_marginalization()**
   - Hand-computes the marginalized likelihood for a toy scenario
   - Compares against the Stan model output
   - Validates the incremental sojourn update trick

3. **test_multistate_ic_prog_deterministic_bypass()**
   - Confirms that deterministic progressors (gap unknown) correctly bypass IC
   - Ensures the fallback flag logic works as intended

These tests run against the standalone Stan model in `tests/testthat/stan/test_multistate_ic.stan`, which isolates the IC logic for fast validation.

## Key Design Decisions

1. **Marginalization over all weeks, not just mean**: We don't approximate s ≈ (T_c + T_d) / 2. Instead, we integrate over the full posterior distribution of s.

2. **log_sum_exp for numerical stability**: The sum of exponentials is computed using `log_sum_exp()` to avoid overflow.

3. **State-specific next events**: The sojourn survival depends on which state the patient enters after progression (state 2 = death with progression, state 3 = other event). The likelihood must condition on the observed next event.

4. **Visit-week granularity for state-3 progressors**: State-3 (other events) are observed only at visit weeks, not continuously. The IC marginalization only sums over s values corresponding to potential visit weeks before the observation week.

## Related Files

- **Model code**: `stan/modules/multistate/transformed_data.stan` (gap computation), `stan/modules/multistate/log_lik.stan` (likelihood evaluation)
- **Tests**: `tests/testthat/test-stan-multistate-ic.R`, `tests/testthat/stan/test_multistate_ic.stan`
- **Targets pipeline**: `targets/sclc_targets.R` — no special handling needed; IC is automatic when `enable_ms_01 = TRUE`
