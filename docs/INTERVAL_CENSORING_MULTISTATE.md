# Interval Censoring for the Multistate Model

**Related**: [GitHub Issue #81 (sclc repo)](https://github.com/azu-oncology-rd/cds-dsi-pioneer-sclc-01-2025/issues/81)

## 1. Problem Statement

The multistate illness-death model estimates hazards for transitions between
states (alive/progression-free → progressed → dead). The key challenge is that
**progression is interval-censored**: the biological event occurs at an unknown
time between the last clean assessment and the detection visit.

After removing interval censoring (IC) from the SLD measurement model (where it
was not needed — the state-space model naturally handles observation gaps), the
multistate model inherited the same "place event at detection week" convention.
This document describes why IC must be restored for the multistate model
specifically, and derives the marginalized likelihood.

### Why IC matters for multistate but not SLD

| Aspect | SLD Model | Multistate Model |
|--------|-----------|-----------------|
| **Nature** | Continuous latent state | Discrete transition event |
| **Between visits** | State-space model propagates the latent trajectory | Event could have occurred at any week |
| **Observation model** | Noisy measurement at assessment | Binary detection (event or no event) |
| **IC needed?** | No — measurement model is sufficient | **Yes** — event time is uncertain |

### Bias from placing events at detection week

When event times are set to the detection week $T_d$, the model asserts that the
patient survived in state 0 through weeks $T_c+1, \ldots, T_d - 1$ (where $T_c$
is the last clean assessment). This:

1. **Inflates survival** in weeks just before common assessment times
2. **Creates GP hazard artifacts** — the GP learns spurious peaks at common
   assessment weeks (6, 12, 18, 24...) because that is where all "events" cluster
3. **Makes predictions schedule-dependent** — the learned hazard shape reflects
   the assessment schedule rather than biology

## 2. Notation

| Symbol | Description |
|--------|-------------|
| $T_c^{(i)}$ | Last clean assessment week for patient $i$ (no progression detected) |
| $T_d^{(i)}$ | Detection week for patient $i$ (progression first detected) |
| $T_\dagger^{(i)}$ | Death week for patient $i$ (exact, not interval-censored) |
| $h_{jk}(i, t)$ | Discrete-time cause-specific hazard for transition $j \to k$, patient $i$, week $t$ |
| $\ell_{jk}(i, t)$ | Log conditional survival: $\log(1 - h_{jk}(i, t))$ — stored as `log_cond_surv_jk[i, t]` |
| $S_{jk}(i, a, b)$ | Interval survival: $\prod_{t=a}^{b}(1 - h_{jk}(i, t)) = \exp\!\bigl(\sum_{t=a}^{b} \ell_{jk}(i, t)\bigr)$ |
| $\Delta^{(i)}$ | IC gap: $T_d^{(i)} - T_c^{(i)}$ (typically 4–6 weeks on a 6-week cycle) |

### Transitions and interval censoring status

| Transition | Event type | IC? | Reason |
|------------|-----------|-----|--------|
| $0 \to 1$ | Progression | **Yes** | Detected at assessment, occurred earlier |
| $0 \to 2$ | Death without PD | No | Death time observed exactly |
| $0 \to 3$ | Dropout | No | Administrative event at last visit |
| $1 \to 2$ | Post-PD death | No | Death time exact, but **sojourn start is uncertain** |
| $3 \to 2$ | Off-trial death | No | Death time exact, sojourn start known |

The $1 \to 2$ sojourn start depends on the true (unknown) $0 \to 1$ transition
time, so IC in the $0 \to 1$ transition **cascades** to the $1 \to 2$ sojourn.

## 3. Current Likelihood (No IC)

### State 1: Progressed, alive

Patient progressed at detection week $T_d$, censored for $1 \to 2$ at time
$T_{\text{cens}}$, with sojourn $\tau = T_{\text{cens}} - T_d$.

$$
\log L_{\text{current}}^{(1)} =
  \underbrace{\sum_{t=1}^{T_d - 1} \ell_{01}(t)}_{\text{survived } 0 \to 1}
+ \underbrace{\log\!\bigl(1 - e^{\ell_{01}(T_d)}\bigr)}_{\text{event at } T_d}
+ \underbrace{\sum_{t=1}^{T_d - 1} \ell_{02}(t)}_{\text{survived } 0 \to 2}
+ \underbrace{\sum_{t=1}^{T_d - 1} \ell_{03}(t)}_{\text{survived } 0 \to 3}
+ \underbrace{\mathcal{L}_{12}(\tau)}_{\text{sojourn (censored)}}
$$

> **Note**: The $0 \to 1$ event hazard is only included for stochastic
> progressions (`prog_deterministic = 0`). For RECIST-determined PD, the SLD
> model pins the event time and the hazard term is omitted.

### State 2 via $0 \to 1 \to 2$: Progressed then died

Patient progressed at $T_d$, died at $T_\dagger$, sojourn $\tau = T_\dagger - T_d$.

$$
\log L_{\text{current}}^{(2 \text{ via } 1)} =
  \sum_{t=1}^{T_d - 1} \bigl[\ell_{01}(t) + \ell_{02}(t) + \ell_{03}(t)\bigr]
+ \log\!\bigl(1 - e^{\ell_{01}(T_d)}\bigr)
+ \mathcal{L}_{12}(\tau, \text{event})
$$

### Sojourn likelihood $\mathcal{L}_{12}$

For semi-Markov (`ms_time_scale_12 = 1`), with sojourn time $\tau$ indexed from 1:

$$
\mathcal{L}_{12}(\tau, \text{censored}) = \sum_{u=1}^{\tau} \ell_{12}^{(s)}(u)
$$

$$
\mathcal{L}_{12}(\tau, \text{event}) = \sum_{u=1}^{\tau - 1} \ell_{12}^{(s)}(u) + \log\!\bigl(1 - e^{\ell_{12}^{(s)}(\tau)}\bigr)
$$

For Markov and extended time scales, the clock-forward component
$\ell_{12}^{(t)}$ is evaluated at calendar weeks $T_d + 1, \ldots, T_d + \tau$.

## 4. IC Likelihood with Marginalization (Option B)

### Core idea

Instead of placing the event at $T_d$, we sum over all possible true transition
weeks $s \in \{T_c + 1, \ldots, T_d\}$, weighting each by its probability and
the downstream $1 \to 2$ likelihood evaluated at the corrected sojourn.

### State 1 (progressed, alive) with IC

Let $T_{\text{cens}}$ be the censoring time for $1 \to 2$.

$$
\boxed{
\log L_{\text{IC}}^{(1)} =
  \underbrace{\sum_{t=1}^{T_c}\bigl[\ell_{01}(t) + \ell_{02}(t) + \ell_{03}(t)\bigr]}_{\text{survived all risks to last clean scan}}
+ \log \sum_{s=T_c+1}^{T_d} \exp\!\bigl(\alpha_s + \beta_s + \gamma_s\bigr)
}
$$

where the three components for each candidate week $s$ are:

$$
\alpha_s = \underbrace{\sum_{t=T_c+1}^{s-1} \ell_{01}(t)}_{\text{survived } 0 \to 1 \text{ from } T_c{+}1 \text{ to } s{-}1}
          + \underbrace{\log\!\bigl(1 - e^{\ell_{01}(s)}\bigr)}_{\text{event at week } s}
$$

$$
\beta_s = \underbrace{\sum_{t=T_c+1}^{s-1} \bigl[\ell_{02}(t) + \ell_{03}(t)\bigr]}_{\text{survived competing risks in state 0}}
$$

$$
\gamma_s = \underbrace{\mathcal{L}_{12}\!\bigl(\tau_s, \text{censored}\bigr)}_{\text{sojourn censored at } \tau_s = T_{\text{cens}} - s}
$$

The $\log\sum\exp(\cdot)$ is computed via `log_sum_exp` for numerical stability.

### State 2 via $0 \to 1 \to 2$ (progressed then died) with IC

$$
\boxed{
\log L_{\text{IC}}^{(2 \text{ via } 1)} =
  \sum_{t=1}^{T_c}\bigl[\ell_{01}(t) + \ell_{02}(t) + \ell_{03}(t)\bigr]
+ \log \sum_{s=T_c+1}^{T_d} \exp\!\bigl(\alpha_s + \beta_s + \gamma_s\bigr)
}
$$

Same $\alpha_s$ and $\beta_s$ as above, but:

$$
\gamma_s = \mathcal{L}_{12}\!\bigl(\tau_s, \text{event}\bigr)
\quad\text{where } \tau_s = T_\dagger - s
$$

### Sojourn $\mathcal{L}_{12}$ with variable start

For semi-Markov, the sojourn likelihood at each candidate $s$ evaluates a
different sojourn duration:

$$
\tau_s = T_{\text{end}} - s
$$

where $T_{\text{end}}$ is $T_\dagger$ (death) or $T_{\text{cens}}$ (censoring).
As $s$ increases by 1, $\tau_s$ decreases by 1, so we peel off one sojourn
interval per iteration.

For Markov (`ms_time_scale_12 = 0`), the clock-forward $1 \to 2$ risk period is
$[s + 1, T_{\text{end}}]$. As $s$ increases, one more week at the start
of the risk window is removed.

For extended (`ms_time_scale_12 = 2`), both sojourn and clock-forward components
shift with $s$.

### When to skip IC

IC is **not applied** when:
1. `ms_prog_deterministic[i] = 1` — progression time pinned by the SLD model
2. `ms_censored_01[i] = 1` — patient is right-censored (no progression detected)
3. `ms_final_state[i] == 2` with `time_01[i] == 0` — death without progression ($0 \to 2$)
4. `ms_final_state[i] == 3` — dropout
5. `ic_gap[i] == 0` — detection at next-week (no interval to marginalize over)

## 5. Stan Pseudocode

### New data fields required

```stan
// In modules/multistate/data.stan:
array[n_patients] int<lower=0> ms_ic_gap_01;
// = T_d - T_c for event patients (the interval censoring gap)
// = 0 for censored patients, deterministic PD, or death-without-PD
```

### Marginalized likelihood for states 1 and 2-via-1

```stan
// Replace the current state-1 and state-2-via-1 blocks:

// --- Shared preamble for IC patients (state 1 or state 2 via 0→1→2) ---
int T_d = time_01[i];               // detection week
int T_c = T_d - ms_ic_gap_01[i];    // last clean assessment
int gap = ms_ic_gap_01[i];          // number of candidate weeks

// Base survival: survived all state-0 risks through T_c
real base_surv = 0;
if (enable_01 && T_c > 0)
  base_surv += sum(log_cond_surv_01[i, 1:T_c]);
if (enable_02 && T_c > 0)
  base_surv += sum(log_cond_surv_02[i, 1:T_c]);
if (enable_03 && T_c > 0)
  base_surv += sum(log_cond_surv_03[i, 1:T_c]);

patient_ll += base_surv;

// --- Marginalize over true transition week s ∈ {T_c+1, ..., T_d} ---

if (gap == 0 || prog_deterministic[i]) {
  // No IC: event pinned at T_d (deterministic PD or gap=0)
  // ... (existing code, but applied only to the T_c+1..T_d interval)
  // 0→1 event hazard at T_d
  if (enable_01 && !prog_deterministic[i])
    patient_ll += log1m_exp(log_cond_surv_01[i, T_d]);
  // 0→2 competing risk survival T_c+1..T_d-1
  if (enable_02 && T_d > T_c + 1)
    patient_ll += sum(log_cond_surv_02[i, (T_c+1):(T_d-1)]);
  if (enable_03 && T_d > T_c + 1)
    patient_ll += sum(log_cond_surv_03[i, (T_c+1):(T_d-1)]);
  // 0→1 survival T_c+1..T_d-1
  if (enable_01 && T_d > T_c + 1)
    patient_ll += sum(log_cond_surv_01[i, (T_c+1):(T_d-1)]);
  // 1→2 sojourn at original time_12
  patient_ll += sojourn_loglik_12(..., time_12[i], ...);

} else {
  // IC marginalization: loop over candidate weeks
  vector[gap] terms;

  // Precompute cumulative sums for efficiency
  // cum_01[k] = Σ_{t=T_c+1}^{T_c+k} ℓ_01(t),  k = 1..gap
  // cum_02[k] = Σ_{t=T_c+1}^{T_c+k} ℓ_02(t)
  // cum_03[k] = Σ_{t=T_c+1}^{T_c+k} ℓ_03(t)
  vector[gap] cum_01 = cumulative_sum(
    to_vector(log_cond_surv_01[i, (T_c+1):T_d]));
  vector[gap] cum_02 = cumulative_sum(
    to_vector(log_cond_surv_02[i, (T_c+1):T_d]));
  vector[gap] cum_03 = cumulative_sum(
    to_vector(log_cond_surv_03[i, (T_c+1):T_d]));

  for (k in 1:gap) {
    int s = T_c + k;  // candidate transition week

    // α: 0→1 survived T_c+1..s-1, event at s
    real alpha;
    if (k == 1) {
      alpha = log1m_exp(log_cond_surv_01[i, s]);
    } else {
      alpha = cum_01[k-1]   // survived T_c+1 .. s-1
            + log1m_exp(log_cond_surv_01[i, s]);  // event at s
    }

    // β: competing risk survival T_c+1..s-1
    real beta = 0;
    if (enable_02 && k > 1)
      beta += cum_02[k-1];
    if (enable_03 && k > 1)
      beta += cum_03[k-1];

    // γ: 1→2 sojourn likelihood with adjusted sojourn
    int tau_s;
    if (final_state[i] == 1) {
      // Censored in state 1: sojourn = T_cens - s
      // T_cens = T_d + time_12[i] (original sojourn was from T_d)
      tau_s = time_12[i] + (T_d - s);  // = time_12[i] + gap - k
    } else {
      // Died: sojourn = T_death - s
      // T_death = T_d + time_12[i] (original sojourn was from T_d)
      tau_s = time_12[i] + (T_d - s);  // same arithmetic
    }

    real gamma = 0;
    if (enable_12 && tau_s > 0) {
      if (final_state[i] == 1) {
        // Censored in state 1
        gamma = sojourn_loglik_12_censored(i, tau_s, s, ...);
      } else {
        // Died in state 2 via state 1
        gamma = sojourn_loglik_12_event(i, tau_s, s, ...);
      }
    }

    terms[k] = alpha + beta + gamma;
  }

  patient_ll += log_sum_exp(terms);
}
```

### Sojourn helper (semi-Markov case)

```stan
// Censored sojourn: survived intervals 1..tau
real sojourn_loglik_12_censored(int i, int tau, int s, ...) {
  return sum(log_cond_surv_12_s[i, 1:tau]);
}

// Event sojourn: survived 1..tau-1, died at tau
real sojourn_loglik_12_event(int i, int tau, int s, ...) {
  real ll = 0;
  if (tau > 1)
    ll += sum(log_cond_surv_12_s[i, 1:(tau - 1)]);
  ll += log1m_exp(log_cond_surv_12_s[i, tau]);
  return ll;
}
```

### Sojourn helper (Markov case, `ms_time_scale_12 = 0`)

```stan
// Censored: survived clock-forward weeks [s+1, s+tau]
real sojourn_loglik_12_censored_markov(int i, int tau, int s, ...) {
  int t_start = s + 1;
  int t_end = s + tau;
  if (t_end <= cols(log_cond_surv_12_t))
    return sum(log_cond_surv_12_t[i, t_start:t_end]);
  return 0;
}

// Event: survived [s+1, s+tau-1], died at s+tau
real sojourn_loglik_12_event_markov(int i, int tau, int s, ...) {
  real ll = 0;
  int t_start = s + 1;
  int t_death = s + tau;
  if (t_death - 1 >= t_start && t_death <= cols(log_cond_surv_12_t))
    ll += sum(log_cond_surv_12_t[i, t_start:(t_death - 1)]);
  if (t_death <= cols(log_cond_surv_12_t))
    ll += log1m_exp(log_cond_surv_12_t[i, t_death]);
  return ll;
}
```

### Sojourn helper (Extended case, `ms_time_scale_12 = 2`)

```stan
// Event: additive hazard from sojourn GP + clock-forward GP
// The combined log-conditional-survival at sojourn step u, clock week s+u is:
//   ℓ_combined(u) = ℓ_12^s(u) + ℓ_12^t(s + u)
real sojourn_loglik_12_event_extended(int i, int tau, int s, ...) {
  real ll = 0;
  // Survived combined hazard for u = 1..tau-1
  for (int u = 1; u < tau; u++) {
    ll += log_cond_surv_12_s[i, u];
    if (s + u <= cols(log_cond_surv_12_t))
      ll += log_cond_surv_12_t[i, s + u];
  }
  // Event at u = tau: must combine BEFORE log1m_exp
  real combined_lcs = log_cond_surv_12_s[i, tau];
  if (s + tau <= cols(log_cond_surv_12_t))
    combined_lcs += log_cond_surv_12_t[i, s + tau];
  ll += log1m_exp(combined_lcs);
  return ll;
}
```

## 6. R Data Preparation Changes

In `r/sclc/prepare_analysis_data.R`, the current code shifts event times to
detection week (lines 189–201). The change:

### Before (current)

```r
mutate(
  # Shift to detection week
  pfs = if_else(!right_censored, pfs + interval_censored + 1L, pfs),
  ms_time_01_full = if_else(!ms_censored_01_full,
    ms_time_01_full + interval_censored + 1L, ms_time_01_full),
  ms_time_12 = case_when(
    ms_time_12 == 0L ~ 0L,
    !ms_censored_01_full ~ pmax(1L, ms_time_12 - interval_censored - 1L),
    TRUE ~ ms_time_12
  ),
)
```

### After (with IC)

```r
mutate(
  # Keep detection week as the event time (unchanged)
  pfs = if_else(!right_censored, pfs + interval_censored + 1L, pfs),
  ms_time_01_full = if_else(!ms_censored_01_full,
    ms_time_01_full + interval_censored + 1L, ms_time_01_full),

  # NEW: IC gap for the multistate model
  # = T_d - T_c = interval_censored + 1 for event patients
  # But only for stochastic (non-RECIST) progressions
  ms_ic_gap_01 = case_when(
    ms_censored_01_full ~ 0L,                  # Censored: no IC
    ms_prog_deterministic == 1L ~ 0L,          # RECIST PD: no IC
    TRUE ~ as.integer(interval_censored + 1L)  # Stochastic PD: IC gap
  ),

  # Sojourn time: keep as detection-to-endpoint
  # The marginalization in Stan adjusts this per candidate week s
  ms_time_12 = case_when(
    ms_time_12 == 0L ~ 0L,
    !ms_censored_01_full ~ pmax(1L, ms_time_12 - interval_censored - 1L),
    TRUE ~ ms_time_12
  ),
)
```

### Stan data list addition

```r
# In the stan_data list:
ms_ic_gap_01 = analysis_data$ms_ic_gap_01,
```

### Stan data declaration

```stan
// In modules/multistate/data.stan:
array[n_patients] int<lower=0> ms_ic_gap_01;
```

## 7. Edge Cases

### `tau_s` could be zero or negative

When $s = T_d$ (the last candidate week equals detection week), $\tau_s$ equals
the original `time_12`. As $s$ decreases toward $T_c + 1$, $\tau_s$ increases.
So $\tau_s$ is always $\geq$ `time_12`, which is already floored at 1. No issue.

Actually: for patients who died in the same week as progression detection
($T_\dagger = T_d$), the original `time_12 = 1` (floored). With earlier $s$,
$\tau_s = 1 + (T_d - s) \geq 2$, which is fine.

### `ms_max_sojourn_t` may need to increase

The maximum sojourn time increases by up to `max(ms_ic_gap_01)` because earlier
progression implies longer sojourn. Ensure `ms_max_sojourn_t` is set to:

```r
ms_max_sojourn_t = max(analysis_data$ms_time_12 + analysis_data$ms_ic_gap_01)
```

### `gap = 1` (progression at next assessment)

When $T_d = T_c + 1$ (detection one week after last clean), the loop has a
single term ($s = T_d$), and the likelihood reduces to the current no-IC form.
The `log_sum_exp` of a single-element vector is just that element.

### Interaction with the SLD model

For `ms_prog_deterministic = 1` patients, the SLD model's RECIST evaluation
determines the progression time. The tumor may have crossed the PD threshold
between visits, but in the joint model the SLD component already handles this
through its latent trajectory. No IC marginalization is needed — the multistate
model receives these as deterministic events.

### Competing risks during the IC interval

During the IC interval $(T_c, T_d]$, the patient is at risk for $0 \to 2$
(death) and $0 \to 3$ (dropout) in addition to $0 \to 1$. The $\beta_s$ term
captures this: for each candidate $s$, the patient survived $0 \to 2$ and
$0 \to 3$ through weeks $T_c + 1, \ldots, s - 1$. At week $s$ itself, the
$0 \to 1$ event fires, which is consistent with the existing convention in the
codebase (competing risk survival runs through `time_01 - 1`).

## 8. Generated Quantities (Prediction)

The generated quantities sampling functions
(`sample_post_progression_death_rng`, etc.) currently receive a single PFS week.
With IC, prediction remains unchanged: when simulating forward, the model
samples a progression time from the hazard and that becomes the exact transition
time. IC is purely an inference-time correction for the observed data likelihood;
it does not affect the generative process.

## 9. Computational Cost

The marginalization loop iterates over `gap` candidate weeks per patient. With a
typical 6-week assessment cycle, `gap` $\leq 7$. The loop body evaluates
cumulative sums (precomputed) and one sojourn likelihood per candidate. Total
additional cost: $O(n_{\text{progressed}} \times \bar{\Delta})$ where
$\bar{\Delta} \approx 5$ — negligible compared to the GP computations.
