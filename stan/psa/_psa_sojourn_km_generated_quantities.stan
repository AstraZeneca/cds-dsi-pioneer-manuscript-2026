// ============================================================================
// SOJOURN KM CURVES — PSA FLAT VERSION (no trial dimension)
// ============================================================================
// Flat version of modules/multistate/sojourn_km_generated_quantities.stan.
// Used by pioneer.stan only (ms-standalone still uses the trial-indexed version).
//
// Requires same contract as the trial-indexed version:
//   sample_ms_pfs, sample_ms_right_censored, sample_os, sample_os_censored,
//   sample_pfs, sample_is_dropout  (and spop_ variants)
//   ms_max_sojourn_t, ms_max_sojourn_t_32, estimate_kaplan_meier()

vector<lower=0, upper=1>[ms_max_sojourn_t + 1]
  sample_km_12 = ones_vector(ms_max_sojourn_t + 1),
  spop_km_12   = ones_vector(ms_max_sojourn_t + 1);
vector<lower=0, upper=1>[ms_max_sojourn_t_32 + 1]
  sample_km_32 = ones_vector(ms_max_sojourn_t_32 + 1),
  spop_km_32   = ones_vector(ms_max_sojourn_t_32 + 1);

// ── 1→2: post-progression sojourn KM ─────────────────────────────────────────
{
  int n_prog = n_forecast_patients - sum(sample_ms_right_censored);
  if (n_prog > 0) {
    array[n_prog] int soj; array[n_prog] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (!sample_ms_right_censored[i]) {
        soj[idx]  = max(1, sample_os[i] - sample_ms_pfs[i]);
        cens[idx] = sample_os_censored[i];
        idx += 1;
      }
    }
    sample_km_12 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
  }
}
{
  int n_prog = n_forecast_patients - sum(spop_ms_right_censored);
  if (n_prog > 0) {
    array[n_prog] int soj; array[n_prog] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (!spop_ms_right_censored[i]) {
        soj[idx]  = max(1, spop_os[i] - spop_ms_pfs[i]);
        cens[idx] = spop_os_censored[i];
        idx += 1;
      }
    }
    spop_km_12 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
  }
}

// ── 3→2: off-trial sojourn KM ─────────────────────────────────────────────────
{
  int n_drop = sum(to_array_1d(sample_is_dropout));
  if (n_drop > 0) {
    array[n_drop] int soj; array[n_drop] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (sample_is_dropout[i]) {
        soj[idx]  = max(1, sample_os[i] - sample_pfs[i]);
        cens[idx] = sample_os_censored[i];
        idx += 1;
      }
    }
    sample_km_32 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
  }
}
{
  int n_drop = sum(spop_is_dropout);
  if (n_drop > 0) {
    array[n_drop] int soj; array[n_drop] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (spop_is_dropout[i]) {
        soj[idx]  = max(1, spop_os[i] - spop_pfs[i]);
        cens[idx] = spop_os_censored[i];
        idx += 1;
      }
    }
    spop_km_32 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
  }
}
