// ============================================================================
// SOJOURN KM CURVES (1→2 post-progression, 3→2 off-trial)
// ============================================================================
// Shared include for both standalone multistate and full PSA models.
//
// Requires (declared and computed before this include):
//   - sample_ms_pfs, sample_ms_right_censored, sample_os, sample_os_censored,
//     sample_pfs, sample_is_dropout  (and spop_ variants)
//   - forecast_trial_patient_pos, n_trials, ms_max_sojourn_t, ms_max_sojourn_t_32
//   - estimate_kaplan_meier() function
//
// sample: conditional on observed progression/dropout (patients who progressed/
//         dropped out in the observed data contribute observed sojourn times;
//         those still in state 0/1 at follow-up contribute censored sojourn times
//         drawn from the posterior).
// spop:   unconditional — every patient contributes a simulated sojourn time
//         drawn from the posterior (whether or not they actually progressed/dropped
//         out in the data).

array[n_trials] vector<lower=0, upper=1>[ms_max_sojourn_t + 1]
  sample_km_12 = rep_array(ones_vector(ms_max_sojourn_t + 1), n_trials),
  spop_km_12   = rep_array(ones_vector(ms_max_sojourn_t + 1), n_trials);
array[n_trials] vector<lower=0, upper=1>[ms_max_sojourn_t_32 + 1]
  sample_km_32 = rep_array(ones_vector(ms_max_sojourn_t_32 + 1), n_trials),
  spop_km_32   = rep_array(ones_vector(ms_max_sojourn_t_32 + 1), n_trials);

for (s in 1:n_trials) {
  int n_tr = get_pos_size(forecast_trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(forecast_trial_patient_pos, s);

    // ── 1→2: post-progression sojourn KM ─────────────────────────────────────
    // sample: patients who progressed (sample_ms_right_censored == 0)
    {
      int n_prog = n_tr - sum(sample_ms_right_censored[tr_start:tr_end]);
      if (n_prog > 0) {
        array[n_prog] int soj; array[n_prog] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (!sample_ms_right_censored[i]) {
            soj[idx]  = max(1, sample_os[i] - sample_ms_pfs[i]);
            cens[idx] = sample_os_censored[i];
            idx += 1;
          }
        }
        sample_km_12[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
      }
    }

    // spop: all patients get a simulated post-progression sojourn
    {
      int n_prog = n_tr - sum(spop_ms_right_censored[tr_start:tr_end]);
      if (n_prog > 0) {
        array[n_prog] int soj; array[n_prog] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (!spop_ms_right_censored[i]) {
            soj[idx]  = max(1, spop_os[i] - spop_ms_pfs[i]);
            cens[idx] = spop_os_censored[i];
            idx += 1;
          }
        }
        spop_km_12[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
      }
    }

    // ── 3→2: off-trial sojourn KM ─────────────────────────────────────────────
    // sample: patients who dropped out (sample_is_dropout == 1)
    {
      int n_drop = sum(to_array_1d(sample_is_dropout[tr_start:tr_end]));
      if (n_drop > 0) {
        array[n_drop] int soj; array[n_drop] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (sample_is_dropout[i]) {
            soj[idx]  = max(1, sample_os[i] - sample_pfs[i]);
            cens[idx] = sample_os_censored[i];
            idx += 1;
          }
        }
        sample_km_32[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
      }
    }

    // spop: all patients get a simulated off-trial sojourn
    {
      int n_drop = sum(spop_is_dropout[tr_start:tr_end]);
      if (n_drop > 0) {
        array[n_drop] int soj; array[n_drop] int cens;
        int idx = 1;
        for (i in tr_start:tr_end) {
          if (spop_is_dropout[i]) {
            soj[idx]  = max(1, spop_os[i] - spop_pfs[i]);
            cens[idx] = spop_os_censored[i];
            idx += 1;
          }
        }
        spop_km_32[s] = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
      }
    }
  }
}
