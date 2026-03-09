// ============================================================================
// Multistate Competing Risks CIF (per-trial, empirical subdistribution)
// ============================================================================
// Requires: spop_ms_pfs, spop_ms_right_censored, spop_os, spop_os_censored,
//           sample_ms_pfs, sample_ms_right_censored, sample_os, sample_os_censored
//           (declared and computed in the enclosing GQ block before this include),
//           plus trial_patient_pos, max_all_t, n_trials from transformed data / data.
//
// Classification rule per draw × patient (from state 0):
//   0→2 direct death : PFS event AND OS death AND pfs_time == os_time
//   0→1 progression  : PFS event AND NOT 0→2
//   0→3 dropout      : PFS censored AND pfs_time <= max_all_t  (within window)
//   fully censored   : PFS censored AND pfs_time > max_all_t   (no state-0 exit)
//
// CIF[s][t] = fraction of trial-s patients with cause j AND exit by week t-1
// (1-indexed vector: index 1 = time 0, matching the KM vector convention)

// ── Declarations ─────────────────────────────────────────────────────────────

array[n_trials] vector<lower=0, upper=1>[max_all_t + 1]
  spop_cif_01   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  spop_cif_02   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  spop_cif_03   = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_01 = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_02 = rep_array(zeros_vector(max_all_t + 1), n_trials),
  sample_cif_03 = rep_array(zeros_vector(max_all_t + 1), n_trials);

// ── Computation ──────────────────────────────────────────────────────────────

for (s in 1:n_trials) {
  int n_tr = get_pos_size(trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(trial_patient_pos, s);

    // Event counts at each time-index position (O(n_patients + max_t) total)
    array[max_all_t + 1] int cnt_s01  = zeros_int_array(max_all_t + 1);
    array[max_all_t + 1] int cnt_s02  = zeros_int_array(max_all_t + 1);
    array[max_all_t + 1] int cnt_s03  = zeros_int_array(max_all_t + 1);
    array[max_all_t + 1] int cnt_sm01 = zeros_int_array(max_all_t + 1);
    array[max_all_t + 1] int cnt_sm02 = zeros_int_array(max_all_t + 1);
    array[max_all_t + 1] int cnt_sm03 = zeros_int_array(max_all_t + 1);

    for (j in tr_start:tr_end) {

      // ── spop (unconditional posterior predictive) ─────────────────────────
      if (spop_ms_right_censored[j] == 0) {
        // PFS event: distinguish direct death (0→2) from progression (0→1)
        if (spop_os_censored[j] == 0 && spop_ms_pfs[j] == spop_os[j]) {
          cnt_s02[spop_ms_pfs[j]] += 1;   // 0→2: died without progressing
        } else {
          cnt_s01[spop_ms_pfs[j]] += 1;   // 0→1: progressed
        }
      } else if (spop_ms_pfs[j] <= max_all_t) {
        cnt_s03[spop_ms_pfs[j]] += 1;     // 0→3: dropped out within window
      }
      // else: fully censored (pfs_time > max_all_t) — no CIF contribution

      // ── sample (conditional on observed data) ─────────────────────────────
      if (sample_ms_right_censored[j] == 0) {
        if (sample_os_censored[j] == 0 && sample_ms_pfs[j] == sample_os[j]) {
          cnt_sm02[sample_ms_pfs[j]] += 1;
        } else {
          cnt_sm01[sample_ms_pfs[j]] += 1;
        }
      } else if (sample_ms_pfs[j] <= max_all_t) {
        cnt_sm03[sample_ms_pfs[j]] += 1;
      }
    }

    // Cumulative sum → normalized CIF
    int c01 = 0; int c02 = 0; int c03 = 0;
    int cm01 = 0; int cm02 = 0; int cm03 = 0;
    for (t in 1:max_all_t + 1) {
      c01  += cnt_s01[t];   c02  += cnt_s02[t];   c03  += cnt_s03[t];
      cm01 += cnt_sm01[t];  cm02 += cnt_sm02[t];  cm03 += cnt_sm03[t];
      spop_cif_01[s][t]   = c01  * 1.0 / n_tr;
      spop_cif_02[s][t]   = c02  * 1.0 / n_tr;
      spop_cif_03[s][t]   = c03  * 1.0 / n_tr;
      sample_cif_01[s][t] = cm01 * 1.0 / n_tr;
      sample_cif_02[s][t] = cm02 * 1.0 / n_tr;
      sample_cif_03[s][t] = cm03 * 1.0 / n_tr;
    }
  }
}
