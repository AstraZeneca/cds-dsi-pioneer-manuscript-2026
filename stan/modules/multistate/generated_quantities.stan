// ============================================================================
// Multistate Competing Risks CIF (per-trial, empirical subdistribution)
// ============================================================================
// Requires: spop_ms_pfs, spop_ms_right_censored, spop_os, spop_os_censored,
//           sample_ms_pfs, sample_ms_right_censored, sample_os, sample_os_censored
//           (declared and computed in the enclosing GQ block before this include),
//           plus trial_patient_pos, max_all_t, n_trials from transformed data / data.
//
// Classification rule per draw × patient — see compute_trial_cif() in pfs.stanfunctions.
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

    (spop_cif_01[s], spop_cif_02[s], spop_cif_03[s]) = compute_trial_cif(
      spop_ms_pfs[tr_start:tr_end], spop_ms_right_censored[tr_start:tr_end],
      spop_os[tr_start:tr_end],  spop_os_censored[tr_start:tr_end],
      max_all_t
    );

    (sample_cif_01[s], sample_cif_02[s], sample_cif_03[s]) = compute_trial_cif(
      sample_ms_pfs[tr_start:tr_end], sample_ms_right_censored[tr_start:tr_end],
      sample_os[tr_start:tr_end],  sample_os_censored[tr_start:tr_end],
      max_all_t
    );
  }
}
