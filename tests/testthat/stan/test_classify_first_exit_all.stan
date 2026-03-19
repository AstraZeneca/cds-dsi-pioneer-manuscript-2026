// Tests for classify-first endpoint functions (Issue 79)
//
// Covers all four pure functions in pfs.stanfunctions and the deterministic
// branches of the two OS RNG functions in multistate.stanfunctions.
//
// Cause enum (shared by all functions):
//   0 = CENSORED     – no state-0 exit within prediction horizon
//   1 = PROGRESSION  – 0→1 fired first
//   2 = DIRECT_DEATH – 0→2 fired first
//   3 = DROPOUT      – 0→3 fired first
//
// Run via test_stan_function() with fixed_param=TRUE, chains=1, iter_sampling=1.
// The R test asserts n_failures == 0.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}

data {
  int<lower=0> dummy;  // Required non-empty data block; pass dummy = 0L from R
}

generated quantities {
  int n_failures = 0;

  // ===========================================================================
  // classify_spop_exit(t_01, c_01, t_02, c_02, t_03, c_03,
  //                    enable_ms_02, enable_ms_03)
  // Priority on ties: 0→1 > 0→2 > 0→3
  // When a transition is disabled (enable=0), caller sets c_xx=1.
  // ===========================================================================

  // S1: Progression wins — only 0→1 fires; 02 and 03 censored
  {
    tuple(int, int) r = classify_spop_exit(10, 0, 101, 1, 101, 1, 1, 1);
    if (r.1 != 1 || r.2 != 10) {
      print("FAIL S1 classify_spop_exit progression wins: got (", r.1, ",", r.2, ") expected (1,10)");
      n_failures += 1;
    }
  }

  // S2: Direct death wins — 0→2 beats 0→1 (t_02=5 < t_01=15)
  {
    tuple(int, int) r = classify_spop_exit(15, 0, 5, 0, 101, 1, 1, 1);
    if (r.1 != 2 || r.2 != 5) {
      print("FAIL S2 classify_spop_exit direct death wins: got (", r.1, ",", r.2, ") expected (2,5)");
      n_failures += 1;
    }
  }

  // S3: Dropout wins — 0→3 beats 0→1 (t_03=8 < t_01=20); 02 censored
  {
    tuple(int, int) r = classify_spop_exit(20, 0, 101, 1, 8, 0, 1, 1);
    if (r.1 != 3 || r.2 != 8) {
      print("FAIL S3 classify_spop_exit dropout wins: got (", r.1, ",", r.2, ") expected (3,8)");
      n_failures += 1;
    }
  }

  // S4: All censored — cause=0, exit=max of all three times
  {
    tuple(int, int) r = classify_spop_exit(50, 1, 101, 1, 101, 1, 1, 1);
    if (r.1 != 0 || r.2 != 101) {
      print("FAIL S4 classify_spop_exit all censored: got (", r.1, ",", r.2, ") expected (0,101)");
      n_failures += 1;
    }
  }

  // S5: 0→1 and 0→2 tie at t=10 — progression wins (01 priority)
  {
    tuple(int, int) r = classify_spop_exit(10, 0, 10, 0, 101, 1, 1, 1);
    if (r.1 != 1 || r.2 != 10) {
      print("FAIL S5 classify_spop_exit 01/02 tie progression wins: got (", r.1, ",", r.2, ") expected (1,10)");
      n_failures += 1;
    }
  }

  // S6: 0→2 and 0→3 tie at t=8; 0→1 censored — death wins (02 priority over 03)
  {
    tuple(int, int) r = classify_spop_exit(50, 1, 8, 0, 8, 0, 1, 1);
    if (r.1 != 2 || r.2 != 8) {
      print("FAIL S6 classify_spop_exit 02/03 tie death wins: got (", r.1, ",", r.2, ") expected (2,8)");
      n_failures += 1;
    }
  }

  // S7: 0→1 censored; 0→2 fires — cause=2
  {
    tuple(int, int) r = classify_spop_exit(50, 1, 5, 0, 101, 1, 1, 1);
    if (r.1 != 2 || r.2 != 5) {
      print("FAIL S7 classify_spop_exit 01 censored 02 fires: got (", r.1, ",", r.2, ") expected (2,5)");
      n_failures += 1;
    }
  }

  // S8: All three fire simultaneously at t=10 — progression wins (triple tie)
  {
    tuple(int, int) r = classify_spop_exit(10, 0, 10, 0, 10, 0, 1, 1);
    if (r.1 != 1 || r.2 != 10) {
      print("FAIL S8 classify_spop_exit triple tie progression wins: got (", r.1, ",", r.2, ") expected (1,10)");
      n_failures += 1;
    }
  }

  // S9: enable_ms_02=0 suppresses direct death — even when c_02=0 (event fired),
  //     died_directly=0; with 0→1 also censored, result is cause=0
  {
    tuple(int, int) r = classify_spop_exit(50, 1, 5, 0, 101, 1, 0, 1);
    if (r.1 != 0 || r.2 != 101) {
      print("FAIL S9 classify_spop_exit enable_ms_02=0 suppresses: got (", r.1, ",", r.2, ") expected (0,101)");
      n_failures += 1;
    }
  }

  // ===========================================================================
  // derive_spop_pfs(cause, exit_time, t_target, c_target,
  //                 t_ms01, c_ms01, t_03,
  //                 enable_ms_02, enable_ms_03)
  // Returns (pfs, right_censored, ms_pfs, ms_right_censored)
  // ===========================================================================

  // D1: cause=1 (progression) — pfs and ms_pfs both equal exit_time; uncensored
  {
    tuple(int, int, int, int) r = derive_spop_pfs(1, 10, 10, 0, 12, 0, 50, 1, 1);
    if (r.1 != 10 || r.2 != 0 || r.3 != 10 || r.4 != 0) {
      print("FAIL D1 derive_spop_pfs cause=1: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (10,0,10,0)");
      n_failures += 1;
    }
  }

  // D2: cause=2 (direct death) — pfs and ms_pfs both equal exit_time; uncensored
  {
    tuple(int, int, int, int) r = derive_spop_pfs(2, 5, 15, 0, 20, 0, 50, 1, 1);
    if (r.1 != 5 || r.2 != 0 || r.3 != 5 || r.4 != 0) {
      print("FAIL D2 derive_spop_pfs cause=2: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (5,0,5,0)");
      n_failures += 1;
    }
  }

  // D3: cause=3, SLD progression before dropout (t_target=25 < t_03=30)
  //     — PFS event at t_target; ms_pfs censored at t_03
  {
    tuple(int, int, int, int) r = derive_spop_pfs(3, 30, 25, 0, 40, 1, 30, 1, 1);
    if (r.1 != 25 || r.2 != 0 || r.3 != 30 || r.4 != 1) {
      print("FAIL D3 derive_spop_pfs cause=3 SLD before dropout: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (25,0,30,1)");
      n_failures += 1;
    }
  }

  // D4: cause=3, SLD progression AT dropout boundary (t_target=30 == t_03=30)
  //     — still counts as a PFS event (<=, not <)
  {
    tuple(int, int, int, int) r = derive_spop_pfs(3, 30, 30, 0, 40, 1, 30, 1, 1);
    if (r.1 != 30 || r.2 != 0 || r.3 != 30 || r.4 != 1) {
      print("FAIL D4 derive_spop_pfs cause=3 SLD at boundary: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (30,0,30,1)");
      n_failures += 1;
    }
  }

  // D5: cause=3, SLD progression AFTER dropout (t_target=40 > t_03=30)
  //     — censored at t_03
  {
    tuple(int, int, int, int) r = derive_spop_pfs(3, 30, 40, 0, 40, 1, 30, 1, 1);
    if (r.1 != 30 || r.2 != 1 || r.3 != 30 || r.4 != 1) {
      print("FAIL D5 derive_spop_pfs cause=3 SLD after dropout: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (30,1,30,1)");
      n_failures += 1;
    }
  }

  // D6: cause=3, SLD censored (c_target=1) — censored at t_03
  {
    tuple(int, int, int, int) r = derive_spop_pfs(3, 30, 25, 1, 40, 1, 30, 1, 1);
    if (r.1 != 30 || r.2 != 1 || r.3 != 30 || r.4 != 1) {
      print("FAIL D6 derive_spop_pfs cause=3 SLD censored: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (30,1,30,1)");
      n_failures += 1;
    }
  }

  // D7: cause=0 (all censored) — pfs=min(t_target, t_ms01); ms_pfs=t_ms01 raw
  {
    tuple(int, int, int, int) r = derive_spop_pfs(0, 50, 20, 1, 25, 1, 50, 1, 1);
    if (r.1 != 20 || r.2 != 1 || r.3 != 25 || r.4 != 1) {
      print("FAIL D7 derive_spop_pfs cause=0: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (20,1,25,1)");
      n_failures += 1;
    }
  }

  // ===========================================================================
  // classify_sample_exit(sample_time_01, sample_cens_01,
  //                      sample_time_02, sample_cens_02,
  //                      ms_final_state_i, ms_censored_02_i, ms_time_03_i,
  //                      enable_ms_02, enable_ms_03)
  // Uses observed ground truth: state=3 overrides all forecasts
  // ===========================================================================

  // C1: Observed dropout (ms_final_state=3) — overrides direct death and progression
  {
    tuple(int, int) r = classify_sample_exit(15, 0, 5, 0, 3, 0, 20, 1, 1);
    if (r.1 != 3 || r.2 != 20) {
      print("FAIL C1 classify_sample_exit dropout overrides: got (", r.1, ",", r.2, ") expected (3,20)");
      n_failures += 1;
    }
  }

  // C2: Observed direct death (ms_censored_02=0, ms_final_state=2) beats progression
  //     (t_02=8 < t_01=15)
  {
    tuple(int, int) r = classify_sample_exit(15, 0, 8, 0, 2, 0, 50, 1, 1);
    if (r.1 != 2 || r.2 != 8) {
      print("FAIL C2 classify_sample_exit observed death: got (", r.1, ",", r.2, ") expected (2,8)");
      n_failures += 1;
    }
  }

  // C3: Forecasted direct death (ms_censored_02=1 but sample_cens_02=0 and
  //     sample_time_02=5 < sample_time_01=15) — died_directly fires
  {
    tuple(int, int) r = classify_sample_exit(15, 0, 5, 0, 0, 1, 50, 1, 1);
    if (r.1 != 2 || r.2 != 5) {
      print("FAIL C3 classify_sample_exit forecasted death beats progression: got (", r.1, ",", r.2, ") expected (2,5)");
      n_failures += 1;
    }
  }

  // C4: Progression first — 0→1 fires before 0→2 (t_01=10 < t_02=20)
  {
    tuple(int, int) r = classify_sample_exit(10, 0, 20, 0, 1, 1, 50, 1, 1);
    if (r.1 != 1 || r.2 != 10) {
      print("FAIL C4 classify_sample_exit progression first: got (", r.1, ",", r.2, ") expected (1,10)");
      n_failures += 1;
    }
  }

  // C5: All censored — cause=0, exit=max(sample_time_01, sample_time_02)
  {
    tuple(int, int) r = classify_sample_exit(20, 1, 25, 1, 0, 1, 50, 1, 1);
    if (r.1 != 0 || r.2 != 25) {
      print("FAIL C5 classify_sample_exit all censored: got (", r.1, ",", r.2, ") expected (0,25)");
      n_failures += 1;
    }
  }

  // C6: Progression wins over later direct death (t_01=10 <= t_02=20)
  //     — no dropout, 02 fires after 01
  {
    tuple(int, int) r = classify_sample_exit(10, 0, 20, 0, 0, 1, 50, 1, 1);
    if (r.1 != 1 || r.2 != 10) {
      print("FAIL C6 classify_sample_exit progression over later death: got (", r.1, ",", r.2, ") expected (1,10)");
      n_failures += 1;
    }
  }

  // C7: enable_ms_02=0 suppresses direct death — progression wins despite simultaneous
  //     0→2 event (sample_cens_02=0), because died_directly is gated by enable_ms_02
  {
    tuple(int, int) r = classify_sample_exit(10, 0, 10, 0, 0, 1, 50, 0, 1);
    if (r.1 != 1 || r.2 != 10) {
      print("FAIL C7 classify_sample_exit enable_ms_02=0 suppresses: got (", r.1, ",", r.2, ") expected (1,10)");
      n_failures += 1;
    }
  }

  // ===========================================================================
  // derive_sample_pfs(cause, exit_time,
  //                   sample_target_pfs_i, sample_target_right_cens_i,
  //                   sample_ms_pfs_raw, sample_ms_right_cens_raw,
  //                   ms_time_03_i, enable_ms_02, enable_ms_03)
  // Returns (pfs, right_censored, ms_pfs, ms_right_censored)
  // ===========================================================================

  // P1: cause=2 (direct death) — both pfs and ms_pfs set to exit; uncensored
  {
    tuple(int, int, int, int) r = derive_sample_pfs(2, 10, 20, 0, 15, 0, 50, 1, 1);
    if (r.1 != 10 || r.2 != 0 || r.3 != 10 || r.4 != 0) {
      print("FAIL P1 derive_sample_pfs cause=2: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (10,0,10,0)");
      n_failures += 1;
    }
  }

  // P2: cause=3, SLD PD before dropout (sample_target_pfs=25 <= ms_time_03=30)
  //     — PFS event at sample_target_pfs; ms_pfs censored at ms_time_03
  {
    tuple(int, int, int, int) r = derive_sample_pfs(3, 30, 25, 0, 30, 1, 30, 1, 1);
    if (r.1 != 25 || r.2 != 0 || r.3 != 30 || r.4 != 1) {
      print("FAIL P2 derive_sample_pfs cause=3 SLD before dropout: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (25,0,30,1)");
      n_failures += 1;
    }
  }

  // P3: cause=3, SLD PD at dropout boundary (==) — still a PFS event
  {
    tuple(int, int, int, int) r = derive_sample_pfs(3, 30, 30, 0, 30, 1, 30, 1, 1);
    if (r.1 != 30 || r.2 != 0 || r.3 != 30 || r.4 != 1) {
      print("FAIL P3 derive_sample_pfs cause=3 SLD at boundary: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (30,0,30,1)");
      n_failures += 1;
    }
  }

  // P4: cause=3, SLD PD after dropout — censored at ms_time_03
  {
    tuple(int, int, int, int) r = derive_sample_pfs(3, 30, 35, 0, 30, 1, 30, 1, 1);
    if (r.1 != 30 || r.2 != 1 || r.3 != 30 || r.4 != 1) {
      print("FAIL P4 derive_sample_pfs cause=3 SLD after dropout: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (30,1,30,1)");
      n_failures += 1;
    }
  }

  // P5: cause=3, SLD censored (sample_target_right_cens=1) — censored at ms_time_03
  {
    tuple(int, int, int, int) r = derive_sample_pfs(3, 30, 25, 1, 30, 1, 30, 1, 1);
    if (r.1 != 30 || r.2 != 1 || r.3 != 30 || r.4 != 1) {
      print("FAIL P5 derive_sample_pfs cause=3 SLD censored: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (30,1,30,1)");
      n_failures += 1;
    }
  }

  // P6: cause=1 (progression) — pfs=exit_time (event), ms_pfs stays at raw sample
  {
    tuple(int, int, int, int) r = derive_sample_pfs(1, 10, 10, 0, 10, 0, 50, 1, 1);
    if (r.1 != 10 || r.2 != 0 || r.3 != 10 || r.4 != 0) {
      print("FAIL P6 derive_sample_pfs cause=1: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (10,0,10,0)");
      n_failures += 1;
    }
  }

  // P7: cause=0 (censored) — pfs=exit_time (censored), ms_pfs stays at raw sample
  {
    tuple(int, int, int, int) r = derive_sample_pfs(0, 20, 20, 1, 25, 1, 50, 1, 1);
    if (r.1 != 20 || r.2 != 1 || r.3 != 25 || r.4 != 1) {
      print("FAIL P7 derive_sample_pfs cause=0: got (", r.1, ",", r.2, ",", r.3, ",", r.4, ") expected (20,1,25,1)");
      n_failures += 1;
    }
  }

  // ===========================================================================
  // derive_spop_os_rng — deterministic branches only (no actual RNG calls)
  //
  // The row_vector arguments are never read in the branches tested here, but the
  // compiler requires them to be valid row_vectors. We pass rep_row_vector(0.0, 1).
  //
  // Branches NOT tested here (they call sample_*_rng):
  //   cause=3 && enable_ms_32=1
  //   cause=1 && enable_ms_12=1
  // ===========================================================================

  {
    row_vector[1] dummy_rv = rep_row_vector(0.0, 1);

    // O1: cause=2 (direct death) — returns (exit_time, 0) immediately
    {
      tuple(int, int) r = derive_spop_os_rng(2, 10, 1, 1, 0,
          dummy_rv, dummy_rv, dummy_rv, 10, 25, 30);
      if (r.1 != 10 || r.2 != 0) {
        print("FAIL O1 derive_spop_os_rng cause=2: got (", r.1, ",", r.2, ") expected (10,0)");
        n_failures += 1;
      }
    }

    // O2: cause=3, enable_ms_32=0 — censor OS at dropout time
    {
      tuple(int, int) r = derive_spop_os_rng(3, 20, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv, 10, 25, 20);
      if (r.1 != 20 || r.2 != 1) {
        print("FAIL O2 derive_spop_os_rng cause=3 no 3->2: got (", r.1, ",", r.2, ") expected (20,1)");
        n_failures += 1;
      }
    }

    // O3: cause=1, enable_ms_12=0 — censor OS at max(t_01, t_02, t_03)
    {
      tuple(int, int) r = derive_spop_os_rng(1, 10, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv, 10, 25, 30);
      if (r.1 != 30 || r.2 != 1) {
        print("FAIL O3 derive_spop_os_rng cause=1 no 1->2: got (", r.1, ",", r.2, ") expected (30,1)");
        n_failures += 1;
      }
    }

    // O4: cause=0 (all censored) — censor OS at max(t_01, t_02, t_03)
    {
      tuple(int, int) r = derive_spop_os_rng(0, 50, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv, 10, 25, 30);
      if (r.1 != 30 || r.2 != 1) {
        print("FAIL O4 derive_spop_os_rng cause=0: got (", r.1, ",", r.2, ") expected (30,1)");
        n_failures += 1;
      }
    }
  }

  // ===========================================================================
  // derive_sample_os_rng — deterministic branches only (no actual RNG calls)
  //
  // Branches NOT tested here (they call sample_*_rng):
  //   cause=3 && ms_censored_32=1 && enable_ms_32=1
  //   cause=1 && ms_censored_12=1 && enable_ms_12=1
  // ===========================================================================

  {
    row_vector[1] dummy_rv = rep_row_vector(0.0, 1);

    // Q1: cause=2 (direct death) — returns (exit_time, 0) immediately
    {
      tuple(int, int) r = derive_sample_os_rng(2, 10, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv,
          10, 25, 10, 1, 0, 1, 0, 0);
      if (r.1 != 10 || r.2 != 0) {
        print("FAIL Q1 derive_sample_os_rng cause=2: got (", r.1, ",", r.2, ") expected (10,0)");
        n_failures += 1;
      }
    }

    // Q2: cause=3, observed off-trial death (!ms_censored_32=0)
    //     — returns (ms_time_03 + ms_time_32, 0) = (20+15, 0) = (35, 0)
    {
      tuple(int, int) r = derive_sample_os_rng(3, 20, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv,
          20, 50, 20, 1, 0, 0, 20, 15);
      if (r.1 != 35 || r.2 != 0) {
        print("FAIL Q2 derive_sample_os_rng cause=3 observed death: got (", r.1, ",", r.2, ") expected (35,0)");
        n_failures += 1;
      }
    }

    // Q3: cause=3, ms_censored_32=1, enable_ms_32=0
    //     — censor at (ms_time_03 + ms_time_32, 1) = (20+10, 1) = (30, 1)
    {
      tuple(int, int) r = derive_sample_os_rng(3, 20, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv,
          20, 50, 20, 1, 0, 1, 20, 10);
      if (r.1 != 30 || r.2 != 1) {
        print("FAIL Q3 derive_sample_os_rng cause=3 censored no 3->2: got (", r.1, ",", r.2, ") expected (30,1)");
        n_failures += 1;
      }
    }

    // Q4: cause=1, observed post-progression death (!ms_censored_12=0)
    //     — returns (pfs_i + ms_time_12, 0) = (10+8, 0) = (18, 0)
    {
      tuple(int, int) r = derive_sample_os_rng(1, 10, 1, 0, 0,
          dummy_rv, dummy_rv, dummy_rv,
          10, 50, 10, 0, 8, 1, 0, 0);
      if (r.1 != 18 || r.2 != 0) {
        print("FAIL Q4 derive_sample_os_rng cause=1 observed PPD: got (", r.1, ",", r.2, ") expected (18,0)");
        n_failures += 1;
      }
    }

    // Q5: cause=1, enable_ms_12=0, ms_censored_12=1
    //     — censor at max(sample_time_01, sample_time_02) = max(10,15) = 15
    {
      tuple(int, int) r = derive_sample_os_rng(1, 10, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv,
          10, 15, 10, 1, 0, 1, 0, 0);
      if (r.1 != 15 || r.2 != 1) {
        print("FAIL Q5 derive_sample_os_rng cause=1 no 1->2: got (", r.1, ",", r.2, ") expected (15,1)");
        n_failures += 1;
      }
    }

    // Q6: cause=0 (all censored) — censor at max(sample_time_01, sample_time_02) = 20
    {
      tuple(int, int) r = derive_sample_os_rng(0, 20, 0, 0, 0,
          dummy_rv, dummy_rv, dummy_rv,
          10, 20, 10, 1, 0, 1, 0, 0);
      if (r.1 != 20 || r.2 != 1) {
        print("FAIL Q6 derive_sample_os_rng cause=0: got (", r.1, ",", r.2, ") expected (20,1)");
        n_failures += 1;
      }
    }
  }
}
