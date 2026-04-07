// Tests for compute_trial_cif() (pfs.stanfunctions)
//
// Classification rule:
//   0→2 direct death : rc=0 AND os_censored=0 AND pfs==os
//   0→1 progression  : rc=0 AND NOT 0→2
//   0→3 dropout      : rc=1 AND is_dropout==1 (explicit cause flag, NOT inferred from pfs<=max_t)
//   admin censored   : rc=1 AND is_dropout==0  (no CIF contribution, even if pfs<=max_t)
//
// Four scenarios, 4 patients each (n=4 → fractions in quarters):
//   CIF-A: one of each cause type (0→1, 0→2, 0→3, fully-censored)
//   CIF-B: two 0→1 events at different times (tests cumulative monotonicity)
//   CIF-C: all fully censored (all CIF vectors stay zero)
//   CIF-D: admin-censored patient within horizon is NOT counted as dropout
//
// Output: n_failures == 0. R test asserts n_failures == 0.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}

data {
  int<lower=0> dummy;
}

generated quantities {
  int n_failures = 0;
  real tol = 1e-10;

  // ===========================================================================
  // CIF-A: 4 patients, one of each cause type
  //
  //   j=1: 0→1 at index 3  (pfs=3, rc=0, os=7, os_cens=0)
  //   j=2: 0→2 at index 2  (pfs=2, rc=0, os=2, os_cens=0 — pfs==os)
  //   j=3: 0→3 at index 5  (pfs=5, rc=1, is_dropout=1)
  //   j=4: fully censored  (pfs=9, rc=1, is_dropout=0 — pfs > max_t=8)
  //
  //   cnt_01[3]=1, cnt_02[2]=1, cnt_03[5]=1
  //   cif_01[t]: 0 for t<3, 1/4 for t>=3
  //   cif_02[t]: 0 for t<2, 1/4 for t>=2
  //   cif_03[t]: 0 for t<5, 1/4 for t>=5
  //   total at max (index 9) = 3/4  (1 patient fully censored)
  // ===========================================================================
  {
    int max_t = 8;
    array[4] int pfs            = {3, 2, 5, 9};
    array[4] int right_censored = {0, 0, 1, 1};
    array[4] int is_dropout     = {0, 0, 1, 0};
    array[4] int os             = {7, 2, 5, 9};
    array[4] int os_censored    = {0, 0, 1, 1};

    vector[max_t + 1] c01; vector[max_t + 1] c02; vector[max_t + 1] c03;
    (c01, c02, c03) = compute_trial_cif(pfs, right_censored, is_dropout, os, os_censored, max_t);

    real q = 0.25;  // 1/4

    // cif_01: zero before index 3, then 1/4
    if (abs(c01[2] - 0.0) > tol) { print("FAIL A cif_01[2]=", c01[2]); n_failures += 1; }
    if (abs(c01[3] - q)   > tol) { print("FAIL A cif_01[3]=", c01[3]); n_failures += 1; }
    if (abs(c01[9] - q)   > tol) { print("FAIL A cif_01[9]=", c01[9]); n_failures += 1; }

    // cif_02: zero before index 2, then 1/4
    if (abs(c02[1] - 0.0) > tol) { print("FAIL A cif_02[1]=", c02[1]); n_failures += 1; }
    if (abs(c02[2] - q)   > tol) { print("FAIL A cif_02[2]=", c02[2]); n_failures += 1; }
    if (abs(c02[9] - q)   > tol) { print("FAIL A cif_02[9]=", c02[9]); n_failures += 1; }

    // cif_03: zero before index 5, then 1/4
    if (abs(c03[4] - 0.0) > tol) { print("FAIL A cif_03[4]=", c03[4]); n_failures += 1; }
    if (abs(c03[5] - q)   > tol) { print("FAIL A cif_03[5]=", c03[5]); n_failures += 1; }
    if (abs(c03[9] - q)   > tol) { print("FAIL A cif_03[9]=", c03[9]); n_failures += 1; }

    // total at final index = 3/4 (1 patient fully censored, no contribution)
    if (abs(c01[9] + c02[9] + c03[9] - 0.75) > tol) {
      print("FAIL A total final=", c01[9] + c02[9] + c03[9], " expected 0.75");
      n_failures += 1;
    }

    // monotone non-decreasing
    for (t in 2:max_t + 1) {
      if (c01[t] < c01[t - 1] - tol || c02[t] < c02[t - 1] - tol || c03[t] < c03[t - 1] - tol) {
        print("FAIL A non-monotone at t=", t);
        n_failures += 1;
      }
    }
  }

  // ===========================================================================
  // CIF-B: two 0→1 events at different times, one 0→2, one 0→3
  //
  //   j=1: 0→1 at index 2  (pfs=2, rc=0, os=6, os_cens=0)
  //   j=2: 0→1 at index 4  (pfs=4, rc=0, os=8, os_cens=0)
  //   j=3: 0→2 at index 4  (pfs=4, rc=0, os=4, os_cens=0 — pfs==os)
  //   j=4: 0→3 at index 6  (pfs=6, rc=1, is_dropout=1)
  //
  //   cnt_01[2]=1, cnt_01[4]=1, cnt_02[4]=1, cnt_03[6]=1
  //   cif_01[2]=1/4, cif_01[4]=2/4=0.5 (cumulative)
  //   cif_02[4]=1/4
  //   cif_03[6]=1/4
  //   total at final = 4/4 = 1.0 (all patients accounted for)
  // ===========================================================================
  {
    int max_t = 8;
    array[4] int pfs            = {2, 4, 4, 6};
    array[4] int right_censored = {0, 0, 0, 1};
    array[4] int is_dropout     = {0, 0, 0, 1};
    array[4] int os             = {6, 8, 4, 6};
    array[4] int os_censored    = {0, 0, 0, 1};

    vector[max_t + 1] c01; vector[max_t + 1] c02; vector[max_t + 1] c03;
    (c01, c02, c03) = compute_trial_cif(pfs, right_censored, is_dropout, os, os_censored, max_t);

    real q = 0.25;
    real h = 0.5;

    // cif_01: 1/4 at index 2, cumulates to 2/4 at index 4
    if (abs(c01[1] - 0.0) > tol) { print("FAIL B cif_01[1]=", c01[1]); n_failures += 1; }
    if (abs(c01[2] - q)   > tol) { print("FAIL B cif_01[2]=", c01[2]); n_failures += 1; }
    if (abs(c01[3] - q)   > tol) { print("FAIL B cif_01[3]=", c01[3]); n_failures += 1; }
    if (abs(c01[4] - h)   > tol) { print("FAIL B cif_01[4]=", c01[4]); n_failures += 1; }
    if (abs(c01[9] - h)   > tol) { print("FAIL B cif_01[9]=", c01[9]); n_failures += 1; }

    // cif_02: 1/4 at index 4 (tie with 0→1 — classified as 0→2 because pfs==os)
    if (abs(c02[3] - 0.0) > tol) { print("FAIL B cif_02[3]=", c02[3]); n_failures += 1; }
    if (abs(c02[4] - q)   > tol) { print("FAIL B cif_02[4]=", c02[4]); n_failures += 1; }

    // cif_03: 1/4 at index 6
    if (abs(c03[5] - 0.0) > tol) { print("FAIL B cif_03[5]=", c03[5]); n_failures += 1; }
    if (abs(c03[6] - q)   > tol) { print("FAIL B cif_03[6]=", c03[6]); n_failures += 1; }

    // total at final = 1.0 (all 4 patients contributed)
    if (abs(c01[9] + c02[9] + c03[9] - 1.0) > tol) {
      print("FAIL B total final=", c01[9] + c02[9] + c03[9], " expected 1.0");
      n_failures += 1;
    }
  }

  // ===========================================================================
  // CIF-C: all 4 patients fully censored — all CIF vectors stay at zero
  //
  //   j=1..4: pfs=9, rc=1, is_dropout=0, pfs > max_t=8 → no contributions
  // ===========================================================================
  {
    int max_t = 8;
    array[4] int pfs            = {9, 9, 9, 9};
    array[4] int right_censored = {1, 1, 1, 1};
    array[4] int is_dropout     = {0, 0, 0, 0};
    array[4] int os             = {9, 9, 9, 9};
    array[4] int os_censored    = {1, 1, 1, 1};

    vector[max_t + 1] c01; vector[max_t + 1] c02; vector[max_t + 1] c03;
    (c01, c02, c03) = compute_trial_cif(pfs, right_censored, is_dropout, os, os_censored, max_t);

    for (t in 1:max_t + 1) {
      if (abs(c01[t]) > tol) { print("FAIL C cif_01[", t, "]=", c01[t]); n_failures += 1; }
      if (abs(c02[t]) > tol) { print("FAIL C cif_02[", t, "]=", c02[t]); n_failures += 1; }
      if (abs(c03[t]) > tol) { print("FAIL C cif_03[", t, "]=", c03[t]); n_failures += 1; }
    }
  }

  // ===========================================================================
  // CIF-D: admin-censored patient with pfs<=max_t is NOT counted as dropout
  //
  //   j=1: 0→1 at index 3  (pfs=3, rc=0, is_dropout=0)
  //   j=2: 0→3 at index 5  (pfs=5, rc=1, is_dropout=1)
  //   j=3: admin-censored within horizon (pfs=4, rc=1, is_dropout=0, pfs<=max_t=8)
  //   j=4: admin-censored beyond horizon (pfs=9, rc=1, is_dropout=0, pfs>max_t=8)
  //
  //   Old buggy code would have classified j=3 as dropout; correct code does NOT.
  //   cnt_01[3]=1, cnt_03[5]=1; j=3 and j=4 contribute nothing.
  //   total at max = 2/4 = 0.5
  // ===========================================================================
  {
    int max_t = 8;
    array[4] int pfs            = {3, 5, 4, 9};
    array[4] int right_censored = {0, 1, 1, 1};
    array[4] int is_dropout     = {0, 1, 0, 0};
    array[4] int os             = {7, 5, 4, 9};
    array[4] int os_censored    = {0, 1, 1, 1};

    vector[max_t + 1] c01; vector[max_t + 1] c02; vector[max_t + 1] c03;
    (c01, c02, c03) = compute_trial_cif(pfs, right_censored, is_dropout, os, os_censored, max_t);

    real q = 0.25;

    // cif_01: 1/4 at index 3
    if (abs(c01[2] - 0.0) > tol) { print("FAIL D cif_01[2]=", c01[2]); n_failures += 1; }
    if (abs(c01[3] - q)   > tol) { print("FAIL D cif_01[3]=", c01[3]); n_failures += 1; }
    if (abs(c01[9] - q)   > tol) { print("FAIL D cif_01[9]=", c01[9]); n_failures += 1; }

    // cif_02: zero everywhere (no direct deaths)
    for (t in 1:max_t + 1) {
      if (abs(c02[t]) > tol) { print("FAIL D cif_02[", t, "]=", c02[t]); n_failures += 1; }
    }

    // cif_03: 1/4 at index 5 (j=2 only; j=3 admin-censored should NOT contribute)
    if (abs(c03[4] - 0.0) > tol) { print("FAIL D cif_03[4]=", c03[4]); n_failures += 1; }
    if (abs(c03[5] - q)   > tol) { print("FAIL D cif_03[5]=", c03[5]); n_failures += 1; }
    if (abs(c03[9] - q)   > tol) { print("FAIL D cif_03[9]=", c03[9]); n_failures += 1; }

    // total at final = 0.5 (j=3 and j=4 admin-censored, no CIF contribution)
    if (abs(c01[9] + c02[9] + c03[9] - 0.5) > tol) {
      print("FAIL D total final=", c01[9] + c02[9] + c03[9], " expected 0.5");
      n_failures += 1;
    }
  }
}
