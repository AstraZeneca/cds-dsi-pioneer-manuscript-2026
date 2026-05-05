// Tests for estimate_kaplan_meier, km_quantiles, and calc_km_pfs_n
//
// All expected values are derived by hand from the KM formula and verified
// against the algorithm in pfs.stanfunctions.
//
// Four KM scenarios + quantile extraction + PFS-n evaluation.
// Output: n_failures (count of failed assertions). R test asserts n_failures == 0.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "pfs.stanfunctions"
}

data {
  int<lower=0> dummy;
}

generated quantities {
  int n_failures = 0;
  real tol = 1e-10;

  // ===========================================================================
  // KM-1: 5 patients, mixed events and censoring
  //   event_time = {2,3,5,4,7}, right_censored = {0,0,0,1,1}, max_t = 7
  //
  //   Sorted order: t=2(evt), t=3(evt), t=4(cens), t=5(evt), t=7(cens)
  //   n=5: S(2)=0.8, S(3)=0.6, S(4)=0.6, S(5)=0.3, S(6..7)=0.3
  // ===========================================================================
  {
    array[5] int event_time = {2, 3, 5, 4, 7};
    array[5] int right_censored = {0, 0, 0, 1, 1};
    int max_t = 7;

    vector[max_t + 1] s;
    array[max_t + 1] int at_risk_arr;
    array[max_t + 1] int nrc;
    array[max_t + 1] int ne;
    (s, at_risk_arr, nrc, ne) = estimate_kaplan_meier(event_time, right_censored, max_t);

    // Expected S(0..7): indices 1..8
    vector[8] exp_s = [1.0, 1.0, 0.8, 0.6, 0.6, 0.3, 0.3, 0.3]';
    for (t in 1:max_t + 1) {
      if (abs(s[t] - exp_s[t]) > tol) {
        print("FAIL KM1 S[", t, "]: got ", s[t], " expected ", exp_s[t]);
        n_failures += 1;
      }
    }

    // at_risk: [5,5,5,4,3,2,1,1] (n at start of each t=0..7)
    array[8] int exp_ar = {5, 5, 5, 4, 3, 2, 1, 1};
    for (t in 1:max_t + 1) {
      if (at_risk_arr[t] != exp_ar[t]) {
        print("FAIL KM1 at_risk[", t, "]: got ", at_risk_arr[t], " expected ", exp_ar[t]);
        n_failures += 1;
      }
    }

    // n_exited: events at each t; censoring at t=4(idx5) and t=7(idx8)
    array[8] int exp_ne  = {0, 0, 1, 1, 0, 1, 0, 0};
    array[8] int exp_nrc = {0, 0, 0, 0, 1, 0, 0, 1};
    for (t in 1:max_t + 1) {
      if (ne[t] != exp_ne[t]) {
        print("FAIL KM1 n_exited[", t, "]: got ", ne[t], " expected ", exp_ne[t]);
        n_failures += 1;
      }
      if (nrc[t] != exp_nrc[t]) {
        print("FAIL KM1 n_right_censored[", t, "]: got ", nrc[t], " expected ", exp_nrc[t]);
        n_failures += 1;
      }
    }

    // km_quantiles: median (p=0.5) on KM-1
    //   S drops from 0.6 at index 5 to 0.3 at index 6 (t=4 to t=5)
    //   weight = (0.6-0.5)/(0.6-0.3) = 1/3  → quantile = (6-2) + 1/3 = 13/3
    {
      vector[1] qv; array[1] int qi;
      (qv, qi) = km_quantiles(s, [0.5]');
      real exp_median = 4.0 + 1.0/3.0;
      if (abs(qv[1] - exp_median) > tol) {
        print("FAIL KM1 median: got ", qv[1], " expected ", exp_median);
        n_failures += 1;
      }
      if (qi[1] != 0) {
        print("FAIL KM1 median cannot_calculate: got ", qi[1], " expected 0");
        n_failures += 1;
      }
    }

    // km_quantiles: Q3 (p=0.75) — crosses between S[3]=0.8 and S[4]=0.6 (t=2 to t=3)
    //   weight = (0.8-0.75)/(0.8-0.6) = 0.25 → quantile = (4-2) + 0.25 = 2.25
    {
      vector[1] qv; array[1] int qi;
      (qv, qi) = km_quantiles(s, [0.75]');
      real exp_q3 = 2.25;
      if (abs(qv[1] - exp_q3) > tol) {
        print("FAIL KM1 Q3: got ", qv[1], " expected ", exp_q3);
        n_failures += 1;
      }
      if (qi[1] != 0) {
        print("FAIL KM1 Q3 cannot_calculate: got ", qi[1], " expected 0");
        n_failures += 1;
      }
    }

    // km_quantiles: Q1 (p=0.25) — S never drops below 0.25 → cannot_calculate=1
    {
      vector[1] qv; array[1] int qi;
      (qv, qi) = km_quantiles(s, [0.25]');
      if (qi[1] != 1) {
        print("FAIL KM1 Q1 cannot_calculate: got ", qi[1], " expected 1 (S stays >= 0.3)");
        n_failures += 1;
      }
    }

    // calc_km_pfs_n at exact integer timepoints and fractional
    {
      real pfs4  = calc_km_pfs_n(s, 4.0);   // S[5]=0.6
      real pfs5  = calc_km_pfs_n(s, 5.0);   // S[6]=0.3
      real pfs45 = calc_km_pfs_n(s, 4.5);   // 0.5*S[5]+0.5*S[6]=0.45
      real pfs0  = calc_km_pfs_n(s, 0.0);   // S[1]=1.0
      real pfs7  = calc_km_pfs_n(s, 7.0);   // S[8]=0.3 (at max_t)
      if (abs(pfs4  - 0.6)  > tol) { print("FAIL KM1 pfs_n(4.0)=", pfs4,  " expected 0.6");  n_failures += 1; }
      if (abs(pfs5  - 0.3)  > tol) { print("FAIL KM1 pfs_n(5.0)=", pfs5,  " expected 0.3");  n_failures += 1; }
      if (abs(pfs45 - 0.45) > tol) { print("FAIL KM1 pfs_n(4.5)=", pfs45, " expected 0.45"); n_failures += 1; }
      if (abs(pfs0  - 1.0)  > tol) { print("FAIL KM1 pfs_n(0.0)=", pfs0,  " expected 1.0");  n_failures += 1; }
      if (abs(pfs7  - 0.3)  > tol) { print("FAIL KM1 pfs_n(7.0)=", pfs7,  " expected 0.3");  n_failures += 1; }
    }
  }

  // ===========================================================================
  // KM-2: 3 patients, all events (no censoring)
  //   event_time = {1,3,6}, right_censored = {0,0,0}, max_t = 7
  //
  //   S(1)=2/3, S(3)=1/3, S(6)=0; S constant between events
  // ===========================================================================
  {
    array[3] int event_time = {1, 3, 6};
    array[3] int right_censored = {0, 0, 0};
    int max_t = 7;

    vector[max_t + 1] s;
    array[max_t + 1] int at_risk_arr;
    array[max_t + 1] int nrc;
    array[max_t + 1] int ne;
    (s, at_risk_arr, nrc, ne) = estimate_kaplan_meier(event_time, right_censored, max_t);

    // Expected S: [1, 2/3, 2/3, 1/3, 1/3, 1/3, 0, 0]
    real two_thirds = 2.0 / 3.0;
    real one_third  = 1.0 / 3.0;
    vector[8] exp_s = [1.0, two_thirds, two_thirds, one_third, one_third, one_third, 0.0, 0.0]';
    for (t in 1:max_t + 1) {
      if (abs(s[t] - exp_s[t]) > tol) {
        print("FAIL KM2 S[", t, "]: got ", s[t], " expected ", exp_s[t]);
        n_failures += 1;
      }
    }

    // S is monotone non-increasing
    for (t in 2:max_t + 1) {
      if (s[t] > s[t - 1] + tol) {
        print("FAIL KM2 non-monotone at t=", t);
        n_failures += 1;
      }
    }

    // km_quantiles: median (p=0.5)
    //   Crosses from S[3]=2/3 to S[4]=1/3 between Stan indices 3→4 (t=2→3 in 0-based)
    //   weight = (2/3 - 0.5) / (2/3 - 1/3) = (1/6) / (1/3) = 0.5
    //   quantile = (4-2) + 0.5 = 2.5
    {
      vector[1] qv; array[1] int qi;
      (qv, qi) = km_quantiles(s, [0.5]');
      real exp_median = 2.5;
      if (abs(qv[1] - exp_median) > tol) {
        print("FAIL KM2 median: got ", qv[1], " expected ", exp_median);
        n_failures += 1;
      }
      if (qi[1] != 0) {
        print("FAIL KM2 median cannot_calculate: got ", qi[1], " expected 0");
        n_failures += 1;
      }
    }

    // calc_km_pfs_n
    {
      real pfs1 = calc_km_pfs_n(s, 1.0);  // S[2]=2/3
      real pfs3 = calc_km_pfs_n(s, 3.0);  // S[4]=1/3
      real pfs6 = calc_km_pfs_n(s, 6.0);  // S[7]=0
      if (abs(pfs1 - two_thirds) > tol) { print("FAIL KM2 pfs_n(1.0)=", pfs1, " expected 2/3"); n_failures += 1; }
      if (abs(pfs3 - one_third)  > tol) { print("FAIL KM2 pfs_n(3.0)=", pfs3, " expected 1/3"); n_failures += 1; }
      if (abs(pfs6 - 0.0)        > tol) { print("FAIL KM2 pfs_n(6.0)=", pfs6, " expected 0");   n_failures += 1; }
    }
  }

  // ===========================================================================
  // KM-3: 3 patients, all right-censored — S stays at 1.0 throughout
  //   event_time = {3,5,7}, right_censored = {1,1,1}, max_t = 7
  // ===========================================================================
  {
    array[3] int event_time = {3, 5, 7};
    array[3] int right_censored = {1, 1, 1};
    int max_t = 7;

    vector[max_t + 1] s;
    array[max_t + 1] int at_risk_arr;
    array[max_t + 1] int nrc;
    array[max_t + 1] int ne;
    (s, at_risk_arr, nrc, ne) = estimate_kaplan_meier(event_time, right_censored, max_t);

    // All S = 1.0
    for (t in 1:max_t + 1) {
      if (abs(s[t] - 1.0) > tol) {
        print("FAIL KM3 S[", t, "] != 1.0: got ", s[t]);
        n_failures += 1;
      }
      if (ne[t] != 0) {
        print("FAIL KM3 n_exited[", t, "] != 0: got ", ne[t]);
        n_failures += 1;
      }
    }

    // km_quantiles: p=0.5 cannot be reached (S never drops below 1.0 - epsilon)
    {
      vector[1] qv; array[1] int qi;
      (qv, qi) = km_quantiles(s, [0.5]');
      if (qi[1] != 1) {
        print("FAIL KM3 cannot_calculate for all-censored: got ", qi[1], " expected 1");
        n_failures += 1;
      }
    }
  }

  // ===========================================================================
  // KM-4: 3 patients, tied events at t=4 — S drops to exactly 0
  //   event_time = {4,4,4}, right_censored = {0,0,0}, max_t = 6
  //
  //   S(0..3)=1.0, S(4)=(3-3)/3=0, S(5..6)=0
  // ===========================================================================
  {
    array[3] int event_time = {4, 4, 4};
    array[3] int right_censored = {0, 0, 0};
    int max_t = 6;

    vector[max_t + 1] s;
    array[max_t + 1] int at_risk_arr;
    array[max_t + 1] int nrc;
    array[max_t + 1] int ne;
    (s, at_risk_arr, nrc, ne) = estimate_kaplan_meier(event_time, right_censored, max_t);

    // Expected S: [1, 1, 1, 1, 0, 0, 0] (indices 1..7, t=0..6)
    // S(t=4) = (3-3)/3 = 0 — all events fire at t=4
    vector[7] exp_s = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]';
    for (t in 1:max_t + 1) {
      if (abs(s[t] - exp_s[t]) > tol) {
        print("FAIL KM4 S[", t, "]: got ", s[t], " expected ", exp_s[t]);
        n_failures += 1;
      }
    }

    // All 3 events occur at t=4 (index 5)
    if (ne[5] != 3) {
      print("FAIL KM4 n_exited at t=4: got ", ne[5], " expected 3");
      n_failures += 1;
    }

    // km_quantiles: median (p=0.5) — S drops from 1 to 0 at index 5 (t=4)
    //   weight = (1-0.5)/(1-0) = 0.5 → quantile = (5-2) + 0.5 = 3.5
    {
      vector[1] qv; array[1] int qi;
      (qv, qi) = km_quantiles(s, [0.5]');
      real exp_median = 3.5;
      if (abs(qv[1] - exp_median) > tol) {
        print("FAIL KM4 median: got ", qv[1], " expected ", exp_median);
        n_failures += 1;
      }
      if (qi[1] != 0) {
        print("FAIL KM4 median cannot_calculate: got ", qi[1], " expected 0");
        n_failures += 1;
      }
    }

    // calc_km_pfs_n at t=4 and t=5
    {
      real pfs4 = calc_km_pfs_n(s, 4.0);  // S[5]=0.0 (all events fired at t=4)
      real pfs5 = calc_km_pfs_n(s, 5.0);  // S[6]=0.0
      if (abs(pfs4 - 0.0) > tol) { print("FAIL KM4 pfs_n(4.0)=", pfs4, " expected 0.0"); n_failures += 1; }
      if (abs(pfs5 - 0.0) > tol) { print("FAIL KM4 pfs_n(5.0)=", pfs5, " expected 0.0"); n_failures += 1; }
    }
  }

  // ===========================================================================
  // km_quantiles: multiple quantiles simultaneously on KM-1 curve
  //   p = [0.25, 0.5, 0.75] on S = [1,1,0.8,0.6,0.6,0.3,0.3,0.3]
  //   Q1 (0.25): cannot_calculate=1 (S never reaches 0.25)
  //   Q2 (0.50): 4 + 1/3
  //   Q3 (0.75): 2.25
  // ===========================================================================
  {
    array[5] int event_time = {2, 3, 5, 4, 7};
    array[5] int right_censored = {0, 0, 0, 1, 1};
    int max_t = 7;

    vector[max_t + 1] s;
    array[max_t + 1] int ar;
    array[max_t + 1] int nrc;
    array[max_t + 1] int ne;
    (s, ar, nrc, ne) = estimate_kaplan_meier(event_time, right_censored, max_t);

    vector[3] qv; array[3] int qi;
    (qv, qi) = km_quantiles(s, [0.25, 0.5, 0.75]');

    // p=0.25: cannot_calculate
    if (qi[1] != 1) {
      print("FAIL KM-multi Q1 cannot_calculate=", qi[1], " expected 1");
      n_failures += 1;
    }
    // p=0.5: quantile = 4+1/3
    if (qi[2] != 0 || abs(qv[2] - (4.0 + 1.0/3.0)) > tol) {
      print("FAIL KM-multi Q2: got (", qv[2], ",", qi[2], ") expected (4.333..,0)");
      n_failures += 1;
    }
    // p=0.75: quantile = 2.25
    if (qi[3] != 0 || abs(qv[3] - 2.25) > tol) {
      print("FAIL KM-multi Q3: got (", qv[3], ",", qi[3], ") expected (2.25,0)");
      n_failures += 1;
    }
  }
}
