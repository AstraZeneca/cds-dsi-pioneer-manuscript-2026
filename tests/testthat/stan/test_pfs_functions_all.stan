functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> N;
  array[N_CASES, N] int pfs;
  array[N_CASES, N] int right_censored;
  array[N_CASES, N] int max_time;
  // expected outputs for each case
  array[N_CASES, N] int expected_pfs;
  array[N_CASES, N] int expected_right_censored;
}

generated quantities {
  array[N_CASES, N] int out_pfs;
  array[N_CASES, N] int out_right_censored;
  array[N_CASES] int n_match;
  for (case in 1:N_CASES) {
    tuple(array[N] int, array[N] int) res = truncate_at_max_time(pfs[case], right_censored[case], max_time[case]);
    out_pfs[case] = res.1;
    out_right_censored[case] = res.2;
    n_match[case] = 1;
    for (i in 1:N) {
      if (out_pfs[case, i] != expected_pfs[case, i] || out_right_censored[case, i] != expected_right_censored[case, i]) {
        n_match[case] = 0;
      }
    }
  }
  // Test cases for truncate_at_max_time will go here
  // Test 1: No truncation needed (all pfs < max_time)
  {
    int pfs[3] = {5, 10, 15};
    int right_censored[3] = {0, 0, 0};
    int max_time[3] = {20, 20, 20};
    array[3] int pfs_trunc;
    array[3] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 1: No truncation needed");
    print(pfs_trunc); // expect {5, 10, 15}
    print(cens_trunc); // expect {0, 0, 0}
  }

  // Test 2: All pfs truncated (all pfs > max_time)
  {
    int pfs[3] = {25, 30, 40};
    int right_censored[3] = {0, 0, 0};
    int max_time[3] = {20, 20, 20};
    array[3] int pfs_trunc;
    array[3] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 2: All pfs truncated");
    print(pfs_trunc); // expect {20, 20, 20}
    print(cens_trunc); // expect {1, 1, 1}
  }

  // Test 3: Mixed truncation and censoring
  {
    int pfs[4] = {5, 15, 25, 35};
    int right_censored[4] = {0, 1, 0, 1};
    int max_time[4] = {10, 20, 20, 30};
    array[4] int pfs_trunc;
    array[4] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 3: Mixed truncation and censoring");
    print(pfs_trunc); // expect {5, 15, 20, 30}
    print(cens_trunc); // expect {0, 1, 1, 1}
  }

  // Test 4: Already censored, no truncation
  {
    int pfs[2] = {8, 12};
    int right_censored[2] = {1, 1};
    int max_time[2] = {20, 20};
    array[2] int pfs_trunc;
    array[2] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 4: Already censored, no truncation");
    print(pfs_trunc); // expect {8, 12}
    print(cens_trunc); // expect {1, 1}
  }

  // Test 5: Already censored, but pfs > max_time
  {
    int pfs[2] = {25, 30};
    int right_censored[2] = {1, 1};
    int max_time[2] = {20, 20};
    array[2] int pfs_trunc;
    array[2] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 5: Already censored, but pfs > max_time");
    print(pfs_trunc); // expect {20, 20}
    print(cens_trunc); // expect {1, 1}
  }

  // Test 6: Edge case, pfs == max_time
  {
    int pfs[3] = {10, 20, 30};
    int right_censored[3] = {0, 0, 0};
    int max_time[3] = {10, 20, 30};
    array[3] int pfs_trunc;
    array[3] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 6: Edge case, pfs == max_time");
    print(pfs_trunc); // expect {10, 20, 30}
    print(cens_trunc); // expect {0, 0, 0}
  }

  // Test 7: Single patient, uncensored, needs truncation
  {
    int pfs[1] = {50};
    int right_censored[1] = {0};
    int max_time[1] = {20};
    array[1] int pfs_trunc;
    array[1] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 7: Single patient, uncensored, needs truncation");
    print(pfs_trunc); // expect {20}
    print(cens_trunc); // expect {1}
  }

  // Test 8: Single patient, already censored, no truncation
  {
    int pfs[1] = {10};
    int right_censored[1] = {1};
    int max_time[1] = {20};
    array[1] int pfs_trunc;
    array[1] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 8: Single patient, already censored, no truncation");
    print(pfs_trunc); // expect {10}
    print(cens_trunc); // expect {1}
  }

  // Test 9: Single patient, already censored, needs truncation
  {
    int pfs[1] = {30};
    int right_censored[1] = {1};
    int max_time[1] = {20};
    array[1] int pfs_trunc;
    array[1] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 9: Single patient, already censored, needs truncation");
    print(pfs_trunc); // expect {20}
    print(cens_trunc); // expect {1}
  }

  // Test 10: All zero times
  {
    int pfs[3] = {0, 0, 0};
    int right_censored[3] = {0, 1, 0};
    int max_time[3] = {0, 0, 0};
    array[3] int pfs_trunc;
    array[3] int cens_trunc;
    truncate_at_max_time(pfs, right_censored, max_time, pfs_trunc, cens_trunc);
    print("Test 10: All zero times");
    print(pfs_trunc); // expect {0, 0, 0}
    print(cens_trunc); // expect {1, 1, 1}
  }
}
