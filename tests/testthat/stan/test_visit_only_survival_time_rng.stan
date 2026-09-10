// Test for visit_only_survival_time_rng (unconditional and conditional overloads)
//
// Exercises:
//   1. Unconditional: events can only occur at visit weeks
//   2. Conditional: respects obs_time, only evaluates future visits
//   3. Conditional with right_censored=0: passes through observed event
//   4. High survival probability: most draws censored at max_all_t
//   5. Low survival probability: events cluster at first visit

functions {
  #include "multistate.stanfunctions"
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

data {
  int<lower=1> N_DRAWS;
}

generated quantities {
  int max_t = 100;

  // ── Case 1: Unconditional, moderate hazard at visit weeks ────────────────
  // Visits at weeks 6, 12, 18, 24, 30 (typical 6-week schedule)
  // log_cond_prob_surv = log(0.8) everywhere → ~20% hazard per visit
  {
    row_vector[max_t] log_surv = rep_row_vector(log(0.8), max_t);
    array[5] int visits = {6, 12, 18, 24, 30};

    array[N_DRAWS] int time_out;
    array[N_DRAWS] int cens_out;
    for (j in 1:N_DRAWS) {
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t);
      time_out[j] = r.1;
      cens_out[j] = r.2;
    }
  }

  // ── Case 2: Events can only land on visit weeks ──────────────────────────
  // Near-zero survival → event at first visit (week 10)
  // Verify all event times == 10
  array[N_DRAWS] int case2_time;
  array[N_DRAWS] int case2_cens;
  {
    row_vector[max_t] log_surv = rep_row_vector(log(1e-5), max_t);
    array[3] int visits = {10, 20, 30};

    for (j in 1:N_DRAWS) {
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t);
      case2_time[j] = r.1;
      case2_cens[j] = r.2;
    }
  }

  // ── Case 3: Very high survival → censored at max_t ──────────────────────
  array[N_DRAWS] int case3_time;
  array[N_DRAWS] int case3_cens;
  {
    row_vector[max_t] log_surv = rep_row_vector(log(1 - 1e-5), max_t);
    array[3] int visits = {10, 20, 30};

    for (j in 1:N_DRAWS) {
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t);
      case3_time[j] = r.1;
      case3_cens[j] = r.2;
    }
  }

  // ── Case 4: Conditional overload — observed event passed through ─────────
  array[N_DRAWS] int case4_time;
  array[N_DRAWS] int case4_cens;
  {
    row_vector[max_t] log_surv = rep_row_vector(log(0.5), max_t);
    array[3] int visits = {10, 20, 30};

    for (j in 1:N_DRAWS) {
      // right_censored = 0 → pass through observed event at week 15
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t, 15, 0);
      case4_time[j] = r.1;
      case4_cens[j] = r.2;
    }
  }

  // ── Case 5: Conditional overload — forecast from obs_time, skip past visits ─
  // obs_time = 12, visits at {6, 12, 18, 24, 30}
  // Should skip visits 6 and 12, only evaluate 18, 24, 30
  // Near-zero survival → event at first evaluated visit (week 18)
  array[N_DRAWS] int case5_time;
  array[N_DRAWS] int case5_cens;
  {
    row_vector[max_t] log_surv = rep_row_vector(log(1e-5), max_t);
    array[5] int visits = {6, 12, 18, 24, 30};

    for (j in 1:N_DRAWS) {
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t, 12, 1);
      case5_time[j] = r.1;
      case5_cens[j] = r.2;
    }
  }

  // ── Case 6: Conditional, high survival, all future visits survived → censored ─
  array[N_DRAWS] int case6_time;
  array[N_DRAWS] int case6_cens;
  {
    row_vector[max_t] log_surv = rep_row_vector(log(1 - 1e-5), max_t);
    array[3] int visits = {10, 20, 30};

    for (j in 1:N_DRAWS) {
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t, 5, 1);
      case6_time[j] = r.1;
      case6_cens[j] = r.2;
    }
  }

  // ── Case 7: Events only at visit weeks (moderate hazard, many draws) ─────
  // Track which weeks have events to verify no non-visit-week events
  array[N_DRAWS] int case7_time;
  array[N_DRAWS] int case7_cens;
  {
    row_vector[max_t] log_surv = rep_row_vector(log(0.7), max_t);
    array[4] int visits = {8, 16, 24, 32};

    for (j in 1:N_DRAWS) {
      tuple(int, int) r = visit_only_survival_time_rng(log_surv, visits, max_t);
      case7_time[j] = r.1;
      case7_cens[j] = r.2;
    }
  }
}
