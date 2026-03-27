// modules/multistate/likelihood.stan
// Multistate loglikelihood. Included inside each model's model block.
if (fit_multistate_data) {
  profile("multistate loglik") {
    ms_final_state[forecast_patient_idx] ~ multistate(
        likelihood_weight[forecast_patient_idx],
        enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
        enable_ms_03, enable_ms_32,
        ms_time_01[forecast_patient_idx], ms_time_02[forecast_patient_idx], ms_time_12[forecast_patient_idx],
        ms_time_03[forecast_patient_idx], ms_time_32[forecast_patient_idx],
        ms_censored_01[forecast_patient_idx], ms_censored_02[forecast_patient_idx], ms_censored_12[forecast_patient_idx],
        ms_censored_32[forecast_patient_idx],
        ms_prog_deterministic[forecast_patient_idx],
        ms_ic_gap_01[forecast_patient_idx],
        t_patient_visits,
        patient_visit_pos,
        log_cond_surv_01,
        log_cond_surv_02,
        log_cond_surv_12_s,
        log_cond_surv_12_t,
        log_cond_surv_03,
        log_cond_surv_32,
        enable_ms_visit_gated_01
      );
  }
}
