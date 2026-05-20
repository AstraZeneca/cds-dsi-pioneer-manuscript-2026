// tr/priors.stan — active priors for total rate module
// Multi-level hierarchy: loop over all levels for unified prior structure
// Uses enabled position arrays for efficient indexing into compacted parameter arrays

// Population intercept
tr_loc_pop ~ normal(tr_loc_pop_mean, tr_loc_pop_sd);

// Population covariate effects
if (enable_pop_cov_tr) {
  tr_coef_qr_pop ~ normal(tr_coef_qr_pop_mean, tr_coef_qr_pop_sd);
}

// ===== UNIFIED LOOP OVER ALL LEVELS =====
{
  int sd_idx = 0;
  for (lv in 1:n_levels) {
    int mode = enable_level_intercept_tr[lv];

    // SD prior for any level that samples its SD as a free parameter (RE or RE_CP)
    if (mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      tr_sd_level_intercept_raw[sd_idx] ~ normal(0, tr_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      tr_sd_level_slope[lv] ~ normal(0, tr_sd_level_slope_sd[lv]);
    }

    // RAW path: FE, RE, RE_GP all sample from std_normal / student_t(nu, 0, 1)
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_tr_intercept[lv];
      int r_hi = raw_level_pos_tr_intercept[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          tr_raw_level_intercept[r_lo:r_hi] ~ student_t(tr_nu_level[lv], 0, 1);
        else
          tr_raw_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }

    // CP path: RE_CP samples directly at scale sd
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_tr_intercept[lv];
      int c_hi = cp_level_pos_tr_intercept[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          tr_cp_level_intercept[c_lo:c_hi]
            ~ student_t(tr_nu_level[lv], 0, tr_sd_intercept_pergroup[lv][1:(c_hi - c_lo + 1)]);
        else
          tr_cp_level_intercept[c_lo:c_hi]
            ~ normal(0, tr_sd_intercept_pergroup[lv][1:(c_hi - c_lo + 1)]);
      }
    }

    // Slope effects — parameterization follows the level's intercept mode
    if (enable_level_cov_tr[lv] && n_covar > 0) {
      // RAW path (FE, RE, RE_GP)
      if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
        int r_lo = raw_level_pos_tr_slope[lv];
        int r_hi = raw_level_pos_tr_slope[lv + 1] - 1;
        if (r_hi >= r_lo) {
          if (enable_student_t_hierarchy)
            to_vector(tr_raw_level_slope[r_lo:r_hi, :]) ~ student_t(tr_nu_level[lv], 0, 1);
          else
            to_vector(tr_raw_level_slope[r_lo:r_hi, :]) ~ std_normal();
        }
      }

      // CP path (RE_CP): sample at scale sd_level_slope[lv]
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_tr_slope[lv];
        int c_hi = cp_level_pos_tr_slope[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_covar) {
            if (enable_student_t_hierarchy)
              tr_cp_level_slope[c_lo:c_hi, k]
                ~ student_t(tr_nu_level[lv], 0, tr_sd_level_slope[lv, k]);
            else
              tr_cp_level_slope[c_lo:c_hi, k]
                ~ normal(0, tr_sd_level_slope[lv, k]);
          }
        }
      }
    }
  }
}

// Student-t nu priors (only when enabled)
if (enable_student_t_hierarchy) {
  for (lv in 1:n_levels)
    tr_nu_level[lv] ~ gamma(tr_nu_level_prior_alpha[lv], tr_nu_level_prior_beta[lv]);
}

// Patient-level process noise priors - only when feature is enabled
if (enable_patient_process_noise_tr) {
	to_vector(tr_raw_patient_process_noise) ~ std_normal();
	tr_log_sd_pop_process_noise[1] ~ normal(tr_log_sd_pop_process_noise_mean, tr_log_sd_pop_process_noise_sd);
	tr_logit_phi_pop_process_noise[1] ~ normal(tr_logit_phi_pop_process_noise_mean, tr_logit_phi_pop_process_noise_sd);
	tr_sd_patient_log_sd_process_noise[1] ~ normal(0, tr_log_sd_patient_process_noise_sd);
	tr_sd_patient_phi_process_noise[1] ~ normal(0, tr_phi_patient_process_noise_sd);
	tr_raw_patient_log_sd_process_noise ~ std_normal();
	tr_raw_patient_phi_process_noise ~ std_normal();
}

// Population-level time-varying process noise priors
if (enable_pop_process_noise_tr) {
	tr_raw_pop_process_noise ~ std_normal();
	tr_log_sd_pop_process_noise_pop[1] ~ normal(tr_log_sd_pop_process_noise_pop_mean, tr_log_sd_pop_process_noise_pop_sd);
	tr_logit_phi_pop_process_noise_pop[1] ~ normal(tr_logit_phi_pop_process_noise_pop_mean, tr_logit_phi_pop_process_noise_pop_sd);
}

// SD sub-hierarchy priors (inert when mode matrix is all-NONE)
#include "modules/tr/_sd_subhierarchy_priors.stan"
