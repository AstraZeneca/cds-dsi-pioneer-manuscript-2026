// ============================================================================
// TUMOR/SLD TRANSFORMED DATA
// ============================================================================
// SLD-specific preprocessing. Visit infrastructure (forecast_visits_pos, etc.)
// is now in _base_transformed_data.stan.

// Limit of detection for SLD measurements (cm)
real log_lod = log(0.1);

vector[sum(n_patient_visits)] log_sum_tumor_size = log(sum_tumor_size); // cm
vector[sum(n_patient_visits) - sum(n_patient_screening_visits)] post_treat_sld;



// Extract post-treatment SLD values
{
  int post_treat_pos = 1;

  for (i in 1:n_patients) {
    int curr_patient_visits_pos, curr_patient_visits_end;
    (curr_patient_visits_pos, curr_patient_visits_end) = get_pos(patient_visit_pos, i);

    for (m in curr_patient_visits_pos:curr_patient_visits_end) {
      if (t_patient_visits[m] > 0) {
        post_treat_sld[post_treat_pos] = sum_tumor_size[m];
        post_treat_pos += 1;
      }
    }
  }
}

// Array of time points used for GP modeling in the multistate module
array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

// ============================================================================
// BASELINE SLD AND NORMALIZATION
// ============================================================================
// Patient-level baseline SLD (in cm) for converting normalized states to absolute SLD
// Used by tumor likelihood and by multistate module for tumor covariates

vector[n_patients] log_baseline_sld;

// ============================================================================
// SLD NORMALIZATION CONSTANTS FOR PROPORTIONAL HAZARDS
// ============================================================================
// Robust normalization statistics (median and IQR) from observed log(SLD)
// Using median/IQR instead of mean/SD for robustness to outliers
// Used by multistate module for time-varying covariates

real median_log_sld_obs;
real iqr_log_sld_obs;

// Compute robust normalization statistics from observed log(SLD)
{
  // Get all observed log(SLD) values across all patients and visits
  // Filter out zeros to avoid -Inf
  vector[sum(n_patient_visits)] log_sld_all_obs;
  int n_positive = 0;

  for (i in 1:sum(n_patient_visits)) {
    if (sum_tumor_size[i] > 0) {
      n_positive += 1;
      log_sld_all_obs[n_positive] = log(sum_tumor_size[i]);
    }
  }

  // Compute robust normalization using median and IQR from positive observations
  {
    array[3] real quantiles_obs = quantile(log_sld_all_obs[1:n_positive], {0.25, 0.5, 0.75});
    median_log_sld_obs = quantiles_obs[2];
    iqr_log_sld_obs = quantiles_obs[3] - quantiles_obs[1];
  }
}

// ============================================================================
// VELOCITY NORMALIZATION CONSTANTS (for the (level, velocity) coupling basis)
// ============================================================================
// Robust median/IQR of OBSERVED per-week log(SLD) deltas between consecutive
// measured visits. Anchors the scale of the latent velocity covariate
// (observed-for-scale, latent-for-signal — mirrors median_log_sld_obs). A visit
// is "measured" when sum_tumor_size > 0; deltas are normalized per-week to match
// the latent weekly-grid central difference (visits are irregularly spaced).
// Always computed; consumed only when enable_ms_velocity_basis = 1.

real median_velocity_obs;
real iqr_velocity_obs;
{
  // Upper bound on inter-visit pairs = total visits (each patient contributes
  // at most n_visits - 1 pairs; sum over patients <= sum(n_patient_visits)).
  vector[sum(n_patient_visits)] vel_all_obs;
  int n_deltas = 0;

  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    // Walk consecutive MEASURED visits; difference adjacent measured pairs.
    int prev_v = 0;  // 0 = no previous measured visit yet
    for (v in visit_start:visit_end) {
      if (sum_tumor_size[v] > 0) {
        if (prev_v > 0) {
          int dw = t_patient_visits[v] - t_patient_visits[prev_v];
          if (dw > 0) {  // guard against duplicate-week visits (dw=0)
            n_deltas += 1;
            vel_all_obs[n_deltas] =
              (log(sum_tumor_size[v]) - log(sum_tumor_size[prev_v])) * 1.0 / dw;
          }
        }
        prev_v = v;
      }
    }
  }

  if (n_deltas < 2) {
    // Degenerate: cannot form a robust IQR. Fall back to unit scale so the
    // standardized velocity is just centered raw velocity; flag loudly.
    print("WARNING: only ", n_deltas, " observed inter-visit velocity deltas; ",
          "median_velocity_obs/iqr_velocity_obs fall back to (0, 1).");
    median_velocity_obs = 0.0;
    iqr_velocity_obs = 1.0;
  } else {
    array[3] real q = quantile(vel_all_obs[1:n_deltas], {0.25, 0.5, 0.75});
    median_velocity_obs = q[2];
    iqr_velocity_obs = q[3] - q[1];
    if (iqr_velocity_obs <= 0)
      iqr_velocity_obs = 1.0;  // all deltas equal — avoid divide-by-zero
  }
}

// Normalize SLD by baseline for each patient (used by observation model)
vector<lower = 0>[sum(n_patient_visits)] normalized_sld;

{
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);
    log_baseline_sld[i] = log(sum_tumor_size[visit_start]);
    normalized_sld[visit_start:visit_end] = sum_tumor_size[visit_start:visit_end] / sum_tumor_size[visit_start];
  }
}

// Biomarker-agnostic baseline for state_space module
vector[n_patients] baseline_obs_per_patient;
{
  for (i in 1:n_patients) {
    int visit_start = patient_visit_pos[i];
    baseline_obs_per_patient[i] = sum_tumor_size[visit_start];
  }
}
