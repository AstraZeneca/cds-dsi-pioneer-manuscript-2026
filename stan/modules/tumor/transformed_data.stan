// ============================================================================
// TUMOR/SLD TRANSFORMED DATA
// ============================================================================

// Limit of detection for SLD measurements (cm)
real log_lod = log(0.1);

// Variables for handling separate baseline and proportional hazards
// int<lower = 1> n_tumor_separate_trials = separate_trial_tumor_gp ? n_trials : 1;
//
// print("n_tumor_separate_trials = ", n_tumor_separate_trials);

vector[sum(n_patient_visits)] log_sum_tumor_size = log(sum_tumor_size); // cm
vector[sum(n_patient_visits) - sum(n_patient_screening_visits)] post_treat_sld;

// Population indices ////

int<lower = 1> last_predict_visit = max_all_t; // max(pop_unique_visits);

print("last_predict_visit = ", last_predict_visit);

// Patient indices ////

// Note: patient2pop_unique_visit_idx maps patient visit indices to population-level unique visit times
// This is needed for hierarchical GP modeling

// Compute unique visits at population level (needed by _sf_transformed_data.stan for visit_cumsum_mat)
int n_pop_unique_visits = num_unique(t_patient_visits, 0);
array[n_pop_unique_visits] int pop_unique_visits = unique(t_patient_visits, 0);

print("n_pop_unique_visits = ", n_pop_unique_visits);
print("pop_unique_visits = ", pop_unique_visits);

// Map patient visits to population unique visits
array[sum(n_patient_visits)] int<lower = 1> patient2pop_unique_visit_idx =
  get_level2level_idx(pop_unique_visits, t_patient_visits, patient_visit_pos);

array[n_patients] int<lower = 1> patient_last_obs_visit = get_max_pos(t_patient_visits, patient_visit_pos);
array[n_patients] int<lower = 0, upper = last_predict_visit> n_patient_forecast_visits;

// Extract post-treatment SLD values
{
  int post_treat_pos = 1;

  for (s in 1:n_trials) {
    int curr_patient_pos, curr_patient_end;
    (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);

    for (i in curr_patient_pos:curr_patient_end) {
      n_patient_forecast_visits[i] = last_predict_visit - patient_last_obs_visit[i];

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
}

array[n_patients + 1] int<lower = 1> forecast_visits_pos = create_pos(n_patient_forecast_visits);

// Array of measurement times used for GP modeling
array[max_all_t] real all_tumor_measure_t = linspaced_array(max_all_t, 1, max_all_t);

// ============================================================================
// BASELINE SLD AND NORMALIZATION
// ============================================================================
// Patient-level baseline SLD (in cm) for converting normalized states to absolute SLD
// Used by tumor likelihood and by other_events module for tumor covariates

vector[n_patients] log_baseline_sld;

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
