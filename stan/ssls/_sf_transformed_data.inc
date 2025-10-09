// Include fragment: transformed data block content for SSLS models
array[sf_rep_T] real rep_time_points;
if (sf_rep_T > 0) rep_time_points = linspaced_array(sf_rep_T, 1, sf_rep_T);
real<lower = 0> lod = 0.1; // cm
vector<lower = 0>[sum(n_patient_visits)] normalized_sld;
for (i in 1:n_patients) {
  int visit_pos, visit_end; (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
  normalized_sld[visit_pos:visit_end] = sum_tumor_size[visit_pos:visit_end] / sum_tumor_size[visit_pos]; 
}
int<lower = 1> n_total_visits_m1 = sum(n_patient_visits) - n_patients;
array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
int<lower = 1, upper = n_trials> n_train_trials = max(patient_trial) - min(patient_trial) + 1;
int<lower = 1> n_total_train_visits = sum(n_patient_visits);
int<lower = 1> n_total_train_visits_m1 = n_total_train_visits - n_patients;
int<lower = 1> n_total_train_forecast_visits = get_pos_total_size(forecast_visits_pos);
int<lower = 0, upper = n_patients> n_right_censored_patients = sum(right_censored);
int<lower = 0, upper = n_patients> n_right_uncensored_patients = n_patients - n_right_censored_patients;
array[n_right_uncensored_patients] int<lower = 1, upper = n_patients> right_uncensored_patients;
array[n_trials + 1] int trial_right_censored_pos, trial_right_uncensored_pos;
{
  int right_uncensored_idx = 1; array[n_trials] int n_trial_censored = zeros_int_array(n_trials), n_trial_uncensored = zeros_int_array(n_trials);
  for (i in 1:n_patients) {
    if (right_censored[i]) n_trial_censored[patient_trial[i]] += 1; else { right_uncensored_patients[right_uncensored_idx] = i; right_uncensored_idx += 1; n_trial_uncensored[patient_trial[i]] += 1; }
  }
  trial_right_censored_pos = create_pos(n_trial_censored);
  trial_right_uncensored_pos = create_pos(n_trial_uncensored);
}
array[n_patients] int<lower = 1> ub_pfs_p1; for (i in 1:n_patients) ub_pfs_p1[i] = pfs[i] + interval_censored[i] + 1;
real log_lod = log(0.1);
int CR = 1; int PR = 2; int SD = 3; int PD = 4;
int NT_CR = 1; int NT_STABLE = 2; int NT_PD = 3;
matrix[n_patients, n_covar] Q_covar_design_matrix = qr_thin_Q(covar_design_matrix) * sqrt(n_patients - 1);
matrix[n_covar, n_covar] R_covar_design_matrix = qr_thin_R(covar_design_matrix) / sqrt(n_patients - 1);
