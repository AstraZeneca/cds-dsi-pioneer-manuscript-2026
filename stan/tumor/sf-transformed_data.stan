array[sf_rep_T] real rep_time_points;

if (sf_rep_T > 0) { 
  rep_time_points = linspaced_array(sf_rep_T, 1, sf_rep_T);
}

real<lower = 0> lod = 0.1; // cm

vector<lower = 0>[sum(n_patient_visits)] normalized_sld;

for (i in 1:n_patients) {
  int visit_pos, visit_end;
  (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
  
  normalized_sld[visit_pos:visit_end] = sum_tumor_size[visit_pos:visit_end] / sum_tumor_size[visit_pos]; 
}

int<lower = 1> n_total_visits_m1 = sum(n_patient_visits) - n_patients;
array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);

int<lower = 1, upper = n_trials> n_train_trials = max(patient_trial) - min(patient_trial) + 1;
int<lower = 1> n_train_patients = train_patients_end - train_patients_pos + 1;
array[n_train_patients] int<lower = 0> n_train_patient_visits = n_patient_visits[train_patients_pos:train_patients_end];
int<lower = 1> n_total_train_visits = sum(n_train_patient_visits);
int<lower = 1> n_total_train_visits_m1 = n_total_train_visits - n_train_patients;

array[n_total_train_visits] int train_patient_visits = get_int_sub_array(t_patient_visits, patient_visit_pos, train_patients_pos, train_patients_end);

array[n_train_patients + 1] int<lower = 1> train_patient_visit_pos = create_pos(n_train_patient_visits);
array[n_train_patients + 1] int<lower = 1> train_patient_visit_m1_pos = create_pos(n_train_patient_visits, -1);
array[n_train_patients + 1] int<lower = 1> train_forecast_visits_pos = create_pos(n_patient_forecast_visits[train_patients_pos:train_patients_end]);

int<lower = 1> n_total_train_forecast_visits = get_pos_total_size(train_forecast_visits_pos);

array[n_total_train_visits] int train_obs_recist = get_int_sub_array(recist, patient_visit_pos, train_patients_pos, train_patients_end);

int<lower = 0, upper = n_patients> n_train_right_censored_patients = sum(right_censored[train_patients_pos:train_patients_end]);
int<lower = 0, upper = n_patients> n_train_right_uncensored_patients = n_patients - n_train_right_censored_patients;
array[n_train_right_uncensored_patients] int<lower = 1, upper = n_patients> train_right_uncensored_patients;

array[n_trials + 1] int train_trial_right_censored_pos, train_trial_right_uncensored_pos;

array[n_trials + 1] int train_trial_patient_pos = resize_pos(trial_patient_pos, train_patients_pos, train_patients_end);
print("train_trial_patient_pos = ", train_trial_patient_pos);

array[n_train_patients] int<lower = 1> train_patient_trial = patient_trial[train_patients_pos:train_patients_end];

{
  int right_uncensored_idx = 1;
  array[n_trials] int n_train_trial_censored = zeros_int_array(n_trials), n_train_trial_uncensored = zeros_int_array(n_trials);

  for (i in train_patients_pos:train_patients_end) {
    if (right_censored[i]) {
      n_train_trial_censored[patient_trial[i]] += 1;
    } else {
      train_right_uncensored_patients[right_uncensored_idx] = i;
      right_uncensored_idx += 1;
      n_train_trial_uncensored[patient_trial[i]] += 1;
    }
  }

  train_trial_right_censored_pos = create_pos(n_train_trial_censored);
  train_trial_right_uncensored_pos = create_pos(n_train_trial_uncensored);
}

array[n_patients] int<lower = 1> ub_pfs_p1;

for (i in 1:n_patients) {
  ub_pfs_p1[i] = pfs[i] + interval_censored[i] + 1;
}

real log_lod = log(0.1);

// Define RECIST categories as integers
int CR = 1;  // Complete Response
int PR = 2;  // Partial Response
int SD = 3;  // Stable Disease
int PD = 4;  // Progressive Disease

// Define non-target status categories
int NT_CR = 1;
int NT_STABLE = 2;  // Non-CR/Non-PD
int NT_PD = 3;

matrix[n_train_patients, n_covar] Q_covar_design_matrix = 
  qr_thin_Q(covar_design_matrix[train_patients_pos:train_patients_end]) * sqrt(n_train_patients - 1);
matrix[n_covar, n_covar] R_covar_design_matrix = qr_thin_R(covar_design_matrix[train_patients_pos:train_patients_end]) / sqrt(n_train_patients - 1);
matrix[n_covar, n_covar] R_inv_covar_design_matrix = inverse(R_covar_design_matrix);

real log_abs_det_R_covar_design_matrix = log(abs(determinant(R_covar_design_matrix)));