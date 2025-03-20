
// Variables for handling separate baseline and proportional hazards
int n_tumor_separate_trials = separate_trial_tumor_gp ? n_trials : 1;

print("n_tumor_separate_trials = ", n_tumor_separate_trials);

vector[sum(n_patient_visits)] log_sum_tumor_size = log(sum_tumor_size); // cm

array[n_trials] int<lower = 0> n_trial_non_measured_tumors = zeros_int_array(n_trials);
array[n_patients] int<lower = 0> n_patient_non_measured_tumors = zeros_int_array(n_patients);

array[n_trials + 1] int<lower = 1> trial_non_measured_tumors_pos;
trial_non_measured_tumors_pos[1] = 1;
array[n_trials + 1] int<lower = 1> trial_measured_tumors_pos;
trial_measured_tumors_pos[1] = 1;

array[n_patients + 1] int<lower = 1> patient_non_measured_tumors_pos;
patient_non_measured_tumors_pos[1] = 1;
array[n_patients + 1] int<lower = 1> patient_measured_tumors_pos;
patient_measured_tumors_pos[1] = 1;

int min_visit = min(t_patient_visits), max_visit = max(t_patient_visits), min_missing_visit = min(t_patient_missing_visits), max_missing_visit = max(t_patient_missing_visits);
int<lower = 1> max_idx = max_visit - min_visit + 1, max_missing_visit_idx = max_missing_visit - min_missing_visit + 1;
int<lower = 1> n_pop_unique_visits = num_unique(t_patient_visits);
int<lower = 1> n_pop_unique_missing_visits = num_unique(t_patient_missing_visits);
array[n_pop_unique_visits] int pop_unique_visits = unique(t_patient_visits);
array[n_pop_unique_missing_visits] int pop_unique_missing_visits = unique(t_patient_missing_visits);
array[n_pop_unique_visits] int<lower = 1, upper = max_idx> pop_unique_visits_idx; 
array[n_pop_unique_missing_visits] int<lower = 1, upper = max_missing_visit_idx> pop_unique_missing_visits_idx; 
array[max_idx] int<lower = 0, upper = n_pop_unique_visits> pop_unique_visits_idx_dict; 
array[max_missing_visit_idx] int<lower = 0, upper = n_pop_unique_missing_visits> pop_unique_missing_visits_idx_dict; 

(pop_unique_visits_idx, pop_unique_visits_idx_dict) = id2idx(pop_unique_visits);
(pop_unique_missing_visits_idx, pop_unique_missing_visits_idx_dict) = id2idx(pop_unique_missing_visits);

for (s in 1:n_trials) {
  int curr_patient_pos, curr_patient_end; 
  (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
  
  for (i in curr_patient_pos:curr_patient_end) {
    int curr_patient_visits_pos, curr_patient_visits_end;
    (curr_patient_visits_pos, curr_patient_visits_end) = get_pos(patient_visit_pos, i);
    
    for (m in curr_patient_visits_pos:curr_patient_visits_end) {
      if (sum_tumor_size[m] <= 0) {
        n_patient_non_measured_tumors[i] += 1;
      }
    }
    
    patient_non_measured_tumors_pos[i + 1] = sum(n_patient_non_measured_tumors[:i]) + 1;
    patient_measured_tumors_pos[i + 1] = sum(n_patient_visits[:i]) - sum(n_patient_non_measured_tumors[:i]) + 1;
  }
  
  n_trial_non_measured_tumors[s] += sum(n_patient_non_measured_tumors[curr_patient_pos:curr_patient_end]);
  trial_non_measured_tumors_pos[s + 1] = sum(n_trial_non_measured_tumors[:s]) + 1;
  trial_measured_tumors_pos[s + 1] = sum(n_patient_visits[:curr_patient_end]) - sum(n_trial_non_measured_tumors[:s]) + 1;
}

array[sum(n_trial_non_measured_tumors)] int<lower = 1, upper = sum(n_patient_visits)> non_measured_tumors;
array[sum(n_patient_visits) - sum(n_trial_non_measured_tumors)] int<lower = 1, upper = sum(n_patient_visits)> measured_tumors;
array[sum(n_trial_non_measured_tumors)] int<lower = 1, upper = sum(n_patient_visits)> global_non_measured_tumors_idx;
array[sum(n_patient_visits) - sum(n_trial_non_measured_tumors)] int<lower = 1, upper = sum(n_patient_visits)> global_measured_tumors_idx;

for (s in 1:n_trials) {
  int curr_patient_pos, curr_patient_end;
  int curr_non_measure_pos, curr_non_measure_end;
  int curr_measure_pos, curr_measure_end;
  
  (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
  (curr_non_measure_pos, curr_non_measure_end) = get_pos(trial_non_measured_tumors_pos, s);
  (curr_measure_pos, curr_measure_end) = get_pos(trial_measured_tumors_pos, s);
  
  int trial_non_measured_offset = 0, trial_measured_offset = 0;
  
  for (i in curr_patient_pos:curr_patient_end) {
    int curr_patient_visits_pos, curr_patient_visits_end;
    int curr_patient_non_measure_pos, curr_patient_non_measure_end;
    int curr_patient_measure_pos, curr_patient_measure_end;
    
    (curr_patient_visits_pos, curr_patient_visits_end) = get_pos(patient_visit_pos, i);
    (curr_patient_non_measure_pos, curr_patient_non_measure_end) = get_pos(patient_non_measured_tumors_pos, i);
    (curr_patient_measure_pos, curr_patient_measure_end) = get_pos(patient_measured_tumors_pos, i);
    
    for (m in curr_patient_visits_pos:curr_patient_visits_end) {
      if (sum_tumor_size[m] <= 0) {
        non_measured_tumors[curr_non_measure_pos + trial_non_measured_offset] = m;
        global_non_measured_tumors_idx[curr_non_measure_pos + trial_non_measured_offset] = global_t_visit_idx[m];
        trial_non_measured_offset += 1;
      } else {
        measured_tumors[curr_measure_pos + trial_measured_offset] = m;
        global_measured_tumors_idx[curr_measure_pos + trial_measured_offset] = global_t_visit_idx[m];
        trial_measured_offset += 1;
      }
    }
  }
}

// real mean_log_sum_tumor_size = mean(log_sum_tumor_size[measured_tumors]);
// real<lower = 0> sd_log_sum_tumor_size = sd(log_sum_tumor_size[measured_tumors]);
// vector[sum(n_patient_visits)] std_log_sum_tumor_size = (log_sum_tumor_size - mean_log_sum_tumor_size) / sd_log_sum_tumor_size;
