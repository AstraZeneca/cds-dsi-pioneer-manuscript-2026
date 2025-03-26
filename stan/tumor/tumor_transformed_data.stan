
// Variables for handling separate baseline and proportional hazards
int n_tumor_separate_trials = separate_trial_tumor_gp ? n_trials : 1;

print("n_tumor_separate_trials = ", n_tumor_separate_trials);

vector[sum(n_patient_visits)] log_sum_tumor_size = log(sum_tumor_size); // cm
vector[sum(n_patient_visits) - sum(n_patient_screening_visits)] post_treat_sld;
array[n_patients] int<lower = 0> n_patient_post_treat_visits = zeros_int_array(n_patients);

// Population indices ////

int<lower = 0> n_pop_unique_visits = num_unique(t_patient_visits);
array[n_pop_unique_visits] int pop_unique_visits = unique(t_patient_visits);
array[n_pop_unique_visits] int<lower = 1> pop_unique_visits_idx = id2idx(pop_unique_visits); 

// Trial indices ////

array[n_trials] int<lower = 0> n_trial_unique_visits = num_unique(t_patient_visits, trial_visit_pos);
array[sum(n_trial_unique_visits)] int trial_unique_visits;
array[n_trials + 1] int<lower = 1> trial_unique_visits_pos;
(trial_unique_visits, trial_unique_visits_pos) = unique_by_pos(t_patient_visits, trial_visit_pos);
array[sum(n_trial_unique_visits)] int<lower = 1> trial_unique_visits_idx = id2idx(trial_unique_visits, trial_unique_visits_pos);
array[sum(n_trial_unique_visits)] int<lower = 1, upper = n_pop_unique_visits> trial2pop_unique_visit_idx = 
  get_level2level_idx(pop_unique_visits, trial_unique_visits, trial_unique_visits_pos); 
 
// Patient indicies ////

// This part kind of repeats what happens in base_transformed_data.stan but here we're excluding screening
array[n_patients] int<lower = 0> n_patient_unique_visits = num_unique(t_patient_visits, patient_visit_pos);
array[sum(n_patient_unique_visits)] int patient_unique_visits;
array[n_patients + 1] int<lower = 1> patient_unique_visits_pos;
(patient_unique_visits, patient_unique_visits_pos) = unique_by_pos(t_patient_visits, patient_visit_pos);
array[sum(n_patient_unique_visits)] int<lower = 1> patient_unique_visits_idx = id2idx(patient_unique_visits, patient_unique_visits_pos);
array[sum(n_patient_unique_visits)] int<lower = 1, upper = max(n_trial_unique_visits)> patient2trial_unique_visit_idx = 
  get_level2level_idx(trial_unique_visits, trial_unique_visits_pos, patient_unique_visits, patient_unique_visits_pos, trial_patient_pos); 
  
array[n_patients] int<lower = 1> patient_max_unique_visits_idx = get_max_pos(patient_unique_visits_idx, patient_unique_visits_pos);
array[n_patients + 1] int<lower = 1> patient_full_visits_pos = create_pos(patient_max_unique_visits_idx);
  
// Measured and non-measured SLD ////

array[n_patients] int<lower = 0> n_patient_non_measured_tumor_visits = zeros_int_array(n_patients);

array[n_patients + 1] int<lower = 1> patient_non_measured_tumor_visits_pos;
patient_non_measured_tumor_visits_pos[1] = 1;
array[n_patients + 1] int<lower = 1> patient_measured_tumor_visits_pos;
patient_measured_tumor_visits_pos[1] = 1;

{
  int post_treat_pos = 1;
  
  for (s in 1:n_trials) {
    int curr_patient_pos, curr_patient_end; 
    (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
    
    for (i in curr_patient_pos:curr_patient_end) {
      int curr_patient_visits_pos, curr_patient_visits_end;
      (curr_patient_visits_pos, curr_patient_visits_end) = get_pos(patient_visit_pos, i);
      
      for (m in curr_patient_visits_pos:curr_patient_visits_end) {
        if (t_patient_visits[m] > 0) {
          if (sum_tumor_size[m] <= 0) {
            n_patient_non_measured_tumor_visits[i] += 1;
          }
          
          post_treat_sld[post_treat_pos] = sum_tumor_size[m];
          post_treat_pos += 1;
        }
      }
      
      patient_non_measured_tumor_visits_pos[i + 1] = sum(n_patient_non_measured_tumor_visits[:i]) + 1;
      n_patient_post_treat_visits[i] = n_patient_visits[i] - n_patient_screening_visits[i];
      patient_measured_tumor_visits_pos[i + 1] = sum(n_patient_post_treat_visits[:i]) - sum(n_patient_non_measured_tumor_visits[:i]) + 1;
    }
  }
}

array[n_patients + 1] int<lower = 1> post_treat_visits_pos = create_pos(n_patient_post_treat_visits);

array[sum(n_patient_non_measured_tumor_visits)] int<lower = 1, upper = sum(n_patient_visits)> non_measured_tumor_visits;
array[sum(n_patient_non_measured_tumor_visits)] int<lower = 1> non_measured2patient_visits_idx;

array[sum(n_patient_post_treat_visits) - sum(n_patient_non_measured_tumor_visits)] int<lower = 1, upper = sum(n_patient_post_treat_visits)> measured_tumor_visits;
array[sum(n_patient_post_treat_visits) - sum(n_patient_non_measured_tumor_visits)] int<lower = 1> measured2patient_visits_idx;

for (s in 1:n_trials) {
  int curr_patient_pos, curr_patient_end;
  (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
  
  for (i in curr_patient_pos:curr_patient_end) {
    int curr_patient_visits_pos, curr_patient_visits_end;
    int curr_patient_non_measure_pos, curr_patient_non_measure_end;
    int curr_patient_measure_pos, curr_patient_measure_end;
    
    (curr_patient_visits_pos, curr_patient_visits_end) = get_pos(patient_visit_pos, i);
    (curr_patient_non_measure_pos, curr_patient_non_measure_end) = get_pos(patient_non_measured_tumor_visits_pos, i);
    (curr_patient_measure_pos, curr_patient_measure_end) = get_pos(patient_measured_tumor_visits_pos, i);
    
    int patient_non_measured_offset = 0, patient_measured_offset = 0, curr_patient_visit = 1, post_treat_offset = 0;
    
    for (m in curr_patient_visits_pos:curr_patient_visits_end) {
      if (t_patient_visits[m] > 0) {
        if (sum_tumor_size[m] <= 0) {
          non_measured_tumor_visits[curr_patient_non_measure_pos + patient_non_measured_offset] = t_patient_visits[m];
          non_measured2patient_visits_idx[curr_patient_non_measure_pos + patient_non_measured_offset] = curr_patient_visit;
          patient_non_measured_offset += 1;
          curr_patient_visit += 1;
        } else {
          measured_tumor_visits[curr_patient_measure_pos + patient_measured_offset] = t_patient_visits[m];
          measured2patient_visits_idx[curr_patient_measure_pos + patient_measured_offset] = curr_patient_visit;
          patient_measured_offset += 1;
          curr_patient_visit += 1;
        }
      } else {
        post_treat_offset += 1;
      }
    }
  }
}
