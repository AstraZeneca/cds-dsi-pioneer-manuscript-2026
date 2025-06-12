// RECIST prediction accuracy metrics
int<lower=0> correct_recist_predictions = 0;
int<lower=0> total_recist_predictions = 0;
array[n_trials] matrix<lower = 0>[PD, PD] recist_confusion_matrix; // rows = observed, cols = predicted
real weighted_recist_accuracy_linear = 0;
real weighted_recist_accuracy_quadratic = 0;
int<lower=0> correct_recist_response_class = 0;
int<lower=0> correct_recist_disease_control = 0;

for (s in 1:n_trials) {
  recist_confusion_matrix[s] = rep_matrix(0, PD, PD);
}

// Per-category metrics
vector[PD] recist_category_sensitivity = zeros_vector(PD); // true positive rate per category
vector[PD] recist_category_precision = zeros_vector(PD);   // positive predictive value per category
vector[PD] recist_category_counts = zeros_vector(PD);      // number of observations per category

for (i in train_patients_pos:train_patients_end) {
  int train_idx = i - train_patients_pos + 1;

  int train_visit_start, train_visit_end;
  (train_visit_start, train_visit_end) = get_pos(train_patient_visit_pos, train_idx);
  
  for (t in train_visit_start:train_visit_end) {
    if (train_obs_recist[t] <= PD) {
      // Update all metrics using the function
      (correct_recist_predictions, recist_confusion_matrix[patient_trial[i]], recist_category_counts,
       weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
       correct_recist_response_class, correct_recist_disease_control) = update_recist_metrics(
        train_obs_recist[t], rep_recist[t],
        correct_recist_predictions, recist_confusion_matrix[patient_trial[i]], recist_category_counts,
        weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
        correct_recist_response_class, correct_recist_disease_control
      );
      
      total_recist_predictions += 1;
    } 
  }
}

// Calculate final metrics
real recist_accuracy;
real recist_response_accuracy;
real recist_disease_control_accuracy;
real recist_response_sensitivity;
real recist_response_specificity;
real recist_progression_sensitivity;
real recist_progression_specificity;

(recist_accuracy, recist_response_accuracy, recist_disease_control_accuracy,
 weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
 recist_category_sensitivity, recist_category_precision,
 recist_response_sensitivity, recist_response_specificity,
 recist_progression_sensitivity, recist_progression_specificity) = calculate_recist_summary_metrics(
  correct_recist_predictions, total_recist_predictions,
  recist_confusion_matrix[1], recist_category_counts, // BUG For now just using s = 1 for all metrics, need to break this down by trial 
  weighted_recist_accuracy_linear, weighted_recist_accuracy_quadratic,
  correct_recist_response_class, correct_recist_disease_control
);