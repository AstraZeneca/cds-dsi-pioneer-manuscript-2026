// Back-transformed coefficients from QR space to original space
#include "modules/tr/generated_quantities.stan"
#include "modules/frac/generated_quantities.stan"
#include "modules/init/generated_quantities.stan"

// RECIST prediction accuracy metrics - trial level
array[n_trials] int<lower=0> correct_recist_predictions = zeros_int_array(n_trials);
array[n_trials] int<lower=0> total_recist_predictions = zeros_int_array(n_trials);
array[n_trials] matrix<lower = 0>[PD, PD] recist_confusion_matrix; // rows = observed, cols = predicted
array[n_trials] real weighted_recist_accuracy_linear = zeros_array(n_trials);
array[n_trials] real weighted_recist_accuracy_quadratic = zeros_array(n_trials);
array[n_trials] int<lower=0> correct_recist_response_class = zeros_int_array(n_trials);
array[n_trials] int<lower=0> correct_recist_disease_control = zeros_int_array(n_trials);

// Per-category metrics - trial level
array[n_trials] vector[PD] recist_category_sensitivity; // true positive rate per category
array[n_trials] vector[PD] recist_category_precision;   // positive predictive value per category
array[n_trials] vector[PD] recist_category_counts;      // number of observations per category

for (s in 1:n_trials) {
  recist_confusion_matrix[s] = rep_matrix(0, PD, PD);
  recist_category_sensitivity[s] = zeros_vector(PD);
  recist_category_precision[s] = zeros_vector(PD);
  recist_category_counts[s] = zeros_vector(PD);
}

for (i in 1:n_patients) {
  int trial_idx = patient_trial[i];

  int visit_start, visit_end;
  (visit_start, visit_end) = get_pos(patient_visit_pos, i);
  
  for (t in (visit_start + n_patient_screening_visits[i]):visit_end) {
    if (recist[t] <= PD) {
      if (rep_recist[t] > PD) {
        reject(i, ": rep_recist = ", rep_recist[visit_start:visit_end], 
                  ", recist = ", recist[visit_start:visit_end]);
      }
      
      // Update confusion matrix inline
      recist_confusion_matrix[trial_idx][recist[t], rep_recist[t]] += 1;
      
      // Update all other metrics using the function
      (correct_recist_predictions[trial_idx], recist_category_counts[trial_idx],
       weighted_recist_accuracy_linear[trial_idx], weighted_recist_accuracy_quadratic[trial_idx],
       correct_recist_response_class[trial_idx], correct_recist_disease_control[trial_idx]) = update_recist_metrics(
        recist[t], rep_recist[t],
        correct_recist_predictions[trial_idx], recist_category_counts[trial_idx],
        weighted_recist_accuracy_linear[trial_idx], weighted_recist_accuracy_quadratic[trial_idx],
        correct_recist_response_class[trial_idx], correct_recist_disease_control[trial_idx]
      );
      
      total_recist_predictions[trial_idx] += 1;
    } 
  }
}

// Calculate final metrics - trial level
array[n_trials] real recist_accuracy;
array[n_trials] real recist_response_accuracy;
array[n_trials] real recist_disease_control_accuracy;
array[n_trials] real recist_response_sensitivity;
array[n_trials] real recist_response_specificity;
array[n_trials] real recist_progression_sensitivity;
array[n_trials] real recist_progression_specificity;

for (s in 1:n_trials) {
  (recist_accuracy[s], recist_response_accuracy[s], recist_disease_control_accuracy[s],
   weighted_recist_accuracy_linear[s], weighted_recist_accuracy_quadratic[s],
   recist_category_sensitivity[s], recist_category_precision[s],
   recist_response_sensitivity[s], recist_response_specificity[s],
   recist_progression_sensitivity[s], recist_progression_specificity[s]) = calculate_recist_summary_metrics(
    correct_recist_predictions[s], total_recist_predictions[s],
    recist_confusion_matrix[s], recist_category_counts[s],
    weighted_recist_accuracy_linear[s], weighted_recist_accuracy_quadratic[s],
    correct_recist_response_class[s], correct_recist_disease_control[s]
  );
}