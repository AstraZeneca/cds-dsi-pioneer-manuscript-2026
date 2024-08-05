
/** How many assessments for each tumor were pre-screening assessments (t <= 0).
 *
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @return Number of pre-screening observed assessments per tumor
 */
array[] int calc_n_screening_t(array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure) {
  int n_patients = size(n_patient_tumors);
  array[sum(n_patient_tumors)] int n_screening_t = rep_array(0, sum(n_patient_tumors));
  
  int tumor_pos = 1;
  int t_measure_pos = 1;
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1;
      
      for (tp in t_measure_pos:t_measure_end) {
        if (t_measure[tp] <= 0) {
          n_screening_t[tumor_pos + j - 1] += 1;
        }
      }
      
      t_measure_pos = t_measure_end + 1;
    }
    
    tumor_pos = tumor_end + 1;
  }
  
  return n_screening_t;
}
