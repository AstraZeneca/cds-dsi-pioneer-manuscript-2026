/**
 * Count the number of positive (>0) values in an integer array.
 *
 * @param arr Array of integers
 * @return Number of elements in arr that are >0
 */
int count_positive(array[] int arr) {
  int n = size(arr);
  int count = 0;
  for (i in 1:n) {
    if (arr[i] > 0) {
      count += 1;
    }
  }
  return count;
}
/** Calculate the last observed measure for each patient. 
 *
 * @param t_measure The week each assessment was done.
 * @param n_measures The number of assessments per tumor.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @return Get the week of last observation per patient.
 */
array[] int get_max_t(array[] int t_measure, array[] int n_measures, array[] int n_patient_tumors) {
  int n_patients = size(n_patient_tumors);
  int t_pos = 1;
  int tumor_pos = 1;
  array[n_patients] int max_t;

  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int n_meas_sum = 0;
    for (j in tumor_pos:tumor_end) n_meas_sum += n_measures[j];
    if (n_patient_tumors[i] == 0 || n_meas_sum == 0) {
      // No tumors or no measures for this patient, set to 0 or a sentinel value
      max_t[i] = 0;
      t_pos += n_meas_sum; // still advance t_pos for consistency
    } else {
      int t_end = t_pos + n_meas_sum - 1;
      max_t[i] = max(t_measure[t_pos:t_end]);
      t_pos = t_end + 1;
    }
    tumor_pos = tumor_end + 1;
  }
  return max_t;
}

tuple(real, real, real) summarize_matrix_eigenvalues(matrix m) {
  int n = rows(m);
  
  vector[n] eigenvalues = eigenvalues_sym(m);
  real min_eigenvalue = min(eigenvalues);
  real max_eigenvalue = max(eigenvalues);
  real condition_number = max_eigenvalue / min_eigenvalue;
  
  return (min_eigenvalue, max_eigenvalue, condition_number);
}

matrix diag_matrix(real x, int n) {
  return diag_matrix(rep_vector(x, n));
}

/** Missing measure is defined as one that lies between a _tumor's_ first assessment to the _patient's_ last assessment. Basically,
 * we're counting how many intervals (weeks) we don't have observed assessments of tumor size, for each tumor.
 *
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @return Number of missing assessments per tumor
 */
array[] int calculate_n_missing_measures(array[] int n_measures, array[] int t_measure, array[] int n_patient_tumors) {
  int tumor_pos = 1;
  int t_measure_pos = 1;
  int n_patients = size(n_patient_tumors);
  int n_tumors = size(n_measures);
  array[n_tumors] int n_missing_measures = zeros_int_array(n_tumors);
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int first_tumor_t_measure_pos = t_measure_pos;
    int last_tumor_t_measure_end = t_measure_pos + sum(n_measures[tumor_pos:tumor_end]) - 1;
    
    int max_patient_t = max(t_measure[first_tumor_t_measure_pos:last_tumor_t_measure_end]);
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
      
      int min_tumor_t = min(t_measure[t_measure_pos:t_measure_end]);
      int full_patient_measure_width = max_patient_t - min_tumor_t + 1;
      
      n_missing_measures[tumor_pos] = full_patient_measure_width - n_measures[tumor_pos];
      
      tumor_pos += 1;  
      t_measure_pos = t_measure_end + 1;
    }
    
  }
  
  return n_missing_measures;
}

tuple(array[] int, array[] int) calculate_n_missing_visits(array[] int visit_pos, array[] int visits) {
  return calculate_n_missing_visits(visit_pos, visits, 0, 0);
}

tuple(array[] int, array[] int) calculate_n_missing_visits(array[] int visit_pos, array[] int visits, int make_unique, int post_treatment) {
  int n = size(visit_pos) - 1;
  array[n] int n_missing = zeros_int_array(n);
  
  for (i in 1:n) {
    int curr_visit_pos, curr_visit_end;
    (curr_visit_pos, curr_visit_end) = get_pos(visit_pos, i);
   
    if (post_treatment) { 
      while (curr_visit_pos <= curr_visit_end && visits[curr_visit_pos] <= 0) {
        curr_visit_pos += 1;
      }
    }
    
    int max_t_width = max(visits[curr_visit_pos:curr_visit_end]) - min(visits[curr_visit_pos:curr_visit_end]) + 1;
    int n_visits = make_unique ? num_unique(visits[curr_visit_pos:curr_visit_end]) : (curr_visit_end - curr_visit_pos + 1);
    
    n_missing[i] = max_t_width - n_visits;
  }
  
  return (n_missing, create_pos(n_missing));
}


int calculate_n_missing_visits(array[] int unique_visits, int n_full) {
  return calculate_n_missing_visits(unique_visits, {1, num_elements(unique_visits) + 1}, n_full)[1];
}

array[] int calculate_n_missing_visits(array[] int unique_visits, array[] int unique_visits_pos, int n_full) {
  int n = size(unique_visits_pos) - 1;
  array[n] int n_unique_missing_visits = zeros_int_array(n);
  
  for (p in 1:n) {
    n_unique_missing_visits[p] = n_full - get_pos_size(unique_visits_pos, p);
  }
  
  return n_unique_missing_visits;
}

tuple(array[] int, array[] int) get_missing_visits(array[] int visit_pos, array[] int visits, array[] int missing_visit_pos) {
  return get_missing_visits(visit_pos, visits, missing_visit_pos, 0);
}

tuple(array[] int, array[] int) get_missing_visits(array[] int visit_pos, array[] int visits, array[] int missing_visit_pos, int post_treatment) {
  int n = size(visit_pos) - 1;
  array[n] int missing_size = get_pos_size(missing_visit_pos);
  int n_total_missing = sum(missing_size);
  array[n_total_missing] int missing_visits, missing_visits_idx;
  
  for (i in 1:n) {
    int missing_pos, missing_end;
    (missing_pos, missing_end) = get_pos(missing_visit_pos, i);
   
    int n_curr_visits = get_pos_size(visit_pos, i), visit_offset = 0; 
    array[n_curr_visits] int sorted_visits = sort_asc(get_int_sub_array(visits, visit_pos, i)); 
    
    if (post_treatment) {
      while (visit_offset < n_curr_visits && sorted_visits[visit_offset + 1] <= 0) {
        visit_offset += 1;
      }
    }
    
    int min_visit = sorted_visits[1 + visit_offset], max_visit = sorted_visits[n_curr_visits];
    
    int visit_count = 0;
    int next_visit = min_visit;
    
    for (m in missing_pos:missing_end) {
      while (next_visit < max_visit && sorted_visits[visit_count + 1 + visit_offset] == next_visit) {
        visit_count += 1;
        next_visit += 1;
      }
      
      missing_visits[m] = next_visit;
      missing_visits_idx[m] = next_visit - min_visit + 1;
      next_visit += 1;
    }
  }
  
  return (missing_visits, missing_visits_idx);
}

array[] int get_missing_visits(array[] int unique_visits, int n_full) {
  return get_missing_visits(unique_visits, {1, num_elements(unique_visits) + 1}, n_full);
}

array[] int get_missing_visits(array[] int unique_visits, array[] int unique_visits_pos, int n_full) {
  int n = size(unique_visits_pos) - 1;
  array[n] int n_unique_visits = get_pos_size(unique_visits_pos);
  array[n] int n_unique_missing_visits = zeros_int_array(n);
  
  for (p in 1:n) {
    n_unique_missing_visits[p] = n_full - get_pos_size(unique_visits_pos, p);
  }
  
  array[sum(n_unique_missing_visits)] int unique_missing_visits;
  int curr_missing_idx = 1;
 
  for (p in 1:n) { 
    array[n_unique_visits[p]] int curr_unique_visits = sort_asc(get_int_sub_array(unique_visits, unique_visits_pos, p));
    int curr_unique_visit_idx = 1;
    
    for (q in 1:n_full) {
      if (curr_unique_visit_idx > n_unique_visits[p] || q < curr_unique_visits[curr_unique_visit_idx]) {
        unique_missing_visits[curr_missing_idx] = q;
        curr_missing_idx += 1;
      } else {
        curr_unique_visit_idx += 1;
      } 
    }
  }
  
  return unique_missing_visits;
}

/** Return the actual t for which we don't have observed tumor size assessments.
 *
 * @param n_measures The number of assessments per tumor.
 * @param n_missing_measures Number of missing assessments per tumor
 * @param t_measure The week each assessment was done.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @return The actual weeks in which assessments are unobserved
 */
array[] int calculate_t_missing_measure(
  array[] int n_measures, array[] int n_missing_measures, array[] int t_measure, array[] int n_patient_tumors 
) { 
  int n_patients = size(n_patient_tumors); 
  int tumor_pos = 1;
  int t_measure_pos = 1;
  int t_missing_measure_pos = 1;
  int n_tumors = size(n_measures);
  array[sum(n_missing_measures)] int t_missing_measure;
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int first_tumor_t_measure_pos = t_measure_pos;
    int last_tumor_t_measure_end = t_measure_pos + sum(n_measures[tumor_pos:tumor_end]) - 1;
    
    int max_patient_t = max(t_measure[first_tumor_t_measure_pos:last_tumor_t_measure_end]);
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
      
      int min_tumor_t = min(t_measure[t_measure_pos:t_measure_end]);
      int measures_checked = 0;
      
      for (k in min_tumor_t:max_patient_t) {
        if (measures_checked >= n_measures[tumor_pos] || t_measure[t_measure_pos] > k) {
          t_missing_measure[t_missing_measure_pos] = k;
          t_missing_measure_pos += 1;
        } else {
          t_measure_pos += 1;
          measures_checked += 1;
        }
      }
      
      tumor_pos += 1;
    }
  }
  
  return t_missing_measure;
}

/** Get number of values in x that are less than or equal to y
 */
int num_leq(array[] int x, int y) {
  int n = 0;
  array[size(x)] int sorted_x = sort_asc(x);
  
  for (i in 1:size(x)) {
    if (sorted_x[i] <= y) {
      n += 1;
    } else {
      break;
    }
  }
  
  return n;
}

/**
 * Return the indices of elements > 0 (or <= 0 if inverse=1) in an indicator array.
 * @param mask Array of integers (0/1 or any integer)
 * @param inverse If 1, return indices where mask <= 0; if 0 (default), return indices where mask > 0
 * @return Array of indices (1-based) where mask > 0 (or <= 0 if inverse=1)
 */
array[] int which(array[] int mask, int inverse) {
  int n = size(mask);
  int count = 0;
  for (i in 1:n) {
    if ((inverse == 0 && mask[i] > 0) || (inverse == 1 && mask[i] <= 0)) count += 1;
  }
  array[count] int idx;
  int pos = 1;
  for (i in 1:n) {
    if ((inverse == 0 && mask[i] > 0) || (inverse == 1 && mask[i] <= 0)) {
      idx[pos] = i;
      pos += 1;
    }
  }
  return idx;
}

/**
 * Overload: Return the indices of elements equal to 1 in a binary indicator array.
 * @param mask Binary array (0/1)
 * @return Array of indices (1-based) where mask == 1
 */
array[] int which(array[] int mask) {
  return which(mask, 0);
}

/** Identify which elements in a binary array are 0 and which are 1.
 * @param mask Binary array
 * @return tuple(indices of 0 elements, indices of 1 elements)
 */
tuple(array[] int, array[] int) get_mask_idx(array[] int mask) {
  return (which(mask, 1), which(mask, 0));
}

/** Repeat each value a specific number of times.
 */
array[] int rep_each(array[] int to_repeat, int repeats) {
  int n = size(to_repeat);
  array[n * repeats] int repeated;
  
  int pos = 1;
  
  for (i in 1:n) {
    int end = pos + repeats - 1;
    repeated[pos:end] = rep_array(to_repeat[i], repeats);
    pos = end + 1;
  }
  
  return(repeated);
}

real months_to_weeks(int mon) {
  return mon * 365.25 / (7 * 12);
}

int calendar_date_to_study_date(int first_calendar_date, int calendar_date) {
  return calendar_date - first_calendar_date + 1;
}

array[] int calendar_date_to_study_date(array[] int first_calendar_date, array[] int calendar_date) {
  int n = size(first_calendar_date);
  array[n] int study_date;
  
  for (i in 1:n) {
    study_date[i] = calendar_date_to_study_date(first_calendar_date[i], calendar_date[i]);
  }
  
  return study_date;
}

array[] int calendar_date_to_study_date(array[] int first_calendar_date, int calendar_date) {
  int n = size(first_calendar_date);
  array[n] int study_date;
  
  for (i in 1:n) {
    study_date[i] = calendar_date_to_study_date(first_calendar_date[i], calendar_date);
  }
  
  return study_date;
}

int study_date_to_calendar_date(int first_calendar_date, int study_date) {
  return first_calendar_date + study_date - 1;
}

array[] int study_date_to_calendar_date(array[] int first_calendar_date, array[] int study_date) {
  int n = size(first_calendar_date);
  array[n] int calendar_date;
  
  for (i in 1:n) {
    calendar_date[i] = study_date_to_calendar_date(first_calendar_date[i], study_date[i]);
  }
  
  return calendar_date;
}

int num_unique(array[] int x) {
  return num_unique(x, 1);
}

int num_unique(array[] int x, int post_treatment) {
  int n = size(x);
  int count = 0, last = min(x) - 1;
  array[n] int sorted_x = sort_asc(x);
  
  for (i in 1:n) {
    if ((!post_treatment || sorted_x[i] > 0) && sorted_x[i] > last) {
      count += 1;
      last = sorted_x[i];
    }
  }
  
  return count;
}


array[] int num_unique(array[] int x, array[] int pos) {
  return num_unique(x, pos, 1);
}

array[] int num_unique(array[] int x, array[] int pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int count = zeros_int_array(n);
  
  for (i in 1:n) {
    count[i] = num_unique(get_int_sub_array(x, pos, i), post_treatment);
  }
  
  return count;
}

array[] int num_unique(array[] int x, array[] int pos, array[] int sub_pos) {
  return num_unique(x, pos, sub_pos, 1);
}
  
array[] int num_unique(array[] int x, array[] int pos, array[] int sub_pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int count = zeros_int_array(n);
  
  for (i in 1:n) {
    int i_pos, i_end;
    (i_pos, i_end) = get_pos(pos, i);
    
    count[i] = num_unique(get_int_sub_array(x, sub_pos, i_pos, i_end), post_treatment);
  }
  
  return count;
}

array[] int unique(array[] int x) {
  return unique(x, 1);
}
  
array[] int unique(array[] int x, int post_treatment) {
  int n = size(x);
  int count = 0, last = min(x) - 1;
  array[n] int sorted_x = sort_asc(x);
  int n_unique = num_unique(x, post_treatment);
  array[n_unique] int unique_x;
  
  for (i in 1:n) {
    
    if ((!post_treatment || sorted_x[i] > 0) && sorted_x[i] > last) {
      count += 1;
      last = sorted_x[i];
      unique_x[count] = last;
    }
  }
  
  return unique_x;
}

// If I call this unique the compiler complains about ambiguity which doesn't make sense ¯\_(ツ)_/¯
tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos) {
  return unique_by_pos(x, pos, 1);
}
  
tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int n_unique_x = num_unique(x, pos, post_treatment);
  array[n + 1] int unique_pos = create_pos(n_unique_x);
  array[sum(n_unique_x)] int unique_x;
  
  for (i in 1:n) {
    int n_i = get_pos_size(pos, i);
    array[n_i] int sorted_x_i = sort_asc(get_int_sub_array(x, pos, i));
    int count = 0, last = sorted_x_i[1] - 1;
    
    for (j in 1:n_i) {
      if ((!post_treatment || sorted_x_i[j] > 0) && sorted_x_i[j] > last) {
        last = sorted_x_i[j];
        unique_x[unique_pos[i] + count] = last; 
        count += 1;
      }
    }
  }
  
  return (unique_x, unique_pos);
}

tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos, array[] int sub_pos) {
  return unique_by_pos(x, pos, sub_pos, 1);
}

tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos, array[] int sub_pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int n_unique_x = num_unique(x, pos, sub_pos, post_treatment);
  array[n + 1] int unique_pos = create_pos(n_unique_x);
  array[sum(n_unique_x)] int unique_x;
  
  for (i in 1:n) {
    int i_pos, i_end;
    (i_pos, i_end) = get_pos(pos, i);
    
    
    int i_unique_pos, i_unique_end;
    (i_unique_pos, i_unique_end) = get_pos(unique_pos, i);
    
    unique_x[i_unique_pos:i_unique_end] = unique(get_int_sub_array(x, sub_pos, i_pos, i_end), post_treatment);
  }
  
  return (unique_x, unique_pos);
}

array[] int id2idx(array[] int id) {
  return id2idx(id, min(id));
}

array[] int id2idx(array[] int id, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int n_p = get_pos_size(pos);
  array[sum(n_p)] int idx;
  
  for (p in 1:n) {
    if (get_pos_size(pos, p) > 0) {
      int p_pos, p_end;
      (p_pos, p_end) = get_pos(pos, p);
    
      idx[p_pos:p_end] = id2idx(id[p_pos:p_end]);
    } 
  }
 
  return idx;
}

array[] int id2idx(array[] int id, int min_id) {
  int n = size(id);
  array[n] int idx;
  
  for (i in 1:n) {
    idx[i] = id[i] - min_id + 1;
  }
  
  return idx;
}

array[] int get_level2level_idx(array[] int hi_level, array[] int low_level) {
  int size_hi = size(hi_level), size_low = size(low_level);
 
  array[size_low] int idx = zeros_int_array(size_low);
  array[size_low] int sorted_low_level = sort_asc(low_level);
  int curr_low_idx = 1;
  
  for (h in 1:size_hi) {
    if (hi_level[h] == sorted_low_level[curr_low_idx]) {
      idx[curr_low_idx] = h;
      
      if (curr_low_idx == size_low) {
        break;
      }
      
      curr_low_idx += 1;
    }
  }
  
  // Check for any unmatched values
  for (i in 1:size_low) {
    if (idx[i] == 0) {
      reject("get_level2level_idx: Value ", low_level[i], " at position ", i, 
             " in low_level not found in hi_level. ",
             "low_level: ", low_level, ", hi_level: ", hi_level);
    }
  }
  
  return idx;
}

array[] int get_level2level_idx(array[] int hi_level, array[] int low_level, array[] int low_pos) {
  int n_low = size(low_pos) - 1;
  int size_low = size(low_level);
  array[size_low] int idx = zeros_int_array(size_low);
  
  for (l in 1:n_low) {
    if (get_pos_size(low_pos, l) > 0) {
      int pos, end;
      (pos, end) = get_pos(low_pos, l);
      
      idx[pos:end] = get_level2level_idx(hi_level, low_level[pos:end]);
    } 
  }
  
  return idx;
}

array[] int get_level2level_idx(array[] int hi_level, array[] int hi_pos, array[] int low_level, array[] int low_pos, array[] int low_hi_pos) {
  int n_low = size(low_pos) - 1, n_hi = size(hi_pos) - 1;
  int size_low = size(low_level);
  array[size_low] int idx = zeros_int_array(size_low);
 
  for (h in 1:n_hi) {
    int low_id_from, low_id_to;
    (low_id_from, low_id_to) = get_pos(low_hi_pos, h);
    int low_idx_start, low_idx_end;
    (low_idx_start, low_idx_end) = get_pos(low_pos, low_id_from, low_id_to);
    
    idx[low_idx_start:low_idx_end] = get_level2level_idx(
      get_int_sub_array(hi_level, hi_pos, h), get_int_sub_array(low_level, low_pos, low_id_from, low_id_to), create_pos(low_pos, low_id_from, low_id_to)
    );
  }
  
  return idx;
}

array[] int get_idx_dict(array[] int idx) {
  int max_idx = size(idx) > 0 ? max(idx) : 0;
  array[max_idx] int idx_dict = zeros_int_array(max_idx);
  for (i in 1:max_idx) {
    for (j in 1:size(idx)) {
      if (idx[j] == i) {
        idx_dict[i] = j;
        break;
      }
    }
  }
  return idx_dict;
}

tuple(array[] int, array[] int) get_idx_dict(array[] int idx, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int max_idx = get_max(idx, pos);
  array[n + 1] int dict_pos = create_pos(max_idx);
  array[sum(max_idx)] int idx_dict;

  for (p in 1:n) {
      int p_pos, p_end;
      (p_pos, p_end) = get_pos(dict_pos, p);
    if (max_idx[p] > 0) {
      idx_dict[p_pos:p_end] = get_idx_dict(get_int_sub_array(idx, pos, p));
    } 
  } 
  
  return (idx_dict, dict_pos);
}




void assert_equal(int x, int y) {
  if (x != y) {
    fatal_error("Equality assertion failed.");
  }
}

void assert_equal(real x, real y) {
  if (x != y) {
    fatal_error("Equality assertion failed.");
  }
}

void assert_greater_than_or_equal(int x, int y) {
  if (x < y) {
    fatal_error("Greater than or equal assertion failed.");
  }
}

void assert_greater(int x, int y) {
  if (x <= y) {
    fatal_error("Greater (strict) assertion failed.");
  }
}

void assert_less_or_equal(int x, int y) {
  if (x > y) {
    fatal_error("Less or equal assertion failed.");
  }
}

void assert_strict_ascending(array[] int x) {
  int n = size(x);
  for (i in 2:n) {
    if (x[i] <= x[i - 1]) {
      fatal_error("Array is not strictly ascending at position ", i, ": ", x[i], " <= ", x[i - 1]);
    }
  }
}

void assert_ascending(array[] int x) {
  int n = size(x);
  for (i in 2:n) {
    if (x[i] < x[i - 1]) {
      fatal_error("Array is not ascending at position ", i, ": ", x[i], " <= ", x[i - 1]);
    }
  }
}

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

/** Scale tumor sizes by the standard deviation of all tumors and demean.
 * 
 * @param tumor_size Observed tumor sizes
 * @return (Mean tumor size, Std deviation of tumor sizes, Standardized tumor sizes)
 */
tuple(real, real, vector) standardize_tumor_sizes(vector tumor_size) {
  real tumor_mean;
  real tumor_sd;
  
  tumor_mean = mean(tumor_size); 
  tumor_sd = sd(tumor_size); 
  
  return (tumor_mean, tumor_sd, (tumor_size - tumor_mean) / tumor_sd); 
}

/**
 * Find the first occurrence of n_succ consecutive elements from the 'what' array
 * within the 'all' array.
 *
 * @param all Array to search within
 * @param what Array containing the values to search for (can contain duplicates)
 * @param n_succ Number of consecutive elements to find
 * 
 * @return Starting index (1-based) of the first sequence of n_succ consecutive 
 *         elements where each element is contained in 'what'. Returns 0 if no 
 *         such sequence is found.
 *
 * @throws fatal_error if 'what' array is empty
 *
 * Example:
 *   all = [1, 2, 3, 2, 2, 4]
 *   what = [2, 3]
 *   n_succ = 3
 *   Returns: 2 (positions 2-4 contain [2, 3, 2], all from 'what')
 *
 * Note: Values from 'what' can be used multiple times in the sequence.
 *       For efficiency, 'what' is sorted internally for binary search.
 */
int find_first(array[] int all, array[] int what, int n_succ) {
  int n = size(all);
  int n_what = size(what);
  
  // Need at least n_succ elements in 'all' to find a match
  if (n < n_succ) {
    return 0;
  }
  
  // Need at least one element in 'what' to match against
  if (n_what == 0) {
    fatal_error("find_first: 'what' array cannot be empty");
  }
  
  // Sort 'what' array for faster searching
  array[n_what] int sorted_what = sort_asc(what);
  
  // Check each possible starting position in 'all'
  for (i in 1:(n - n_succ + 1)) {
    int matches = 0;
    
    // Check if we can find n_succ consecutive matches starting at position i
    for (j in 1:n_succ) {
      int found_match = 0;
      
      // Binary search in sorted_what for all[i + j - 1]
      int left = 1;
      int right = n_what;
      int search_val = all[i + j - 1];
      
      while (left <= right) {
        int mid = left + (right - left) %/% 2;
        if (sorted_what[mid] == search_val) {
          found_match = 1;
          break;
        } else if (sorted_what[mid] < search_val) {
          left = mid + 1;
        } else {
          right = mid - 1;
        }
      }
      
      if (found_match) {
        matches += 1;
      } else {
        break;  // No match at this position, move to next starting position
      }
    }
    
    // If we found n_succ consecutive matches, return the starting index
    if (matches == n_succ) {
      return i;
    }
  }
  
  // No sequence of n_succ consecutive matches found
  return 0;
}

int find_first(array[] int all, int what) {
  return find_first(all, { what }, 1);
}

/**
 * Create a compact array of observed patients and a mapping from original IDs to compact indices.
 * Patients are "observed" if their mask value is 1.
 *
 * @param observed_mask Array of 0/1 values indicating which patients are observed (sized n_patients)
 * @param n_observed Number of observed patients (pre-computed)
 * @return Tuple of (observed_patients array, patient_to_compact_idx mapping)
 *   - observed_patients[1:n_observed]: Original patient IDs for observed patients
 *   - patient_to_compact_idx[1:n_patients]: Maps original patient ID to compact index (0 if not observed)
 */
tuple(array[] int, array[] int) create_compact_patient_mapping(
  array[] int observed_mask,
  int n_observed
) {
  int n_patients = size(observed_mask);
  array[n_observed] int observed_patients;
  array[n_patients] int patient_to_compact_idx = zeros_int_array(n_patients);
  
  int obs_idx = 1;
  for (i in 1:n_patients) {
    if (observed_mask[i] == 1) {
      observed_patients[obs_idx] = i;
      patient_to_compact_idx[i] = obs_idx;
      obs_idx += 1;
    }
  }
  
  return (observed_patients, patient_to_compact_idx);
}

/**
 * Create trial/group position array for compact (observed-only) patients.
 * Counts how many observed patients belong to each trial/group.
 *
 * @param compact_patients Array of original patient IDs for observed patients
 * @param patient_group Array mapping original patient ID to group/trial ID (sized n_patients)
 * @param n_groups Number of groups/trials
 * @return Position array (sized n_groups + 1) created by create_pos
 */
array[] int create_compact_group_pos(
  array[] int compact_patients,
  array[] int patient_group,
  int n_groups
) {
  int n_compact = size(compact_patients);
  array[n_groups] int group_size = zeros_int_array(n_groups);
  
  for (obs_idx in 1:n_compact) {
    int i = compact_patients[obs_idx];
    int group_id = patient_group[i];
    group_size[group_id] += 1;
  }
  
  return create_pos(group_size);
}

/**
 * Remap a grouped array to use compact patient indices instead of original IDs.
 * This handles the common pattern of remapping cond_group or similar arrays.
 *
 * @param group Array of original patient IDs organized by group (sized total_entries)
 * @param group_pos Position array for the groups (sized n_groups + 1)
 * @param observed_mask Mask indicating which patients are observed (sized n_patients)
 * @param patient_to_compact_idx Mapping from original ID to compact index (sized n_patients)
 * @param n_groups Number of groups
 * @return Tuple of (compact_group array, compact_group_pos array)
 *   - compact_group: Remapped array with compact patient indices
 *   - compact_group_pos: New position array for compact groups
 */
tuple(array[] int, array[] int) remap_group_to_compact(
  array[] int group,
  array[] int group_pos,
  array[] int observed_mask,
  array[] int patient_to_compact_idx,
  int n_groups
) {
  int total_entries = size(group);
  
  // First pass: count how many observed patients are in each group
  array[n_groups] int compact_group_size = zeros_int_array(n_groups);
  
  for (g_idx in 1:total_entries) {
    int patient_id = group[g_idx];
    if (observed_mask[patient_id] == 1) {
      // Determine which group this entry belongs to
      int group_id = 1;
      while (group_id <= n_groups && g_idx >= group_pos[group_id + 1]) {
        group_id += 1;
      }
      compact_group_size[group_id] += 1;
    }
  }
  
  array[n_groups + 1] int compact_group_pos = create_pos(compact_group_size);
  int n_compact_entries = compact_group_pos[n_groups + 1] - 1;
  
  // Handle case when all groups are empty (no observed patients)
  if (n_compact_entries == 0) {
    return (rep_array(0, 0), compact_group_pos);  // Return empty array and position array
  }
  
  array[n_compact_entries] int compact_group;
  
  // Second pass: populate with compact patient indices
  array[n_groups] int group_fill_idx = compact_group_pos[1:n_groups];
  
  for (g_idx in 1:total_entries) {
    int patient_id = group[g_idx];
    if (observed_mask[patient_id] == 1) {
      // Determine which group this entry belongs to
      int group_id = 1;
      while (group_id <= n_groups && g_idx >= group_pos[group_id + 1]) {
        group_id += 1;
      }
      // Store the compact patient index
      compact_group[group_fill_idx[group_id]] = patient_to_compact_idx[patient_id];
      group_fill_idx[group_id] += 1;
    }
  }
  
  return (compact_group, compact_group_pos);
}

