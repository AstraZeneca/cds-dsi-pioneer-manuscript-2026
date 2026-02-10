/**
 * Position Array Utilities for Stan
 * 
 * This file provides utilities for working with ragged arrays and position-based indexing
 * in Stan models. It implements a position array pattern where data for multiple groups
 * (e.g., patients, trials) with varying sizes can be stored in flat arrays and accessed
 * using position indices.
 * 
 * Position arrays follow the convention where pos[i] gives the starting index for group i,
 * and pos[i+1]-1 gives the ending index. This allows efficient storage and retrieval of
 * variable-length data structures.
 * 
 * Example: For 3 groups with sizes [2, 0, 3], the position array would be [1, 3, 3, 6]
 * - Group 1: indices 1-2 (size 2)
 * - Group 2: indices 3-2 (size 0, empty group)
 * - Group 3: indices 3-5 (size 3)
 */

/**
 * Get maximum values for each group in a ragged array
 * 
 * @param id Array of values to find maxima in
 * @param pos Position array defining group boundaries
 * @return Array of maximum values for each group
 */
array[] int get_max(array[] int id, array[] int pos) {
  return get_max(id, pos, 0);
}

/**
 * Get maximum indices (relative to group start) for each group
 * 
 * @param id Array of values to find maxima in
 * @param pos Position array defining group boundaries
 * @return Array of indices where maximum values occur within each group
 */
array[] int get_max_idx(array[] int id, array[] int pos) {
  return get_max(id, pos, 1);
}

/**
 * Internal function to get maximum values or indices for each group
 * 
 * @param id Array of values to find maxima in
 * @param pos Position array defining group boundaries
 * @param of_idx If 0, return max values; if 1, return max indices
 * @return Array of maximum values or indices for each group
 */
array[] int get_max(array[] int id, array[] int pos, int of_idx) {
  int n = size(pos) - 1;  // Number of groups
  array[n] int p_max = zeros_int_array(n);
  
  for (i in 1:n) {
    int n_i = get_pos_size(pos, i);  // Size of group i
    
    if (n_i > 0) {
      array[n_i] int id_i = get_int_sub_array(id, pos, i);
      // Calculate maximum value or its index relative to group start
      p_max[i] = max(id_i) + (of_idx ? 1 - min(id_i) : 0);
    }
  }
  
  return p_max;
}

/**
 * Create position array from group sizes
 * 
 * @param n_x Array of group sizes
 * @return Position array where pos[i] is the starting index for group i
 */
array[] int create_pos(array[] int n_x) {
  return create_pos(n_x, 0);
}

/**
 * Create position array from group sizes with optional increment
 * 
 * @param n_x Array of group sizes
 * @param inc Additional increment to add to each group size
 * @return Position array with adjusted group sizes
 */
array[] int create_pos(array[] int n_x, int inc) {
  int n = size(n_x);
  array[n + 1] int pos;
  pos[1] = 1;  // Stan uses 1-based indexing
  
  // Build cumulative position array
  for (i in 1:n) {
    pos[i + 1] = pos[i] + n_x[i] + inc;  
  } 
  
  // Verify total size matches expectation
  // Cast to int to avoid type mismatch with optimized fma() function
  int expected_size = sum(n_x) + n * inc;
  assert_equal(pos[n + 1] - 1, expected_size);
  
  return pos;
}

/**
 * Create position array from subset of groups
 * 
 * @param n_x Array of all group sizes
 * @param sub_pos Position array defining which groups to include
 * @return Position array for the subset of groups
 */
array[] int create_pos(array[] int n_x, array[] int sub_pos) {
  int n = size(sub_pos) - 1;
  array[n + 1] int pos;
  pos[1] = 1;
  
  // Build position array by summing sizes of selected groups
  for (i in 1:n) {
    pos[i + 1] = pos[i] + sum(get_int_sub_array(n_x, sub_pos, i));  
  } 
  
  assert_equal(pos[n + 1] - 1, sum(n_x));
  
  return pos;
}

/**
 * Create position array for a range of groups
 * 
 * @param pos Original position array
 * @param from Starting group index
 * @param to Ending group index (inclusive)
 * @return Position array for groups from 'from' to 'to'
 */
array[] int create_pos(array[] int pos, int from, int to) {
  return create_pos(get_pos_size(pos)[from:to]);
}

tuple(int, array[] int) create_double_pos(int n_outer, array[] int inner) {
  return (n_outer, create_pos(inner));
}

/**
 * Get start and end indices for a range of groups
 * 
 * @param pos Position array
 * @param from Starting group index
 * @param to Ending group index (inclusive)
 * @return Tuple of (start_index, end_index) for the range
 */
tuple(int, int) get_pos(array[] int pos, int from, int to) {
  return (pos[from], pos[to + 1] - 1);
}

/**
 * Get start and end indices for a single group
 * 
 * @param pos Position array
 * @param n Group index
 * @return Tuple of (start_index, end_index) for group n
 */
tuple(int, int) get_pos(array[] int pos, int n) {
  return get_pos(pos, n, n);
}

tuple(int, int) get_offset_pos(array[] int pos, int n, int inner_offset) {
  tuple(int, int) start_end_pos = get_pos(pos, n);
  start_end_pos.1 += inner_offset;
  start_end_pos.2 += inner_offset;
  
  return start_end_pos;
}

tuple(int, int) get_pos(tuple(int, array[] int) double_pos, int m, int n) {
  int output_start_m1 = double_pos.1 * (m - 1); 
  tuple(int, int) pos = get_pos(double_pos.2, n);
  
  pos.1 += output_start_m1;
  pos.2 += output_start_m1;
  
  return pos;
} 

tuple(int, int, int, int) get_visit_pos(array[] int pos, int i, int n_screening) {
  int visit_start, screening_visit_end, treat_visit_start, visit_end;
  (visit_start, visit_end) = get_pos(pos, i);
  treat_visit_start = visit_start + n_screening; 
  screening_visit_end = treat_visit_start - 1;

  if (visit_end >= treat_visit_start && treat_visit_start >= screening_visit_end && screening_visit_end >= visit_start) {
    return(visit_start, screening_visit_end, treat_visit_start, visit_end);
  } else {
    fatal_error("Unexpected order of positions");
  }
} 
 
/**
 * Get size of a specific group
 * 
 * @param pos Position array
 * @param i Group index
 * @return Number of elements in group i
 */
int get_pos_size(array[] int pos, int i) {
  return pos[i + 1] - pos[i];
}

/**
 * Get sizes of all groups
 * 
 * @param pos Position array
 * @return Array of group sizes
 */
array[] int get_pos_size(array[] int pos) {
  int n = size(pos) - 1;
  array[n] int sizes;
  
  for (i in 1:n) {
    sizes[i] = get_pos_size(pos, i);
  }
  
  return sizes;
}

/**
 * Get total number of elements across all groups
 * 
 * @param pos Position array
 * @return Total number of elements
 */

int get_pos_total_size(array[] int pos) {
  return pos[size(pos)] - 1;
}

/**
 * Resize a ragged array by adding/removing elements from each group
 * 
 * @param full Original array
 * @param pos Position array for original array
 * @param inc Number of elements to add (positive) or remove (negative) from each group
 * @return Resized array with adjusted group sizes
 */
array[] int resize_int_array(array[] int full, array[] int pos, int inc) {
  int n = size(pos) - 1;
  array[n + 1] int new_pos = create_pos(pos, inc);
  int new_size = new_pos[n + 1] - 1;
  array[new_size] int new_array = zeros_int_array(new_size);
  
  // Copy data from old array to new array, respecting group boundaries
  for (i in 1:n) {
    int old_start, old_end;
    (old_start, old_end) = get_pos(pos, i);
    int new_start, new_end;
    (new_start, new_end) = get_pos(new_pos, i);
    
    // Copy data, handling cases where new size is smaller than old
    new_array[new_start:new_end] = full[old_start:max(min(old_end - old_start + 1, old_end + inc), old_end)];
  } 
  
  return new_array;
}

/**
 * Create a new position array for a subset of elements within a range
 * 
 * @param pos Original position array
 * @param from Starting index in the flat array
 * @param to Ending index in the flat array
 * @return Position array for groups that contain elements in the range [from, to]
 */
array[] int resize_pos(array[] int pos, int from, int to) {
  int n = size(pos) - 1;
  assert_greater_than_or_equal(to, from);  // Ensure valid range
  int first_group = 0, last_group = 0;
  int first_group_size, last_group_size;
  array[n] int new_pos_size = zeros_int_array(n);

  // Find which groups contain the range boundaries
  for (p in 1:n) {
    int start, end;
    (start, end) = get_pos(pos, p);
    
    // Find first group containing 'from'
    if (first_group == 0 && start <= from) {
      first_group = p;
      new_pos_size[p] = end - from + 1;  // Partial group from 'from' to end
    } 
    
    // Find last group containing 'to'
    if (last_group == 0 && end > to) {
      assert_greater_than_or_equal(first_group, 0);
      last_group = p - 1;
      new_pos_size[p] = to - start + 1;  // Partial group from start to 'to'
      break;
    } else if (first_group != 0) {
      new_pos_size[p] = end - start + 1;  // Full group size
    } 
  }
  
  return create_pos(new_pos_size);
}

/**
 * Extract integer sub-array for a specific group
 * 
 * @param full Full array containing all groups
 * @param pos Position array
 * @param n Group index
 * @return Sub-array containing only elements from group n
 */
array[] int get_int_sub_array(array[] int full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);
  
  return full[start:end];
}

/**
 * Extract integer sub-array for a range of groups
 * 
 * @param full Full array containing all groups
 * @param pos Position array
 * @param from Starting group index
 * @param to Ending group index (inclusive)
 * @return Sub-array containing elements from groups 'from' to 'to'
 */
array[] int get_int_sub_array(array[] int full, array[] int pos, int from, int to) {
  int from_start, from_end, to_start, to_end;
  (from_start, from_end) = get_pos(pos, from);
  (to_start, to_end) = get_pos(pos, to);
  
  return full[from_start:to_end];
}

/**
 * Extract real sub-array for a specific group
 * 
 * @param full Full array containing all groups
 * @param pos Position array
 * @param n Group index
 * @return Sub-array containing only elements from group n
 */
array[] real get_real_sub_array(array[] real full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);
  
  return full[start:end];
}

/**
 * Extract real sub-array for a range of groups
 * 
 * @param full Full array containing all groups
 * @param pos Position array
 * @param from Starting group index
 * @param to Ending group index (inclusive)
 * @return Sub-array containing elements from groups 'from' to 'to'
 */
array[] real get_real_sub_array(array[] real full, array[] int pos, int from, int to) {
  int from_start, from_end, to_start, to_end;
  (from_start, from_end) = get_pos(pos, from);
  (to_start, to_end) = get_pos(pos, to);
  
  return full[from_start:to_end];
}

/**
 * Extract matrix rows for a specific group
 * 
 * @param full Full matrix containing all groups
 * @param pos Position array (indexes rows)
 * @param n Group index
 * @return Sub-matrix containing only rows from group n
 */
matrix get_sub_vert_matrix(matrix full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);
  
  return full[start:end];
}

/**
 * Get minimum values for each group
 * 
 * @param x Array of values
 * @param pos Position array
 * @return Array of minimum values for each group
 */
array[] int get_min_pos(array[] int x, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int min_pos;
  
  for (i in 1:n) {
    min_pos[i] = get_min_pos(x, pos, i);
  }
  
  return min_pos;
}

/**
 * Get maximum values for each group
 * 
 * @param x Array of values
 * @param pos Position array
 * @return Array of maximum values for each group
 */
array[] int get_max_pos(array[] int x, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int max_pos;
  
  for (i in 1:n) {
    max_pos[i] = get_max_pos(x, pos, i);
  }
  
  return max_pos;
}

/**
 * Get minimum value for a specific group
 * 
 * @param x Array of values
 * @param pos Position array
 * @param n Group index
 * @return Minimum value in group n
 */
int get_min_pos(array[] int x, array[] int pos, int n) {
  return min(get_int_sub_array(x, pos, n));
}

/**
 * Get maximum value for a specific group
 * 
 * @param x Array of values
 * @param pos Position array
 * @param n Group index
 * @return Maximum value in group n
 */
int get_max_pos(array[] int x, array[] int pos, int n) {
  return max(get_int_sub_array(x, pos, n));
}

/**
 * Get specific element from a group by index
 * 
 * @param x Array of values
 * @param pos Position array
 * @param p Group index
 * @param n Element index within group (1-based)
 * @return Element at position n within group p
 */
int get_int(array[] int x, array[] int pos, int p, int n) { 
  int idx = pos[p] + n - 1;
  
  // Bounds checking
  if (idx >= pos[p + 1] || n < 1) {
    fatal_error("Unexpected index: ", n, " in group ", p, " of size ", get_pos_size(pos, p), " with pos: ", pos);
  }
  
  return x[pos[p] + n - 1]; 
}

array[] int validate_pos(array[] int pos) {
  int n = size(pos) - 1;
  array[n + 1] int sort_idx = sort_indices_asc(pos);

  for (i in 1:(n + 1)) {
    if (sort_idx[i] != i) {
      fatal_error("Invalid pos: ", pos);
    }
  }

  return pos;
}

/**
 * Get the global group index from level and local group ID
 *
 * For a flattened array containing all groups across all levels,
 * converts (level, local_group_id) to the global index.
 *
 * Example: With level_pos = [1, 6, 16, 166] (5 trials, 10 regions, 150 patients):
 *   get_global_group_idx(level_pos, 2, 3) = 6 + 3 - 1 = 8 (region 3)
 *
 * @param level_pos Position array for levels (from create_pos(n_groups_per_level))
 * @param level Level index (1-based)
 * @param group_id Local group ID within level (1-based)
 * @return Global group index in flattened group array
 */
int get_global_group_idx(array[] int level_pos, int level, int group_id) {
  return level_pos[level] + group_id - 1;
}

/**
 * Create position array for enabled levels only
 *
 * Creates a cumulative position array where disabled levels contribute 0 to the
 * positions. This allows parameters to be sized exactly for enabled levels only,
 * avoiding wasted memory and sampling for disabled levels.
 *
 * Example: With n_groups = [2, 732] and enabled = [1, 0]:
 *   Returns [1, 3, 3] - only level 1 contributes, level 2 contributes 0
 *
 * @param n_groups Array of group counts per level (e.g., [n_trials, n_patients])
 * @param enabled Array of 0/1 flags indicating which levels are enabled
 * @return Position array where pos[i] is the starting index for level i in the
 *         compacted (enabled-only) parameter array
 */
array[] int create_enabled_pos(array[] int n_groups, array[] int enabled) {
  int n = size(n_groups);
  array[n + 1] int pos;
  pos[1] = 1;
  for (i in 1:n) {
    pos[i + 1] = pos[i] + (enabled[i] ? n_groups[i] : 0);
  }
  return pos;
}

/**
 * Compute total number of enabled groups across all levels
 *
 * @param n_groups Array of group counts per level
 * @param enabled Array of 0/1 flags indicating which levels are enabled
 * @return Sum of group counts for enabled levels only
 */
int compute_n_enabled_groups(array[] int n_groups, array[] int enabled) {
  int n = size(n_groups);
  int total = 0;
  for (i in 1:n) {
    if (enabled[i]) {
      total += n_groups[i];
    }
  }
  return total;
}