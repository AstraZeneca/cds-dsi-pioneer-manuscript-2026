array[] int get_max(array[] int id, array[] int pos) {
  return get_max(id, pos, 0);
}

array[] int get_max_idx(array[] int id, array[] int pos) {
  return get_max(id, pos, 1);
}

array[] int get_max(array[] int id, array[] int pos, int of_idx) {
  int n = size(pos) - 1;
  array[n] int p_max = zeros_int_array(n);
  
  for (i in 1:n) {
    int n_i = get_pos_size(pos, i);
    
    if (n_i > 0) {
      array[n_i] int id_i = get_int_sub_array(id, pos, i);
      p_max[i] = max(id_i) + (of_idx ? 1 - min(id_i) : 0);
    }
  }
  
  return p_max;
}

array[] int create_pos(array[] int n_x) {
  int n = size(n_x);
  array[n + 1] int pos;
  pos[1] = 1;
  
  for (i in 1:n) {
    pos[i + 1] = pos[i] + n_x[i];  
  } 
  
  assert_equal(pos[n + 1] - 1, sum(n_x));
  
  return pos;
}

array[] int create_pos(array[] int n_x, array[] int sub_pos) {
  int n = size(sub_pos) - 1;
  array[n + 1] int pos;
  pos[1] = 1;
  
  for (i in 1:n) {
    pos[i + 1] = pos[i] + sum(get_int_sub_array(n_x, sub_pos, i));  
  } 
  
  assert_equal(pos[n + 1] - 1, sum(n_x));
  
  return pos;
}

array[] int create_pos(array[] int pos, int from, int to) {
  return create_pos(get_pos_size(pos)[from:to]);
}

tuple(int, int) get_pos(array[] int pos, int from, int to) {
  return (pos[from], pos[to + 1] - 1);
}

tuple(int, int) get_pos(array[] int pos, int n) {
  return get_pos(pos, n, n);
}

int get_pos_size(array[] int pos, int i) {
  return pos[i + 1] - pos[i];
}

array[] int get_pos_size(array[] int pos) {
  int n = size(pos) - 1;
  array[n] int sizes;
  
  for (i in 1:n) {
    sizes[i] = get_pos_size(pos, i);
  }
  
  return sizes;
} 

array[] int get_int_sub_array(array[] int full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);
  
  return full[start:end];
}

array[] int get_int_sub_array(array[] int full, array[] int pos, int from, int to) {
  int from_start, from_end, to_start, to_end;
  (from_start, from_end) = get_pos(pos, from);
  (to_start, to_end) = get_pos(pos, to);
  
  return full[from_start:to_end];
}

vector get_sub_vector(vector full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);
  
  return full[start:end];
}

row_vector get_sub_row_vector(vector full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);
  
  return full[start:end]';
}

array[] int get_min_pos(array[] int x, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int min_pos;
  
  for (i in 1:n) {
    min_pos[i] = get_min_pos(x, pos, i);
  }
  
  return min_pos;
}

array[] int get_max_pos(array[] int x, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int max_pos;
  
  for (i in 1:n) {
    max_pos[i] = get_max_pos(x, pos, i);
  }
  
  return max_pos;
}

int get_min_pos(array[] int x, array[] int pos, int n) {
  return min(get_int_sub_array(x, pos, n));
}

int get_max_pos(array[] int x, array[] int pos, int n) {
  return max(get_int_sub_array(x, pos, n));
}

int get_int(array[] int x, array[] int pos, int p, int n) { 
  int idx = pos[p] + n - 1;
  
  if (idx >= pos[p + 1] || n < 1) {
    fatal_error("Unexpected index: ", n);
  }
  
  return x[pos[p] + n - 1]; 
}

int get_last_int(array[] int x, array[] int pos, int p) {
  return get_int(x, pos, p, pos[p + 1] - 1);
}