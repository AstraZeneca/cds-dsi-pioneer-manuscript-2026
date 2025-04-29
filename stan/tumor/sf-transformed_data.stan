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
