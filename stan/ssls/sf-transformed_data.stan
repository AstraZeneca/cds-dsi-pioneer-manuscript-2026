// Relocated from tumor/sf-transformed_data.stan
array[sf_rep_T] real rep_time_points;
if (sf_rep_T > 0) rep_time_points = linspaced_array(sf_rep_T, 1, sf_rep_T);
real<lower = 0> lod = 0.1; // cm
vector<lower = 0>[sum(n_patient_visits)] normalized_sld;
for (i in 1:n_patients) {
  int visit_pos, visit_end; (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
  normalized_sld[visit_pos:visit_end] = sum_tumor_size[visit_pos:visit_end] / sum_tumor_size[visit_pos]; 
}
int<lower = 1> n_total_visits_m1 = sum(n_patient_visits) - n_patients;
array[n_patients + 1] int<lower = 1> patient_visit_m1_pos = create_pos(n_patient_visits, -1);
array[n_patients] int<lower = 1> ub_pfs_p1; for (i in 1:n_patients) ub_pfs_p1[i] = pfs[i] + interval_censored[i] + 1;
real log_lod = log(0.1);
int CR = 1; int PR = 2; int SD = 3; int PD = 4;
int NT_CR = 1; int NT_STABLE = 2; int NT_PD = 3;
// Handle n_covar = 0 case (no covariates)
matrix[n_patients, n_covar] Q_covar_design_matrix;
matrix[n_covar, n_covar] R_covar_design_matrix;

if (n_covar > 0) {
  Q_covar_design_matrix = qr_thin_Q(covar_design_matrix) * sqrt(n_patients - 1);
  R_covar_design_matrix = qr_thin_R(covar_design_matrix) / sqrt(n_patients - 1);
}
