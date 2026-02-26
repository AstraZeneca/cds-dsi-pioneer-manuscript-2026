// ============================================================================
// Base Hierarchy Transformed Data
// Shared across all models (full joint model, standalone multistate, etc.)
// Requires: _base_hierarchy_data.stan declarations in scope.
// ============================================================================

// Multi-level hierarchy: computed totals and position arrays
int n_total_groups = sum(n_groups_per_level);
array[n_levels + 1] int level_pos = create_pos(n_groups_per_level);

// Validate multi-level hierarchy structure
for (p in 1:n_patients) {
  for (lv in 1:n_levels) {
    if (patient_level_groups[p, lv] < 1 ||
        patient_level_groups[p, lv] > n_groups_per_level[lv]) {
      fatal_error("Patient ", p, " has invalid group assignment at level ", lv,
                  ": got ", patient_level_groups[p, lv],
                  ", expected 1-", n_groups_per_level[lv]);
    }
  }
  // Verify patient level is identity mapping
  if (patient_level_groups[p, n_levels] != p) {
    fatal_error("Patient ", p, " must have identity mapping at patient level (level ",
                n_levels, "), got ", patient_level_groups[p, n_levels]);
  }
}

// Validate patient level has correct group count
if (n_groups_per_level[n_levels] != n_patients) {
  fatal_error("Patient level must have n_groups = n_patients, got ",
              n_groups_per_level[n_levels], " vs ", n_patients);
}

// Validate backward compatibility: patient_level_groups[,1] must match patient_trial
for (p in 1:n_patients) {
  if (patient_level_groups[p, 1] != patient_trial[p]) {
    fatal_error("patient_level_groups[", p, ", 1] = ", patient_level_groups[p, 1],
                " must match patient_trial[", p, "] = ", patient_trial[p]);
  }
}

print("Multi-level hierarchy validated:");
print("  n_levels = ", n_levels);
print("  n_groups_per_level = ", n_groups_per_level);

// Starting position of patients for each trial in a flattened patient array
// Diagram for trial_patient_pos:
// [1, 40, 70, 120, ...]
//  ^  ^   ^   ^
//  |  |   |   |
//  |  |   |   Start of patients in trial 4
//  |  |   Start of patients in trial 3
//  |  Start of patients in trial 2
//  Start of patients in trial 1
array[n_trials + 1] int<lower=1, upper=n_patients + 1> trial_patient_pos;
{
  array[n_trials] int n_trial_patients = rep_array(0, n_trials);
  for (i in 1:n_patients) {
    n_trial_patients[patient_trial[i]] += 1;
  }
  print("n_trial_patients = ", n_trial_patients);
  trial_patient_pos = create_pos(n_trial_patients);
}

real delta = 1e-5; // Small value used for GP modeling to avoid numerical issues
