// ============================================================================
// Endpoint Transformed Data
// Shared between full model and standalone multistate model.
// ============================================================================

// Position array for conditional group patient indices
array[n_cond_group + 1] int cond_group_pos = create_pos(cond_group_size);
