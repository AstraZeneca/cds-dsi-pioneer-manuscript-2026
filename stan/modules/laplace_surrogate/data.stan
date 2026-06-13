// Fixed calendar anchors (weeks) for the surrogate quadratic. The first anchor
// MUST be 0 (baseline normalization pins g(0)=0; the surrogate intercept relies
// on it). Tune the other two to the backgrounded trials' observation window.
vector[3] surrogate_anchor_times;

// Inner Laplace-solver knobs, data-driven so they are tunable without a recompile
// (perf diagnostics + production). A non-positive value selects the built-in
// default (max_num_steps=100, tolerance=1e-8, solver=1). Assembled with defaults
// in prepare_tumor_stan_data(), so callers that don't set them are unaffected.
int surrogate_max_num_steps_in;
real surrogate_tolerance_in;
int surrogate_solver_in;
