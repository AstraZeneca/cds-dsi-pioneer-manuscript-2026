// gr_decay/generated_quantities.stan — diagnostics for prior-vs-posterior contraction.
// pop_log_growth_rate is in scope (declared earlier in the model's GQ block).

real gr_decay_kappa_pop = enable_gr_decay ? exp(gr_decay_log_loc_pop[1]) : 0.0;
// Implied long-horizon plateau OFFSET on the log-growth arm: growth_rate / kappa.
real gr_decay_plateau_offset = enable_gr_decay
  ? exp(pop_log_growth_rate) / gr_decay_kappa_pop
  : positive_infinity();
