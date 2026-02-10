// modules_aggregator.stan
// Central include orchestrating module components.
// Intended usage pattern in main model:
// data    block:   include flags & shared covariate structures
// parameters block: include per-module parameter snippets
// transformed parameters: include per-module derivations
// model block: include per-module priors
// generated quantities: optional per-module diagnostics

// NOTE: Each included file must be syntactically valid inside its target block. Current skeleton
// uses commented placeholders until migration.

// Example (in main model file):
// data {
//   #include "ssls/modules/tr/flags.stan"
//   #include "ssls/modules/frac/flags.stan"
//   #include "ssls/modules/init/flags.stan"
// }
// parameters {
//   #include "ssls/modules/tr/parameters.stan"
//   #include "ssls/modules/frac/parameters.stan"
//   #include "ssls/modules/init/parameters.stan"
// }
// transformed parameters {
//   #include "ssls/modules/tr/transformed_parameters.stan"
//   #include "ssls/modules/frac/transformed_parameters.stan"
//   #include "ssls/modules/init/transformed_parameters.stan"
// }
// model {
//   #include "ssls/modules/tr/priors.stan"
//   #include "ssls/modules/frac/priors.stan"
//   #include "ssls/modules/init/priors.stan"
// }
// generated quantities {
//   #include "ssls/modules/tr/generated_quantities.stan"
//   #include "ssls/modules/frac/generated_quantities.stan"
//   #include "ssls/modules/init/generated_quantities.stan"
// }
