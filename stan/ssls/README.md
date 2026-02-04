# SSLS Model Refactor Skeleton

This directory houses refactored scaffolding for the SSLS tumor model components.

## Layout
```
stan/ssls/
  modules/
    tr/    (total rate)
      flags.stan
      parameters.stan
      transformed_parameters.stan
      priors.stan
      generated_quantities.stan
    frac/  (fraction mix)
      flags.stan
      parameters.stan
      transformed_parameters.stan
      priors.stan
      generated_quantities.stan
    init/  (initial proportions)
      flags.stan
      parameters.stan
      transformed_parameters.stan
      priors.stan
      generated_quantities.stan
  modules_aggregator.stan (usage examples & include ordering)
```

## Documentation

- **[MODULE_DESIGN.md](MODULE_DESIGN.md)** - Comprehensive inventory and design specification for tr (total rate), frac (fraction mix), and init (initial proportions) modules. See Sections 2.1–2.4 for canonical naming patterns and flag semantics.
- **[NAMING_CONVENTION.md](NAMING_CONVENTION.md)** - Naming convention compliance tracking and migration status.

## Migration Strategy
1. Populate data flags in existing main model (add includes).
2. Move current parameter declarations into the appropriate module parameter files (rename simultaneously).
3. Shift transformed parameter logic piecemeal (verify identical log density after each move).
4. Relocate priors to module priors files (ensuring hyperparam usage audit).
5. Introduce slope hierarchies only when enabled by flags; keep disabled code absent to avoid orphan parameters.

## Conventions
- `flags.stan` (singular) per module keeps naming concise (`frac/flags.stan` acceptable; kept plural for clarity—can rename to `flag.stan` if preferred globally).
- Stan built-in functions (`log_inv_logit`, `log1m_inv_logit`) are used for numerically stable logit transformations.

## Next Actions
Proceed with Task 6: hyperparameter alignment, then incrementally migrate code into these stubs.
