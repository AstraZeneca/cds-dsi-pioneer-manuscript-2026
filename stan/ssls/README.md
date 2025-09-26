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
    shared/
      helpers.stan
      covar_qr.stan
  modules_aggregator.stan (usage examples & include ordering)
```

## Naming & Flags
See `docs/tumor_total_fraction_init_inventory.md` Sections 2.1–2.4 for the canonical source of naming patterns and flag semantics.

## Migration Strategy
1. Populate data flags in existing main model (add includes).
2. Move current parameter declarations into the appropriate module parameter files (rename simultaneously).
3. Shift transformed parameter logic piecemeal (verify identical log density after each move).
4. Relocate priors to module priors files (ensuring hyperparam usage audit).
5. Introduce slope hierarchies only when enabled by flags; keep disabled code absent to avoid orphan parameters.

## Conventions
- `flags.stan` (singular) per module keeps naming concise (`frac/flags.stan` acceptable; kept plural for clarity—can rename to `flag.stan` if preferred globally).
- All function-style helpers go in `shared/helpers.stan` (will wrap in a `functions {}` block once populated in use sites).
- QR decomposition objects reside in `shared/covar_qr.stan` once extracted from legacy `base_data.stan`.

## Next Actions
Proceed with Task 6: hyperparameter alignment, then incrementally migrate code into these stubs.
