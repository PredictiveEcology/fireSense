# fireSense

fireSense is a parent SpaDES module: it has no events, parameters, inputs or outputs of its own. It lists the nine modules of the FireSense family of fire models as `childModules`, so that listing it in `setupProject(modules = )` fetches and runs all of them. See `fireSense.Rmd` for an overview of what each child does.

Children are listed by name, so each comes from the parent's GitHub account and branch: `PredictiveEcology/fireSense@development` fetches every child's `development` branch, and `PredictiveEcology/fireSense@main` every child's `main`. At a release tag, `PredictiveEcology/fireSense@v1.1.0`, each child comes at the version this module's `version` list gives it (`fireSense_burn = "2.2.0"` gives `fireSense_burn@v2.2.0`); this needs SpaDES.project >= 1.2.0.9013. A child written with a branch (`fireSense_summary@modsForFireSense`) or an account keeps it. A module the user lists in `setupProject(modules = )` overrides the parent's entry for it, e.g. `c("PredictiveEcology/fireSense@development", "PredictiveEcology/fireSense_burn@testing")`.

## Use

```r
SpaDES.project::setupProject(
  modules = "PredictiveEcology/fireSense@development",
  ...
)
```

## Migrating from the individual modules

To migrate: list `PredictiveEcology/fireSense@development` in `setupProject(modules = ...)` instead of the individual fireSense modules, and rename these `params` keys (and `whichModulesToPrepare` values). `fireSense` now means the whole family; the burn module is `fireSense_burn`.

| old module name | new module name |
|---|---|
| fireSense (the burn module) | fireSense_burn |
| fireSense_IgnitionFit | fireSense_ignitionFit |
| fireSense_SpreadFit | fireSense_spreadFit |
| fireSense_IgnitionPredict | fireSense_ignitionPredict |
| fireSense_SpreadPredict | fireSense_spreadPredict |
| fireSense_ELFs, fireSense_dataPrepFit, fireSense_dataPrepPredict, fireSense_summary | unchanged |
