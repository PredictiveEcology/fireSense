# fireSense

fireSense is a parent SpaDES module: it has no events, parameters, inputs or outputs of its own. It lists the nine modules of the FireSense family of fire models as `childModules`, so that listing it in `setupProject(modules = )` fetches and runs all of them. See `fireSense.Rmd` for an overview of what each child does.

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
