---
title: "fireSense Manual"
subtitle: "v.1.1.0"
date: "Last updated: 2026-10-08"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
always_allow_html: true
---

# fireSense Module



#### Authors:

Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre]; Jean Marchal <jean.d.marchal@gmail.com> [aut]; Alex M. Chubaty <achubaty@for-cast.ca> [ctb]

## Module Overview

fireSense is a parent module. It has no events, parameters, inputs or outputs of its own; it lists the nine modules below as its children, and listing it in `setupProject(modules = )` fetches and runs all of them. Set parameters on the children, by their own names.

The family first fits, then predicts, then burns and summarises:

- `fireSense_ELFs`: defines the ELF study areas and keeps a ledger of which ELFs have been fitted, so an ELF is fitted once and reused.
- `fireSense_dataPrepFit`: prepares fire, vegetation and climate data for fitting.
- `fireSense_ignitionFit`: fits the ignition and escape models.
- `fireSense_spreadFit`: fits the spread model.
- `fireSense_dataPrepPredict`: prepares the covariates for prediction each year.
- `fireSense_ignitionPredict`: predicts ignitions and escapes from the fitted model.
- `fireSense_spreadPredict`: predicts spread probability from the fitted model.
- `fireSense_burn`: burns fires across the landscape using those predictions.
- `fireSense_summary`: summarises the fitted models and the burns.

Each module's own manual documents its inputs, outputs and parameters.

### Usage


``` r
SpaDES.project::setupProject(
  modules = "PredictiveEcology/fireSense@development",
  paths = list(projectPath = "myProject")
)
```

### Child modules


|child                              |
|:----------------------------------|
|fireSense_ELFs                     |
|fireSense_dataPrepFit              |
|fireSense_ignitionFit              |
|fireSense_spreadFit                |
|fireSense_dataPrepPredict          |
|fireSense_ignitionPredict          |
|fireSense_spreadPredict            |
|fireSense_burn                     |
|fireSense_summary@modsForFireSense |

### Migrating from the individual modules

Children are listed by name, so each comes from the parent's GitHub account and branch: `PredictiveEcology/fireSense@development` fetches every child's `development` branch, and `PredictiveEcology/fireSense@main` every child's `main`. A child written with a branch (`fireSense_summary@modsForFireSense`) or an account keeps it. A module the user lists in `setupProject(modules = )` overrides the parent's entry for it, e.g. `c("PredictiveEcology/fireSense@development", "PredictiveEcology/fireSense_burn@testing")`.

To migrate: list `PredictiveEcology/fireSense@development` in `setupProject(modules = ...)` instead of the individual fireSense modules, and rename these `params` keys (and `whichModulesToPrepare` values). `fireSense` now means the whole family; the burn module is `fireSense_burn`.

| old module name | new module name |
|---|---|
| fireSense (the burn module) | fireSense_burn |
| fireSense_IgnitionFit | fireSense_ignitionFit |
| fireSense_SpreadFit | fireSense_spreadFit |
| fireSense_IgnitionPredict | fireSense_ignitionPredict |
| fireSense_SpreadPredict | fireSense_spreadPredict |
| fireSense_ELFs, fireSense_dataPrepFit, fireSense_dataPrepPredict, fireSense_summary | unchanged |

### Getting help

Open an issue at <https://github.com/PredictiveEcology/fireSense/issues>.
