## Parent module: it has no events, parameters, inputs or outputs of its own.
## Listing it in `setupProject(modules = )` fetches and runs its child modules.
defineModule(sim, list(
  name = "fireSense",
  description = paste("The FireSense family of fire modules: ELF study areas, fit data prep,",
                      "ignition/escape fit, spread fit, predict data prep, ignition/escape predict,",
                      "spread predict, burning and summary. Listing this parent in",
                      "`setupProject(modules = )` fetches and runs all nine child modules."),
  keywords = c("fire", "fireSense", "parent module"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut"),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = c(
    "PredictiveEcology/fireSense_ELFs@development",
    "PredictiveEcology/fireSense_dataPrepFit@development",
    "PredictiveEcology/fireSense_ignitionFit@development",
    "PredictiveEcology/fireSense_spreadFit@development",
    "PredictiveEcology/fireSense_dataPrepPredict@development",
    "PredictiveEcology/fireSense_ignitionPredict@development",
    "PredictiveEcology/fireSense_spreadPredict@development",
    "PredictiveEcology/fireSense_burn@development",
    "PredictiveEcology/fireSense_summary@modsForFireSense"
  ),
  version = list(
    fireSense = "1.0.0",
    fireSense_ELFs = "1.1.7",
    fireSense_dataPrepFit = "1.2.0.9020",
    fireSense_ignitionFit = "1.1.0",
    fireSense_spreadFit = "1.1.0",
    fireSense_dataPrepPredict = "1.0.4.9009",
    fireSense_ignitionPredict = "1.1.0",
    fireSense_spreadPredict = "1.1.1",
    fireSense_burn = "2.1.0",
    fireSense_summary = "1.0.1.9002"
  ),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "fireSense.Rmd")
  ## a parent module has no reqdPkgs, parameters, inputObjects, nor outputObjects
))

## a parent module has no events.
