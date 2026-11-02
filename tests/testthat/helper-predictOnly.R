## End-to-end predict-only runs (test-predictOnly.R): a user adds fireSense to a project only to
## predict, from the fits stored in the production ledgers. They download real data and read the
## ledgers on Google Drive, so they run only when FIRESENSE_PREDICT_E2E is "true" (see
## predictOnlySkip()).

## The two ELFs: 6.6.1 and 6.6.2 are neighbours (a 535 km border) and both have stored spread and ignition
## fits; the default ledger, "latest" with no fitFileSuffix, takes each from the newest ledger file that
## has it. The study area is a 60 x 60 km box across their border, so one run
## predicts with both ELFs' fits (and blends them), as a user's own study area would.
predictOnlyELFs <- c("6.6.1", "6.6.2")
predictOnlyCRS <- "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
predictOnlyBox <- function() {
  terra::vect(terra::ext(1331940, 1391940, 6414900, 6474900), crs = predictOnlyCRS)
}

## The gate. The runs need the Google Drive login of `~/googledriveAuthentication.R` (or one in
## FIRESENSE_GOOGLE_AUTH) and the canClimateData and climateYear modules beside this one.
predictOnlySkip <- function() {
  testthat::skip_if_not(identical(Sys.getenv("FIRESENSE_PREDICT_E2E"), "true"),
                        "FIRESENSE_PREDICT_E2E is not \"true\" (long test: real data, Google Drive)")
  testthat::skip_if_not_installed("SpaDES.project")
}

## The module list of each variant, as a user writes it. "parent": the fireSense parent, so the fit modules
## come too (their default doFit = FALSE predicts from the stored fits). "noFit": only the modules
## prediction needs, with Biomass_borealDataPrep for the vegetation that fireSense_dataPrepFit would
## otherwise build. A module already in `modulePath` is used as it is (no download).
predictOnlyModules <- function(variant = c("parent", "noFit")) {
  pe <- function(m, branch = "development") paste0("PredictiveEcology/", m, "@", branch)
  climate <- c(pe("canClimateData", "modsForFireSense"), pe("climateYear", "modsForFireSense"))
  switch(match.arg(variant),
         parent = c(climate, pe("fireSense")),
         noFit = c(climate, pe(c("Biomass_borealDataPrep", "fireSense_ELFs", "fireSense_dataPrepPredict",
                                 "fireSense_ignitionPredict", "fireSense_spreadPredict", "fireSense_burn"))))
}

## setupProject() for one variant, in `root` (inputs, outputs, cache and scratch below it). Only what a
## predict-only user must set: the study area and the years.
predictOnlySetup <- function(variant, root, modulePath, times = list(start = 2021, end = 2023),
                             authFile = Sys.getenv("FIRESENSE_GOOGLE_AUTH", "~/googledriveAuthentication.R")) {
  sal <- predictOnlyBox()
  SpaDES.project::setupProject(
    paths = list(projectPath = root, modulePath = modulePath,
                 packagePath = .libPaths()[1],
                 inputPath = file.path(root, "inputs"),
                 outputPath = file.path(root, "outputs", variant),
                 cachePath = file.path(root, "cache"),
                 scratchPath = file.path(root, "scratch")),
    modules = predictOnlyModules(variant),
    times = times,
    studyAreaLarge = sal,
    options = list(authFile,
                   reproducible.cacheSaveFormat = "qs2",
                   reproducible.shapefileRead = "terra::vect",
                   reproducible.useMemoise = FALSE,
                   spades.useRequire = FALSE,
                   spades.moduleCodeChecks = FALSE,
                   spades.allowInitDuringSimInit = TRUE),
    params = list(.globals = list(.plots = NA),
                  ## only the simulated years' projected climate (the default is 2017:2100)
                  canClimateData = list(projectedClimateYears = times$start:times$end)),
    standAlone = FALSE, updateRprofile = FALSE, useGit = FALSE, Restart = FALSE
  )
}
