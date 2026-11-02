## End-to-end: a user adds fireSense to a project only to predict, from the stored production fits of
## two neighbour ELFs (helper-predictOnly.R). Long (real data, Google Drive): runs only with
## FIRESENSE_PREDICT_E2E=true. FIRESENSE_PREDICT_E2E_DIR keeps inputs and cache between runs (default: a
## temporary folder, so every run downloads everything again). The modules must be in one folder, with
## canClimateData and climateYear beside this one (FIRESENSE_PREDICT_E2E_MODULEPATH, default `modulePath`).

predictOnlyRoot <- function() {
  root <- Sys.getenv("FIRESENSE_PREDICT_E2E_DIR")
  if (!nzchar(root)) root <- withr::local_tempdir(.local_envir = testthat::teardown_env())
  root
}

predictOnlyRun <- function(variant) {
  mp <- Sys.getenv("FIRESENSE_PREDICT_E2E_MODULEPATH", modulePath)
  out <- predictOnlySetup(variant, root = file.path(predictOnlyRoot(), variant), modulePath = mp)
  SpaDES.core::simInitAndSpades2(out)
}

## what every predict-only run must show: each simulated year burned from predictions made with both
## ELFs' stored fits
expectPredictedWithBothELFs <- function(sim, times = list(start = 2021, end = 2023)) {
  ev <- as.data.frame(SpaDES.core::completed(sim))
  burnYears <- ev$eventTime[ev$moduleName == "fireSense_burn" & ev$eventType != "init"]
  expect_setequal(unique(floor(burnYears)), times$start:times$end)
  expect_setequal(intersect(as.character(sim$studyAreaWithSpreadParams$polygonID), predictOnlyELFs),
                  predictOnlyELFs)
  expect_setequal(intersect(names(sim$fireSense_IgnitionFittedList), predictOnlyELFs), predictOnlyELFs)
  expect_s4_class(sim$fireSense_SpreadPredicted, "SpatRaster")
  expect_gt(sum(!is.na(terra::values(sim$fireSense_SpreadPredicted, mat = FALSE))), 0)
  expect_s4_class(sim$fireSense_IgAndEscapeProbRas, "SpatRaster")
  expect_s4_class(sim$burnMap, "SpatRaster")
  expect_true(data.table::is.data.table(sim$burnSummary))
}

test_that("the fireSense parent with the default doFit predicts from the stored fits, and fits nothing", {
  predictOnlySkip()
  sim <- predictOnlyRun("parent")
  expectPredictedWithBothELFs(sim)
  ## doFit = FALSE (the default) never fits: no threshold calibration and no DEoptim run, no xgboost fit
  ev <- as.data.frame(SpaDES.core::completed(sim))
  expect_false(any(ev$moduleName == "fireSense_spreadFit" & ev$eventType %in% c("estimateThreshold", "run")))
  expect_null(sim$DE)
})

test_that("without the fit modules (no fireSense_dataPrepFit), prediction reads the stored fits itself", {
  predictOnlySkip()
  sim <- predictOnlyRun("noFit")
  expectPredictedWithBothELFs(sim)
  ## the fits' climate variables, not fireSense_dataPrepPredict's former 'MDC' default
  expect_false("MDC" %in% unlist(sim$climateVariablesForFire))
})
