## fireSense's `burn` must run after the same year's `run` of fireSense_IgnitionPredict
## (makes ignitionsAndEscapes) and fireSense_SpreadPredict (makes fireSense_SpreadPredicted).
## With equal time and priority, SpaDES breaks ties by module load order, and load order
## with no `loadOrder` metadata comes only from the dependency graph -- which has no edge
## from fireSense to these two here, since the stubs below declare no inputObjects or
## outputObjects. Only fireSense's own `loadOrder` metadata can fix that.

stubModulesPath <- testthat::test_path("stubModules")

loadOrderPaths <- function() {
  root <- withr::local_tempdir(.local_envir = testthat::teardown_env())
  paths <- list(cachePath = file.path(root, "cache"), inputPath = file.path(root, "inputs"),
                modulePath = c(dirname(toyModuleRoot()), stubModulesPath),
                outputPath = file.path(root, "outputs"))
  for (p in paths[c("cachePath", "inputPath", "outputPath")])
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  paths
}

test_that("fireSense loads after fireSense_IgnitionPredict and fireSense_SpreadPredict", {
  rtm <- toyRTM()
  objs <- list(rasterToMatch = rtm, flammableRTM = toyFlammable(rtm),
               fireSense_SpreadPredicted = toySpreadProb(rtm),
               ignitionsAndEscapes = data.table::data.table(pixelID = 1L, escapes = 0L, escaped = FALSE))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  sim <- suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 1),
    modules = c("fireSense", "fireSense_IgnitionPredict", "fireSense_SpreadPredict"),
    objects = objs, paths = loadOrderPaths()
  ))

  modLoadOrder <- unlist(sim@modules)
  fsPos <- which(modLoadOrder == "fireSense")
  ignPos <- which(modLoadOrder == "fireSense_IgnitionPredict")
  sprPos <- which(modLoadOrder == "fireSense_SpreadPredict")

  expect_gt(fsPos, ignPos)
  expect_gt(fsPos, sprPos)
})
