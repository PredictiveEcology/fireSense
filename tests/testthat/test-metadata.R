## The parent's metadata is its contract: the children it fetches, and their versions.
## Needs no network and none of the children.

md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)

childSpecs <- c(
  "PredictiveEcology/fireSense_ELFs@development",
  "PredictiveEcology/fireSense_dataPrepFit@development",
  "PredictiveEcology/fireSense_ignitionFit@development",
  "PredictiveEcology/fireSense_spreadFit@development",
  "PredictiveEcology/fireSense_dataPrepPredict@development",
  "PredictiveEcology/fireSense_ignitionPredict@development",
  "PredictiveEcology/fireSense_spreadPredict@development",
  "PredictiveEcology/fireSense_burn@development",
  "PredictiveEcology/fireSense_summary@modsForFireSense"
)
childNames <- c("fireSense_ELFs", "fireSense_dataPrepFit", "fireSense_ignitionFit",
                "fireSense_spreadFit", "fireSense_dataPrepPredict", "fireSense_ignitionPredict",
                "fireSense_spreadPredict", "fireSense_burn", "fireSense_summary")

test_that("metadata parses and names the parent", {
  expect_identical(md$name, "fireSense")
  expect_identical(md$timeunit, "year")
  expect_identical(md$version$fireSense, "1.0.0")
})

test_that("a parent has no parameters, inputs or outputs", {
  expect_equal(NROW(md$parameters), 0L)
  expect_equal(NROW(md$inputObjects), 0L)
  expect_equal(NROW(md$outputObjects), 0L)
})

test_that("childModules are the nine child specs", {
  expect_identical(md$childModules, childSpecs)
})

test_that("the names derived from the specs are the nine module names", {
  skip_if_not_installed("Require")
  expect_identical(unname(Require::extractPkgName(md$childModules)), childNames)
})

test_that("version has an entry for the parent and each child", {
  expect_setequal(names(md$version), c("fireSense", childNames))
})
