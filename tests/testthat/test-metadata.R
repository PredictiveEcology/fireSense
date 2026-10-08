## The parent's metadata is its contract: the children it fetches, and their versions.
## Needs no network and none of the children.

md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)

childSpecs <- c(
  "fireSense_ELFs",
  "fireSense_dataPrepFit",
  "fireSense_ignitionFit",
  "fireSense_spreadFit",
  "fireSense_dataPrepPredict",
  "fireSense_ignitionPredict",
  "fireSense_spreadPredict",
  "fireSense_burn",
  "fireSense_summary@modsForFireSense"
)
childNames <- c("fireSense_ELFs", "fireSense_dataPrepFit", "fireSense_ignitionFit",
                "fireSense_spreadFit", "fireSense_dataPrepPredict", "fireSense_ignitionPredict",
                "fireSense_spreadPredict", "fireSense_burn", "fireSense_summary")

test_that("metadata parses and names the parent", {
  expect_identical(md$name, "fireSense")
  expect_identical(md$timeunit, "year")
  expect_false(is.na(package_version(md$version$fireSense, strict = FALSE)))
})

test_that("the version list names a valid version for every child", {
  ## SpaDES.project fetches a child of fireSense@v<x> at v<its version here>
  expect_true(all(childNames %in% names(md$version)))
  expect_false(anyNA(package_version(unlist(md$version[childNames]), strict = FALSE)))
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
