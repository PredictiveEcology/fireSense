## burn(): with spread probability 1 inside a block and 0 on the barrier, spread2() is
## deterministic -- a fire burns its whole block and nothing else -- so every expected
## value below is derived by hand from the toy landscape in helper-toy.R.

ig <- function(pixelID, escapes) data.table::data.table(pixelID = pixelID, escapes = escapes)

test_that("one escaped fire burns exactly its own block", {
  sim <- runFireSense(ig(1L, 1L)) # cell 1 is in the west block (cols 1-4)
  expect_identical(which(vals(sim$rstCurrentBurn) == 1), westCells())   # 40 cells
  expect_identical(sort(sim$burnDT$pixels), westCells())
  expect_identical(unique(sim$burnDT$initialPixels), 1L)
  ## nothing east of the barrier, and nothing on it
  expect_true(all(is.na(vals(sim$rstCurrentBurn)[c(barrierCells(), eastCells())])))
})

test_that("burnSummary, burnDT and the rasters agree with each other", {
  sim <- runFireSense(ig(c(1L, 10L), c(1L, 1L))) # cell 10 is in the east block (cols 6-10)
  bs <- sim$burnSummary
  expect_identical(names(bs), c("igLoc", "N", "year", "areaBurnedHa"))
  expect_identical(bs$igLoc, c(1L, 10L))
  expect_identical(bs$N, c(40L, 50L))                 # 4 cols x 10 rows; 5 cols x 10 rows
  expect_identical(as.numeric(bs$year), c(1, 1))
  expect_equal(bs$areaBurnedHa, c(250, 312.5))        # 40 * 6.25 ha; 50 * 6.25 ha

  ## the same numbers from the rasters and from burnDT
  expect_identical(sum(vals(sim$rstCurrentBurn) == 1, na.rm = TRUE), sum(bs$N))
  expect_identical(nrow(sim$burnDT), sum(bs$N))
  expect_equal(sum(bs$areaBurnedHa),
               sum(vals(sim$rstCurrentBurn) == 1, na.rm = TRUE) *
                 prod(terra::res(sim$rasterToMatch)) / 1e4)
  expect_identical(as.vector(table(sim$burnDT$initialPixels)), c(40L, 50L))
})

test_that("rstAnnualBurnID labels each fire with its fire_id", {
  sim <- runFireSense(ig(c(1L, 10L), c(1L, 1L)))
  id <- vals(sim$rstAnnualBurnID)
  expect_true(all(id[westCells()] == 1))  # fire_id is .GRP by initialPixels, in order of appearance
  expect_true(all(id[eastCells()] == 2))
  expect_true(all(is.na(id[barrierCells()])))
  expect_identical(sort(unique(sim$burnDT$fire_id)), 1:2)
  expect_identical(unique(sim$burnDT[initialPixels == 10L]$fire_id), 2L)
})

test_that("non-flammable pixels never burn and stay NA in burnMap", {
  sim <- runFireSense(ig(c(1L, 10L), c(1L, 1L)))
  expect_false(any(barrierCells() %in% sim$burnDT$pixels))
  expect_true(all(is.na(vals(sim$burnMap)[barrierCells()])))
  expect_identical(vals(sim$burnMap)[-barrierCells()], rep(1, 90))
})

test_that("a barrier of NA spread probability also stops fire", {
  sim <- runFireSense(ig(1L, 1L),
                      objects = list(fireSense_SpreadPredicted = toySpreadProb(barrier = NA)))
  expect_identical(sort(sim$burnDT$pixels), westCells())
})

test_that("ignitions that did not escape do not burn", {
  ## cell 10 ignited but did not escape: only the west block burns
  sim <- runFireSense(ig(c(1L, 10L), c(1L, 0L)))
  expect_identical(sim$burnSummary$igLoc, 1L)
  expect_identical(which(vals(sim$rstCurrentBurn) == 1), westCells())
  expect_identical(vals(sim$burnMap)[eastCells()], rep(0, 50))
})

test_that("no escapes at all leaves every output as init made it", {
  sim <- runFireSense(ig(c(1L, 10L), c(0L, 0L)))
  expect_null(sim$burnSummary)
  expect_null(sim$burnDT)
  expect_null(sim$rstAnnualBurnID)
  expect_false(terra::hasValues(sim$rstCurrentBurn)) # rast(rasterToMatch): geometry only
  expect_identical(sum(vals(sim$burnMap), na.rm = TRUE), 0)
})

test_that("an empty or all-NA ignitions table burns nothing", {
  empty <- runFireSense(ig(integer(), integer()))
  expect_null(empty$burnSummary)
  expect_identical(sum(vals(empty$burnMap), na.rm = TRUE), 0)

  allNA <- runFireSense(ig(1L, NA_integer_)) # sum(..., na.rm = TRUE) is 0
  expect_null(allNA$burnSummary)
})

test_that("nothing spreads unless whichModulesToPrepare includes fireSense_SpreadPredict", {
  sim <- runFireSense(ig(1L, 1L),
                      params = list(whichModulesToPrepare = "fireSense_IgnitionPredict"))
  expect_null(sim$burnSummary)
  expect_identical(sum(vals(sim$burnMap), na.rm = TRUE), 0)
})

test_that("zero spread probability burns only the ignition pixel", {
  sim <- runFireSense(ig(c(12L, 78L), c(1L, 1L)),
                      objects = list(fireSense_SpreadPredicted = toySpreadProb(p = 0)))
  expect_identical(sort(sim$burnDT$pixels), c(12L, 78L))
  expect_identical(sim$burnSummary$N, c(1L, 1L))
  expect_equal(sim$burnSummary$areaBurnedHa, c(6.25, 6.25))   # 1 pixel = 250 * 250 / 1e4 ha
})

test_that("fire spreads in 8 directions", {
  ## probability 1 only on the main diagonal: the 10 diagonal cells touch only at corners,
  ## so all 10 burn with 8 neighbours and only the ignition cell would with 4
  rtm <- toyRTM()
  sp <- rtm
  terra::values(sp) <- 0
  diagCells <- terra::cellFromRowCol(rtm, 1:10, 1:10)
  sp[diagCells] <- 1
  flam <- rtm # all flammable
  sim <- runFireSense(ig(diagCells[1], 1L),
                      objects = list(fireSense_SpreadPredicted = sp, flammableRTM = flam))
  expect_identical(sort(sim$burnDT$pixels), as.integer(diagCells))
  expect_equal(sim$burnSummary$areaBurnedHa, 62.5)            # 10 * 6.25 ha
})

test_that("two fires in one block share it without overlap", {
  sim <- runFireSense(ig(c(1L, 94L), c(1L, 1L))) # cells 1 and 94 are both west (cols 1 and 4)
  expect_identical(sort(sim$burnDT$pixels), westCells())      # each pixel burns once
  expect_identical(sum(sim$burnSummary$N), 40L)
  expect_true(all(sim$burnSummary$N >= 1L))
  expect_setequal(sim$burnSummary$igLoc, c(1L, 94L))
  ## a pixel belongs to the fire recorded in rstAnnualBurnID
  expect_identical(vals(sim$rstAnnualBurnID)[sim$burnDT$pixels], as.numeric(sim$burnDT$fire_id))
})

test_that("two escapes on one pixel burn the block once", {
  ## spread2() cannot start twice from one pixel, so the duplicate gets its own call and
  ## identical rows are then dropped. (With 3+ escapes on a pixel spread2() stops with
  ## "start has duplicates"; that is not pinned here.)
  sim <- runFireSense(ig(1L, 2L))
  expect_identical(sort(sim$burnDT$pixels), westCells())
  expect_identical(sim$burnSummary$N, 40L)
  expect_identical(vals(sim$burnMap)[westCells()], rep(1, 40)) # not 2
})

test_that("burnMap and burnSummary accumulate over years", {
  sim <- runFireSense(ig(1L, 1L), times = list(start = 1, end = 3))
  bm <- vals(sim$burnMap)
  expect_identical(bm[westCells()], rep(3, 40))                # burned in years 1, 2, 3
  expect_identical(bm[eastCells()], rep(0, 50))
  expect_identical(as.numeric(sim$burnSummary$year), c(1, 2, 3))
  expect_identical(sim$burnSummary$N, rep(40L, 3))
  expect_equal(sum(sim$burnSummary$areaBurnedHa), 750)         # 3 * 40 * 6.25 ha
  expect_identical(nrow(sim$burnDT), 40L)                      # burnDT holds the last year only
})

test_that("the pixel size sets the area burned", {
  rtm <- toyRTM(res = 100) # 100 * 100 / 1e4 = 1 ha per pixel
  sim <- runFireSense(ig(1L, 1L),
                      objects = list(rasterToMatch = rtm, flammableRTM = toyFlammable(rtm),
                                     fireSense_SpreadPredicted = toySpreadProb(rtm)))
  expect_equal(sim$burnSummary$areaBurnedHa, 40)
})

test_that("stochastic spread is reproducible and stays within the rules", {
  run <- function() {
    set.seed(4321)
    runFireSense(ig(c(1L, 10L), c(1L, 1L)),
                 objects = list(fireSense_SpreadPredicted = toySpreadProb(p = 0.3)))
  }
  a <- run()
  b <- run()
  expect_identical(a$burnDT$pixels, b$burnDT$pixels)
  expect_identical(a$burnSummary, b$burnSummary)
  ## a fire never leaves its block, never burns the barrier, always burns its ignition pixel
  expect_true(all(a$burnDT[initialPixels == 1L]$pixels %in% westCells()))
  expect_true(all(a$burnDT[initialPixels == 10L]$pixels %in% eastCells()))
  expect_true(all(c(1L, 10L) %in% a$burnDT$pixels))
  expect_identical(anyDuplicated(a$burnDT$pixels), 0L)
  expect_lt(nrow(a$burnDT), 90L) # p = 0.3 is far below the ~0.5 needed to fill a block
})
