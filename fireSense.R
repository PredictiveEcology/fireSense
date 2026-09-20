# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense",
  description = "A landscape fire model, sensitive to environmental changes (e.g.
                 weather and land-cover).",
  keywords = c("fire", "percolation", "environmental control", "feedback",
               "weather", "vegetation", "land-cover"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut"),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = numeric_version("2.0.2.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "fireSense.Rmd"), ## same file
  reqdPkgs = list("data.table", "ggplot2", "ggspatial", "terra"),
  parameters = rbind(
    defineParameter(".plots", "character|logical", default = NULL, ## TODO: use .plotInitialTime etc.
                    desc = "Passed to `types` in `Plots()`, e.g. \"screen\", \"png\". `NULL` or `NA` for no plots."),
    defineParameter("plotIgnitions", "logical", FALSE, NA, NA,
                    "Currently unused."),
    defineParameter(".plotInterval", "numeric", 10, NA, NA,
                    "Years between plots of annual fire IDs, cumulative burns and spread probability."),
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    "Time of the first `burn` event."),
    defineParameter(".runInterval", "numeric", 1, NA, NA,
                    "Years between `burn` events. `NA` burns once only."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Currently unused."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "Currently unused."),
    defineParameter("whichModulesToPrepare", "character",
                    default = c("fireSense_SpreadPredict", "fireSense_IgnitionPredict", "fireSense_EscapePredict"),
                    NA, NA,
                    "Fires spread only if this includes `fireSense_SpreadPredict`. Other values are ignored.")
  ),
  inputObjects = rbind(
    expectsInput("fireSense_SpreadPredicted", "SpatRaster",
                 "Per-pixel spread probability for the current year."),
    expectsInput("flammableRTM", "list", 
                 "Binary SpatRaster (1 = flammable, 0 = not). Non-flammable pixels are `NA` in `burnMap`."),
    expectsInput("ignitionsAndEscapes", "data.table",
                 "One row per ignited pixel, with `pixelID` and `escapes`, the number of escaped fires there."),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "Template raster for the study area, ideally buffered to limit fire edge effects.")
  ),
  outputObjects = rbind(
    createsOutput("burnDT", "data.table",
                  "`spread2()` output for the most recent fire year: one row per burned pixel, plus `fire_id`."),
    createsOutput("burnMap", "SpatRaster",
                  "Number of times each pixel has burned. `NA` where not flammable."),
    createsOutput("burnSummary", "data.table",
                  "One row per fire: `igLoc` (ignition pixel), `N` (pixels burned), `year`, `areaBurnedHa`."),
    createsOutput("rstAnnualBurnID", "SpatRaster",
                  "Fire ID of each pixel burned this year; `NA` elsewhere."),
    createsOutput("rstCurrentBurn", "SpatRaster",
                  "1 where burned this year; `NA` elsewhere.")
  )
))

#' Event dispatcher
#'
#' `init` creates `burnMap` and schedules the first `burn`; `burn` spreads this
#' year's escaped fires and reschedules itself.
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType `"init"` or `"burn"`.
#' @param debug Unused.
#'
#' @return The `simList`, invisibly.
doEvent.fireSense = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      ## bail early if there's a problem with ignition, escape, or spread rasters

      if (!is.null(sim$fireSense_SpreadPredicted))
        stopifnot(length(na.omit(sim$fireSense_SpreadPredicted[])) > 0)

      ## trying to avoid the raster warning no non-missing arguments to max
      sim$burnMap <- rast(sim$rasterToMatch)
      sim$rstCurrentBurn <- rast(sim$rasterToMatch)
      
      burnVals <- as.vector(sim$flammableRTM)
      sim$burnMap[burnVals == 0] <- NA
      sim$burnMap[burnVals == 1] <- 0
      rm(burnVals)
      gc()

      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "burn", 
                           eventPriority = 5.13)
    },
    burn = {
      sim <- burn(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "burn", 
                             eventPriority = 5.13)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  invisible(sim)
}

#' Spread this year's escaped fires
#'
#' Spreads fires with `SpaDES.tools::spread2()` from every pixel in
#' `sim$ignitionsAndEscapes` with `escapes > 0`, then updates the burn outputs.
#' Does nothing if there are no escapes.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly, with `burnDT`, `burnMap`, `burnSummary`,
#'   `rstAnnualBurnID` and `rstCurrentBurn` updated.
burn <- function(sim) {

  moduleName <- current(sim)$moduleName

  escaped <- sum(sim$ignitionsAndEscapes$escapes, na.rm = TRUE)

  if (escaped > 0L) {
    if ("fireSense_SpreadPredict" %in% P(sim)$whichModulesToPrepare) {
      ## Spread
      # Note: if none of the cells are active SpaDES.tools::spread2() returns spreadState unchanged
      successfulEscapes <- sim$ignitionsAndEscapes[escapes > 0]
      igLocs <- rep(successfulEscapes$pixelID, times = successfulEscapes$escapes)
      igLocsList <- list(igLocs)
      ## spread2 fails with duplicated start pixels, so duplicates get their own spread2 call
      if (any(duplicated(tail(igLocsList, 1)[[1]]))) { 
        len <- length(igLocsList)
        igLocsList[[len + 1]] <- 
          igLocsList[[len]][duplicated(igLocsList[[len]])]
        igLocsList[[len]] <- unique(igLocsList[[len]])
      }
      
      spreadStates <- Map(igLocs = igLocsList, function(igLocs) {
        spreadState <- SpaDES.tools::spread2(
          landscape = sim$fireSense_SpreadPredicted,
          spreadProb = sim$fireSense_SpreadPredicted,
          directions = 8L,
          start = igLocs,
          asRaster = FALSE)  
      })
      
      spreadState <- rbindlist(spreadStates) |> unique()
      spreadState[ , fire_id := .GRP, by = "initialPixels"] # Add an fire_id column
      sim$rstAnnualBurnID <- rast(sim$fireSense_SpreadPredicted)
      sim$rstCurrentBurn <- rast(sim$fireSense_SpreadPredicted)
      sim$rstAnnualBurnID[spreadState$pixels] <- spreadState$fire_id
      sim$rstCurrentBurn[spreadState$pixels] <- 1
      sim$burnMap[spreadState$pixels] <- sim$burnMap[spreadState$pixels] + 1
      par("pin" = pmax(par()$pin, 0)) # not sure why par$pin is negative
      if ((time(sim) - start(sim) ) %% P(sim)$.plotInterval < 1) {
        Plots(c(sim$rstAnnualBurnID |> setNames(paste0("Annual Fire IDs ", time(sim))),
                sim$burnMap |> setNames(paste0("Cumulative Burn Map ", time(sim))),
                sim$fireSense_SpreadPredicted |> setNames(paste0("Spread Probability Map ", time(sim)))),
              types = Par$.plots, filename = paste0("Annual Fire Maps ", time(sim)),
              deviceArgs = list(width = 10, height = 8, units = "in", res = 144))
      }
      
      #get fire year, pixels burned, area burned, poly ID of all burned pixels
      # Make burnSummary --> similar to SCFM
      sim$burnDT <- spreadState

      tempDT <- sim$burnDT[, .(.N), by = "initialPixels"]
      tempDT$year <- time(sim)
      tempDT$areaBurnedHa <- tempDT$N * prod(res(sim$fireSense_SpreadPredicted)) * 1e-4
      setnames(tempDT, c("initialPixels"), c("igLoc"))
      sim$burnSummary <- rbind(sim$burnSummary, tempDT)
    }
  }

  invisible(sim)
}
