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
  version = numeric_version("2.0.2.9002"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "fireSense.Rmd"), ## same file
  reqdPkgs = list("data.table", "ggplot2", "ggspatial", "terra"),
  parameters = rbind(
    defineParameter(".plots", "character|logical", default = NULL, ## TODO: use .plotInitialTime etc.
                    desc = "Passed to `types` in `Plots()`, e.g. \"screen\", \"png\". `NULL` or `NA` for no plots."),
    defineParameter(".plotInterval", "numeric", 10, NA, NA,
                    "Years between plots of annual fire IDs, cumulative burns and spread probability."),
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    "Time of the first `burn` event."),
    defineParameter(".runInterval", "numeric", 1, NA, NA,
                    "Years between `burn` events. `NA` burns once only."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used."),
    defineParameter("whichModulesToPrepare", "character",
                    default = c("fireSense_SpreadPredict", "fireSense_IgnitionPredict", "fireSense_EscapePredict"),
                    NA, NA,
                    "Fires spread only if this includes `fireSense_SpreadPredict`. Other values are ignored.")
  ),
  inputObjects = rbind(
    expectsInput("fireSense_SpreadPredicted", "SpatRaster",
                 "Per-pixel spread probability for the current year."),
    expectsInput("fireSense_SpreadSD", "SpatRaster|numeric",
                 paste("Sd of the per-year random effect on logit spread probability, as fitted by",
                       "fireSense_SpreadFit (`yearSpreadSD`): a raster aligned with `fireSense_SpreadPredicted`",
                       "(from fireSense_SpreadPredict, per ELF) or one number. Each year draws one z ~ N(0, 1);",
                       "all of that year's fires spread with plogis(qlogis(p) + z * sd). NULL or 0: no effect.")),
    expectsInput("flammableRTM", "SpatRaster", 
                 "Binary SpatRaster (1 = flammable, 0 = not). Non-flammable pixels are `NA` in `burnMap`."),
    expectsInput("ignitionsAndEscapes", "data.table",
                 "One row per ignited pixel, with `pixelID` and `escapes`, the number of escaped fires there."),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "Template raster for the study area, ideally buffered to limit fire edge effects.")
  ),
  outputObjects = rbind(
    createsOutput("burnDT", "data.table",
                  paste("The most recent fire year's burned pixels: `initialPixels` (the fire's ignition pixel),",
                        "`pixels`, and `fire_id`.")),
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
#' Spreads fires with `SpaDES.tools::spreadCpp()` from every pixel in
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
      successfulEscapes <- sim$ignitionsAndEscapes[escapes > 0]
      igLocs <- rep(successfulEscapes$pixelID, times = successfulEscapes$escapes)
      ## this year's spread probability, with the year's random effect if the fit has one
      spreadProbYear <- yearSpreadProb(sim$fireSense_SpreadPredicted, sim$fireSense_SpreadSD)
      ## spreadCpp(), the spread the fit uses (fireSenseUtils' objective), so a forecast spreads fires
      ## as the fitted parameters assume. Several escapes on one pixel are one fire: the first start
      ## burns the pixel and the others cannot.
      spreadState <- SpaDES.tools::spreadCpp(
        landscape = sim$fireSense_SpreadPredicted,
        loci = igLocs,
        spreadProb = spreadProbYear,
        directions = 8L)
      spreadState <- data.table(initialPixels = spreadState$initialLocus, pixels = spreadState$indices)
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

#' This year's spread probability, with the per-year random effect
#'
#' The fit (`fireSenseUtils::.objfunSpreadFit()` with `yearSpreadSD`) gives each year one eps on logit
#' spread probability, shared by all that year's fires. Here each year draws z ~ N(0, 1) and uses
#' plogis(qlogis(p) + z * sd). With a raster `sd` (several ELFs, each with its own fitted sd) the ELFs
#' share the year's z, each scaled by its own sd.
#'
#' @param spreadProbRas `SpatRaster` of spread probability.
#' @param sd `NULL`, one number, or a `SpatRaster` aligned with `spreadProbRas`.
#' @return Numeric vector of spread probability, one per pixel.
yearSpreadProb <- function(spreadProbRas, sd) {
  p <- terra::values(spreadProbRas, mat = FALSE)
  if (is.null(sd)) return(p)
  s <- if (inherits(sd, "SpatRaster")) terra::values(sd, mat = FALSE) else rep_len(sd, length(p))
  s[is.na(s)] <- 0
  if (!any(s > 0)) return(p)
  z <- stats::rnorm(1)
  stats::plogis(stats::qlogis(p) + z * s)
}
