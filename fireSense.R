# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense",
  description = "A landscape fire model, sensitive to environmental changes (e.g.
                 weather and land-cover).",
  keywords = c("fire", "percolation", "environmental control", "feedback",
               "weather", "vegetation", "land-cover"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense.Rmd"),
  reqdPkgs = list("data.table", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc ="time to simulate initial fire"),
    defineParameter(name = "whichModulesToPrepare", class = "character",
                    default = c("fireSense_SpreadPredict", "fireSense_IgnitionPredict", "fireSense_EscapePredict"),
                    NA, NA, desc = "Which fireSense fit modules to prep? defaults to all 3. Must include ignition"),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time. By default, 1 year."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA,
                    desc = "optional. Interval between save events."),
    defineParameter(name = "plotIgnitions", class = 'logical', FALSE, NA, NA,
                    desc = "plot ignitions and escapes as SpatialPoints"),
    defineParameter(name = ".plotInitialTime", class = "numeric", default = NA,
                    desc = "optional. When to start plotting."),
    defineParameter(name = ".plotInterval", class = "numeric", default = NA,
                    desc = "optional. Interval between plot events.")
  ),
  inputObjects = rbind(
    expectsInput(objectName = "fireSense_IgnitionPredicted", objectClass = "data.frame",
                 desc = "A RasterLayer of ignition probabilities."),
    expectsInput(objectName = "fireSense_EscapePredicted", objectClass = "RasterLayer",
                 desc = "A RasterLayer of escape probabilities."),
    expectsInput(objectName = "fireSense_SpreadPredicted", objectClass = "RasterLayer",
                 desc = "A RasterLayer of spread probabilities.")
  ),
  outputObjects = rbind(
    createsOutput(objectName = "rstCurrentBurn", objectClass = "RasterLayer",
                  desc = "A binary raster with 1 values representing burned pixels."),
    createsOutput(objectName = "rstAnnualBurnID", objectClass = "RasterLayer",
                  desc = "annual raster whose values distinguish individual fires"),
    createsOutput(objectName = "burnMap", objectClass = "RasterLayer",
                  desc = "A raster of cumulative burns"),
    createsOutput(objectName = "burnDT", objectClass = "data.table",
                  desc = "Data table with pixel IDs of most recent burn."),
    createsOutput(objectName = "burnSummary", objectClass = "data.table", desc = "Describes details of all burned pixels.")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {

      #trying to avoid the raster warning no non-missing arguments to max
      sim$burnMap <- sim$flammableRTM
      sim$burnMap[getValues(sim$burnMap) == 0] <- NA #make a map of flammable pixels with value 0
      sim$burnMap[getValues(sim$burnMap) == 1] <- 0

      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "burn", eventPriority = 5.13)

      if (!is.na(P(sim)$.plotInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, moduleName, "plot", eventPriority = .last())
    },
    burn = {
      sim <- burn(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "burn", eventPriority = 5.13)
    },
    plot = {
      sim <- plot(sim)

      if (!is.na(P(sim)$.plotInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, moduleName, "plot", eventPriority = .last())
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  invisible(sim)
}

burn <- function(sim) {

  moduleName <- current(sim)$moduleName

  ## Ignite
  notNA <- which(!is.na(sim$fireSense_IgnitionPredicted[]))
  ignitionProbs <- sim$fireSense_IgnitionPredicted[notNA]

  ignited <- notNA[which(
    rbinom(n = length(ignitionProbs),
           size = 1,
           prob = pmin(ignitionProbs, 1)
    ) > 0
  )]

  ignited <- sample(ignited) # Randomize order

  rm(ignitionProbs)
  if (P(sim)$plotIgnitions) {
    mod$ignitions <- ignited
  }

  if (length(ignited) > 0L) {
    if ("fireSense_EscapePredict" %in% P(sim)$whichModulesToPrepare) {
      ## Escape
      adjacent <- SpaDES.tools::adj(
        x = sim$fireSense_EscapePredicted,
        cells = ignited,
        directions = 8,
        returnDT = TRUE
      )

      if (is.matrix(adjacent))
        adjacent <- as.data.table(adjacent)

      from <- unique(adjacent, by = "from")
      from[, `:=` (probEscape = sim$fireSense_EscapePredicted[from], to = NULL)]

      # Update probEscape to get p0
      from <- from[!is.na(probEscape),]
      p0 <- with(
        data = from[adjacent, on = "from"][, probEscape := (1 - (1 - probEscape)^(1 / .N)), by = "from"],
        expr = {
          p0 <- sim$fireSense_EscapePredicted
          p0[to] <- probEscape
          p0
        }
      )

      mod$spreadState <- SpaDES.tools::spread2(
        landscape = sim$fireSense_EscapePredicted,
        start = ignited,
        iterations = 1,
        spreadProb = p0,
        directions = 8L,
        asRaster = FALSE
      )

      if (P(sim)$plotIgnitions) {
        mod$escapes <- unique(mod$spreadState[state == "activeSource",]$initialPixels)
      }
    }
    if ("fireSense_SpreadPredict" %in% P(sim)$whichModulesToPrepare) {
      ## Spread
      # Note: if none of the cells are active SpaDES.tools::spread2() returns spreadState unchanged
      mod$spreadState <- SpaDES.tools::spread2(
        landscape = sim$fireSense_SpreadPredicted,
        spreadProb = sim$fireSense_SpreadPredicted,
        directions = 8L,
        start = mod$spreadState,
        asRaster = FALSE)

      mod$spreadState[ , fire_id := .GRP, by = "initialPixels"] # Add an fire_id column

      sim$rstAnnualBurnID <- raster(sim$fireSense_SpreadPredicted)
      sim$rstCurrentBurn <- raster(sim$fireSense_SpreadPredicted)

      sim$rstAnnualBurnID[mod$spreadState$pixels] <- mod$spreadState$fire_id
      sim$rstCurrentBurn[mod$spreadState$pixels] <- 1
      sim$burnMap[mod$spreadState$pixels] <- sim$burnMap[mod$spreadState$pixels] + 1

      #get fire year, pixels burned, area burned, poly ID of all burned pixels
      # Make burnSummary --> similar to SCFM
      sim$burnDT <- mod$spreadState

      tempDT <- sim$burnDT[, .(.N), by = "initialPixels"]
      tempDT$year <- time(sim)
      tempDT$areaBurnedHa <- tempDT$N * prod(res(sim$fireSense_SpreadPredicted))*1e-4
      setnames(tempDT, c("initialPixels"), c("igLoc"))
      sim$burnSummary <- rbind(sim$burnSummary, tempDT)
    }
  }

  invisible(sim)
}

plot <- function(sim) {

  if (P(sim)$plotIgnitions){
   #TODO: make this figure a ggplot instead

    #for potting
    ignitions <- mod$ignitions[!mod$ignitions %in% mod$escapes]

    ignitions <- raster::xyFromCell(sim$flammableRTM, cell = mod$ignitions)
    escapes <- raster::xyFromCell(sim$flammableRTM, cell = mod$escapes)
    ignitions <- SpatialPointsDataFrame(coords = ignitions,
                                        proj4string = crs(sim$flammableRTM),
                                        data = data.frame(color = rep("yellow",
                                                                      times = nrow(ignitions))))

    escapes <- SpatialPointsDataFrame(coords = escapes,
                                      proj4string = crs(sim$flammableRTM),
                                      data = data.frame(color = rep("red",
                                                                    times = nrow(escapes))))
    flamMap <- sim$flammableRTM
    IandE <- rbind(ignitions, escapes)

    clearPlot()
    mapTitle <- paste0("ignitions and escapes in ", time(sim))
    raster::plot(flamMap, col = c("#8DA0CB", "#66C2A5"), main = mapTitle)
    raster::plot(IandE, add = TRUE, col = IandE$color, pch = 19, cex = 0.5)

    mod$ignitions <- NULL
    mod$escapes <- NULL
  }
  mapTitle <- paste0("Burn map ", time(sim))
  clearPlot()
  Plot(sim$rstCurrentBurn, new = TRUE,  title = mapTitle)
  Plot(sim$burnMap, title = "Cumulative burn map")

  invisible(sim)
}
