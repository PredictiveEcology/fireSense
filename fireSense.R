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
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense.Rmd"),
  reqdPkgs = list("data.table", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings 
                            mapping the names of data objects required by the 
                            module to those in the simList environment."),
    defineParameter(name = "initialRunTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the 
                            simulation."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, 
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time.")
  ),
  inputObjects = rbind(
    expectsInput(
      objectName = "ignitProb",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer or RasterStack (time series) describing spatial
              variations in ignition probabilities."
    ),
    expectsInput(
      objectName = "escapeProb",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer or RasterStack (time series) describing spatial
              variations in escape probabilities."
    ),
    expectsInput(
      objectName = "spreadProb",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer or RasterStack describing spatial variations in the
              spread probabilities."
    ),
    expectsInput(
      objectName = "ageMap",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer describing spatial variations in the age of forest 
              stands at the start of the simulation."
    ),
    expectsInput(
      objectName = "vegMap",
      objectClass = "RasterStack",
      sourceURL = NA_character_,
      desc = "A RasterStack describing the spatial distribution of land-cover 
              classes at the start of the simulation. There should be one layer
              per land-cover class. Each layer should describe the proportion of
              each cell covered by a specific land-cover class."
    )
  ),
  outputObjects = rbind(
    createsOutput(
      objectName = "ageMap",
      objectClass = "RasterLayer",
      desc = "A RasterLayer describing spatial variations in the age of forest
              stands at the end of the simulation."
    ),
    createsOutput(
      objectName = "vegMap",
      objectClass = "RasterStack",
      desc = "A RasterStack describing the spatial distribution of land-cover
              classes at the end of the simulation. Each layer describes the
              proportion of each cell covered by a specific land-cover class."
    )
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- sim$fireSenseInit(sim)

  } else if (eventType == "plot") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event

    #Plot(objectFromModule) # uncomment this, replace with object to plot
    # schedule future event(s)

    # e.g.,
    #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "fireSense", "plot")

    # ! ----- STOP EDITING ----- ! #
  } else if (eventType == "save") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event

    # e.g., call your custom functions/methods here
    # you can define your own methods below this `doEvent` function

    # schedule future event(s)

    # e.g.,
    # sim <- scheduleEvent(sim, time(sim) + increment, "fireSense", "save")

    # ! ----- STOP EDITING ----- ! #
  } else if (eventType == "burn") {
    sim <- sim$fireSenseBurn(sim)
      
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  invisible(sim)
}


fireSenseInit <- function(sim) {
  
  sim <- scheduleEvent(sim, eventTime = P(sim)$initialRunTime, current(sim)$moduleName, "burn")
  sim
  
}

### template for save events
fireSenseSave <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
fireSensePlot <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot("object")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

fireSenseBurn <- function(sim) {
  
  moduleName <- current(sim)$moduleName
  
  ## Mapping
    iP <- sim[[
      if (tryCatch(!is.null(P(sim)$mapping[["ignitProb"]]), error = function(e) FALSE)) {
        P(sim)$mapping[["ignitProb"]]
      } else "ignitProb"
    ]]
  
    eP <- sim[[
      if (tryCatch(!is.null(P(sim)$mapping[["escapeProb"]]), error = function(e) FALSE)) {
        P(sim)$mapping[["escapeProb"]]
      } else "escapeProb"
    ]]
    
    sP <- sim[[
      if (tryCatch(!is.null(P(sim)$mapping[["spreadProb"]]), error = function(e) FALSE)) {
        P(sim)$mapping[["spreadProb"]]
      } else "spreadProb"
    ]]
    
    AM <- if (tryCatch(!is.null(P(sim)$mapping[["ageMap"]]), error = function(e) FALSE)) {
      P(sim)$mapping[["ageMap"]]
    } else "ageMap"
    
  ## Ignite
  ignited <- which(as.logical(rbinom(n = ncell(iP), size = 1, prob = iP[])))
  
  ## Escape
  loci <- ignited[eP[ignited] > runif(length(ignited))]
  
  if (length(loci) > 0L) {
  
    ## Spread
    fires <- SpaDES::spread(sP, loci = loci, spreadProb = sP, returnIndices = TRUE)
  
    ## Update age map
      if (is(sim[[AM]], "RasterLayer")) {

        sim[[AM]][fires[["indices"]]] <- 0
        
      } else if (is.data.table(sim[[AM]])) {
        
        sim[[AM]][px_id %in% fires[["indices"]], age := 0L] 
        
      }
    
    #sim$fireSize[[time(sim) - start(sim) + 1L]] <- tabulate(fires[["id"]])
  }
  
  if (!is.na(P(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + P(sim)$intervalRunModule, moduleName, "burn")
  
  sim
  
}
