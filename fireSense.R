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
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense.Rmd"),
  reqdPkgs = list("data.table", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
      desc = "optional. Named character vector or list mapping the names of data
              objects required by the module to those in the simList environment."),
    defineParameter(name = "initialRunTime", class = "numeric", default = start(sim),
      desc = "optional. Simulation time at which to start this module. Defaults 
              to simulation start time."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, 
      desc = "optional. Interval in simulation time units between two runs of 
              this module.")
  ),
  inputObjects = data.frame(
    objectName = c("ignitProb", "escapeProb", "spreadProb", "ageMap"),
    objectClass = c("raster", "raster", "raster", "raster, data.table"),
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = NA_character_,
    objectClass = NA_character_,
    other = NA_character_,
    stringsAsFactors = FALSE
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
    #sim <- scheduleEvent(sim, p(sim)$.plotInitialTime, "fireSense", "plot")

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
  
  sim <- scheduleEvent(sim, eventTime = p(sim)$initialRunTime, "fireSense", "burn")
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
  
  ## Mapping
    iP <- sim[[
      if (tryCatch(!is.null(p(sim)$mapping[["ignitProb"]]), error = function(e) FALSE)) {
        p(sim)$mapping[["ignitProb"]]
      } else "ignitProb"
    ]]
  
    eP <- sim[[
      if (tryCatch(!is.null(p(sim)$mapping[["escapeProb"]]), error = function(e) FALSE)) {
        p(sim)$mapping[["escapeProb"]]
      } else "escapeProb"
    ]]
    
    sP <- sim[[
      if (tryCatch(!is.null(p(sim)$mapping[["spreadProb"]]), error = function(e) FALSE)) {
        p(sim)$mapping[["spreadProb"]]
      } else "spreadProb"
    ]]
    
    AM <- if (tryCatch(!is.null(p(sim)$mapping[["ageMap"]]), error = function(e) FALSE)) {
      p(sim)$mapping[["ageMap"]]
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
  
  if (!is.na(p(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + p(sim)$intervalRunModule, "fireSense", "burn")
  
  sim
  
}
