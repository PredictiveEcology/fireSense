### Specify module (and dependencies) definitions:
###
### name:         fireSense
###
### description:  a percolation-based fire spread model that is sensitive to climate and vegetation changes.
###               The number of fires to be spread are sampled from a user defined distribution.
###               The spread probabilities are derived from maps of climate and vegetation cover
###               #- fires may be classified in various ways (e.g. lightning vs man-caused) or considered without distinctions
###
### keywords:     fire; climate sensitive, vegetation sensitive; percolation; spread algorithm
###
### authors:      Jean Marchal <jean.d.marchal@gmail.com>
###
### version:      0.1.0
###
### spatialExtent: NA
###
### timeframe:    NA
###
### timestep:     31557600 (1 year)
###
### citation:     NA
###
### reqdPkgs:     gsl; raster
###
### parameters:   paramName: timestep
###               paramClass: numeric
###               default: 1
###
###               paramName: startTime
###               paramClass: numeric
###               default: 1
###
###               paramName: scalingFireFreq
###               paramClass: numeric
###               default: 1.0
###
###               paramName: persistence
###               paramClass: numeric
###               default: 0.00
###
###               paramName: .plotInitialTime
###               paramClass: numeric
###               default: 0
###
###               paramName: .plotInterval
###               paramClass: numeric
###               default: 1
###
###               paramName: .saveInitialTime
###               paramClass: numeric
###               default: NA
###
###               paramName: .saveInterval
###               paramClass: numeric
###               default: NA
###
### inputObjects: objectName: vegMap
###               objectClass: RasterLayer
###               other: NA
###
###               objectName: fireWeatherMap
###               objectClass: RasterLayer
###               other: NA
###
###               objectName: freqModelObj
###               objectClass: NA
###               other: NA
###
### outputObjects: objectName: firesMap
###                objectClass: RasterLayer
###                other: NA
###
###                objectName: firesCumulMap
###                objectClass: RasterLayer
###                other= NA
###
### fireSense module metadata
defineModule(sim, list(
  name="fireSense",
  description="insert module description here",
  keywords=c("fire", "fire behaviour", "climate sensitive", "vegetation sensitive", "percolation", "spread"),
  authors=person("Jean", "Marchal", email="jean.d.marchal@gmail.com", role=c("aut", "cre")),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c("2000-01-01", NA)),
  timestep=31557600,
  citation=list(),
  reqdPkgs=list("gsl", "raster"),
  parameters=rbind(
    defineParameter("timestep", "numeric", 1.0),
    defineParameter("startTime", "numeric", 1.0),
    defineParameter("scalingFireFreq", "numeric", 1.0),
    defineParameter("persistence", "numeric", 0.0),
    defineParameter(".plotInitialTime", "numeric", 0.0),
    defineParameter(".plotInterval", "numeric", 0.0)),
  inputObjects=data.frame(objectName="FireWeather",
                          objectClass="RasterLayer", 
                          other=NA_character_,stringsAsFactors=FALSE),
  outputObjects=data.frame(objectName=c("firesSizes", "firesMap","firesCumulMap"),
                           objectClass=c("integer", "RasterLayer","RasterLayer"),
                           other=rep_len(NA_character_,3),stringsAsFactors=FALSE)
))

## Toolbox: internal set of functions needed by the module
  ## Define function to compute the mean of the tapered Pareto
  meanTP <- function(lambda, theta, a){
    a + a ^ lambda * theta ^ (1L - lambda) * exp(a/theta) * gsl::gamma_inc(1L - lambda, a/theta) ## gamma_inc computes the incomplete gamma function
  }

  ## Predict method for a fireFreq object
  predict.fireFrequencyOptim <- function(object, newdata){
    mm.cl <- quote(model.matrix(object$formula, data = newdata))
    if(missing(newdata) || is.null(newdata)){
      mm.cl$data <- quote(cbind.data.frame(rbind(object$knots), object$data))
    }
    
    fit <- drop(with(as.list(object$knots), eval(mm.cl) %*% object$coefficients))
    
    object$family$linkinv(fit)
  }


  ## Predict fireSize object

### template event
doEvent.fireSense = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)
    #checkObject(name="vegMap")
    #checkObject(name="fireWeatherMap")
    
    # do stuff for this event
    sim <- fireSenseInit(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, simParams(sim)$fireSense$.plotInitialTime, "fireSense", "plot")
    #sim <- scheduleEvent(sim, simParams(sim)$fireSense$.saveInitialTime, "fireSense", "save")
  } else if (eventType=="burn") {
    # do stuff for this event
    sim <- fireSenseBurn(sim)
    
    # schedule future event(s)
    sim <- scheduleEvent(sim, simCurrentTime(sim) + simParams(sim)$fireSense$timestep, "fireSense", "burn")
  } else if (eventType=="plot") {
    # do stuff for this event
    
    # schedule future event(s)
    sim <- scheduleEvent(sim, simCurrentTime(sim) + simParams(sim)$fireSense$.plotInterval, "fireSense", "plot")
  } else {
      warning(paste("Undefined event type: '", simEvents(sim)[1L, "eventType", with=FALSE],
                    "' in module '", simEvents(sim)[1L, "moduleName", with=FALSE], "'", sep=""))
    }
  return(invisible(sim))
}

### template initilization
fireSenseInit = function(sim) {
  ### create burn map that tracks fire locations over time
  firesMap <- setValues(raster(getGlobal("fireSizes")),0L)
  
#   setColors(Fires,n=simParams(sim)$fireSense$nFires+1) <-
#     c("#FFFFFF", rev(heat.colors(simParams(sim)$fireSense$nFires)))
  assignGlobal("firesMap")
  assignGlobal("firesCumulMap", firesMap) ### Cumulative burn map
  assignGlobal("firesSizes",vector("integer",0))
  return(invisible(sim))
}

fireSenseBurn = function(sim) {
  ## Number of fires to be spread
  meanFF <- predict.fireFrequencyOptim(simParams(sim)$fireSense$freqModelObj)#, newdata = ) ## need additional parameters to get the spatial resolution used to fit this model
  
  meanFF <- meanFF * simParams(sim)$fireSense$scalingFireFreq ## Rescale the mean of the fitted distribution to the data resolution
  
  if(simParams(sim)$fireSense$freqModelObj$family$family == "poisson"){
    #resFireFreqObj[1L] * resFireFreqObj[2L] #rescale
    firesLoc <- as.logical(rpois(ncell(getGlobal("firesMap")), lambda = meanFF))
  } else {
    theta <- 
      if(!exists("theta", simParams(sim)$fireSense$freqModelObj)){
        as.numeric(gsub("[\\(\\)]", "", regmatches(
          simParams(sim)$fireSense$freqModelObj$family$family,
          gregexpr("\\(.*?\\)", simParams(sim)$fireSense$freqModelObj$family$family))[[1]]))  
      } else {
        simParams(sim)$fireSense$freqModelObj$theta
      }
    
    firesLoc <- as.logical(rnbinom(ncell(getGlobal("firesMap")), size = theta, mu = meanFF))
  }
    
  nFires <- sum(firesLoc)
  
  ## Compute Pspreads given weather and vegetation conditions
  #simParams(sim)$fireSense$sizeModel
  spreadProbMap <- 0.08 + log(getGlobal(paste0("fireWeatherMap",simCurrentTime()))) * 0.02
  
  firesMap <- spread(spreadProbMap,
                  loci=(1L:ncell(getGlobal("firesMap")))[firesLoc],
                  spreadProb=spreadProbMap,
                  persistence=simParams(sim)$fireSense$persistence,
                  mask=NULL,
                  plot.it=FALSE,
                  mapID=TRUE)
  
  firesSizes <- c(getGlobal("firesSizes"), tabulate(firesMap[]))
  
  #firesMap[is.na(getGlobal("vegMap"))] <- NA
  names(firesMap) <- "firesMap"
  setColors(firesMap,n=nFires + 1L) <- c("#FFFFFF", rev(heat.colors(nFires)))
  
  firesCumulMap <- getGlobal("firesCumulMap") + as.logical(firesMap)
  setColors(firesCumulMap) <- colorRampPalette(c("orange", "darkred"))(maxValue(firesCumulMap)+1)#c("#FFFFFF", colorRampPalette(c("orange", "darkred"))(maxValue(firesCumulMap)+1))
  
  assignGlobal("firesMap")
  assignGlobal("firesCumulMap")
  assignGlobal("firesSizes")
  
  return(invisible(sim))
}

### template for save events
fireSenseSave = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  saveFiles(sim)

  # schedule future event(s)
  sim <- scheduleEvent(sim, simCurrentTime(sim) + simParams(sim)$fireSense$.saveInterval, "fireSense", "save")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
fireSensePlot = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  Plot(firesMap, new = TRUE)
  Plot(firesCumulMap)

  # schedule future event(s)
  sim <- scheduleEvent(sim, simCurrentTime(sim) + simParams(sim)$fireSense$.plotInterval, "fireSense", "plot")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
fireSenseEvent1 = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event

  # schedule future event(s)
  sim <- scheduleEvent(sim, simCurrentTime(sim), "fireSense", "event1")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
fireSenseEvent2 = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event

  # schedule future event(s)
  sim <- scheduleEvent(sim, simCurrentTime(sim), "fireSense", "event2")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
