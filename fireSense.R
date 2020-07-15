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
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings 
                            mapping one or more inputs required by the module to
                            objects loaded in the simList environment."),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start
                            time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1, 
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time. By default, 1 year."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA, 
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA, 
                    desc = "optional. Interval between save events."),
    defineParameter(name = ".plotInitialTime", class = "numeric", default = NA, 
                    desc = "optional. When to start plotting."),
    defineParameter(name = ".plotInterval", class = "numeric", default = NA, 
                    desc = "optional. Interval between plot events.")
  ),
  inputObjects = rbind(
    expectsInput(
      objectName = "ignitionProbRaster",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer or RasterStack (time series) describing spatial
              variations in ignition probabilities."
    ),
    expectsInput(
      objectName = "escapeProbRaster",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer or RasterStack (time series) describing spatial
              variations in escape probabilities."
    ),
    expectsInput(
      objectName = "spreadProbRaster",
      objectClass = "RasterLayer",
      sourceURL = NA_character_,
      desc = "A RasterLayer or RasterStack describing spatial variations in the
              spread probabilities."
    )
  ),
  outputObjects = rbind(
    createsOutput(
      objectName = "rstCurrentBurn",
      objectClass = "RasterLayer",
      desc = "A RasterLayer describing how which pixels burned this timestep."
    ),
    createsOutput(
      objectName = "burnMap",
      objectClass = "RasterLayer",
      desc = "A RasterLayer describing how many times each pixel burned over the
              course of the simulation."
    )
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense = function(sim, eventTime, eventType, debug = FALSE) 
{
  moduleName <- current(sim)$moduleName
  
  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "burn")
      
      if (!is.na(P(sim)$.plotInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "plot", .last())
    },
    burn = { 
      sim <- burn(sim)
      
      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "burn")
    },
    plot = { 
      sim <- plot(sim)
      
      if (!is.na(P(sim)$.plotInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, moduleName, "plot")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  invisible(sim)
}

burn <- function(sim) 
{
  moduleName <- current(sim)$moduleName
  ## Mapping
  mod[["ignitionProbRaster"]] <- if (!is.null(P(sim)[["mapping"]][["ignitionProbRaster"]])){
    sim[[P(sim)[["mapping"]][["ignitionProbRaster"]]]]
  } else {
    if (!is.null(sim[["ignitionProbRaster"]])){
      sim[["ignitionProbRaster"]]
    } else {
      if (!is.null(sim[["fireSense_FrequencyPredicted"]])){
        sim[["fireSense_FrequencyPredicted"]]
      } else {
        if (!is.null(sim[["fireSense_IgnitionPredicted"]])){
          sim[["fireSense_IgnitionPredicted"]]
        } else {
          stop("Neither `fireSense_FrequencyPredicted` (i.e. being deprecated), `fireSense_IgnitionPredicted` nor 
                 `ignitionProbRaster` were found. Please provide one of these")
        }
      }        
    }
  }
  
  mod[["escapeProbRaster"]] <- if (!is.null(P(sim)[["mapping"]][["escapeProbRaster"]])){
    sim[[P(sim)[["mapping"]][["escapeProbRaster"]]]]
  } else {
    if (!is.null(sim[["escapeProbRaster"]])){
      sim[["escapeProbRaster"]]
    } else {
      if (!is.null(sim[["fireSense_EscapePredicted"]])){
        sim[["fireSense_EscapePredicted"]]
      } else {
        stop("Neither `fireSense_EscapePredicted` nor `escapeProb` were found. Please provide one of these")
      }        
    }
  }
  
  mod[["spreadProbRaster"]] <- if (!is.null(P(sim)[["mapping"]][["spreadProbRaster"]])) {
    sim[[P(sim)[["mapping"]][["spreadProbRaster"]]]]
  } else {
    if (!is.null(sim[["spreadProbRaster"]])){
      sim[["spreadProbRaster"]]
    } else {
      if (!is.null(sim[["fireSense_SpreadPredicted"]])){
        sim[["fireSense_SpreadPredicted"]]
      } else {
        stop("Neither `fireSense_SpreadPredicted` nor `spreadProb` were found. Please provide one of these")
      }        
    }
  }
  
  if (is.null(sim$burnMap))
  {
    sim$burnMap <- mod[["spreadProbRaster"]]
    sim$burnMap[!is.na(sim$burnMap[])] <- 0
  }
  
  ## Ignite
  notNA <- which(!is.na(mod[["ignitionProbRaster"]][]))
  ignitionProbs <- mod[["ignitionProbRaster"]][notNA]

  ignited <- notNA[which(
    rbinom(n = length(ignitionProbs),
           size = 1,
           prob = pmin(ignitionProbs, 1)
    ) > 0
  )]
  
  ignited <- sample(ignited) # Randomize order
  
  rm(ignitionProbs)
  
  if (length(ignited) > 0L)
  {
    ## Escape
    adjacent <- SpaDES.tools::adj(
      x = mod[["escapeProbRaster"]],
      cells = ignited,
      directions = 8,
      returnDT = TRUE
    )
    
    if (is.matrix(adjacent))
      adjacent <- as.data.table(adjacent)
    
    from <- unique(adjacent, by = "from")
    from[, `:=` (probEscape = mod[["escapeProbRaster"]][from], to = NULL)]
    
    # Update probEscape to get p0
    p0 <- with(
      from[adjacent, on = "from"][!is.na(probEscape)][
        , 
        probEscape := (1 - (1 - probEscape)^(1 / .N)),
        by = "from"
        ],
      {
        p0 <- mod[["escapeProbRaster"]]
        p0[to] <- probEscape
        p0
      }
    )
    
    mod$spreadState <- SpaDES.tools::spread2(
      landscape = mod[["escapeProbRaster"]],
      start = ignited,
      iterations = 1,
      spreadProb = p0,
      directions = 8L, 
      asRaster = FALSE
    )
    
    ## Spread
    # Note: if none of the cells are active SpaDES.tools::spread2() returns spreadState unchanged
    mod$spreadState <- SpaDES.tools::spread2(
      landscape = mod[["spreadProbRaster"]],
      spreadProb = mod[["spreadProbRaster"]],
      directions = 8L,
      start = mod$spreadState,
      asRaster = FALSE
    )
    
    mod$spreadState[ , fire_id := .GRP, by = "initialPixels"] # Add an fire_id column
    
    sim$rstCurrentBurn <- raster(mod[["spreadProbRaster"]])
    sim$rstCurrentBurn[mod$spreadState$pixels] <- mod$spreadState$fire_id
    sim$rstCurrentBurn <- sim$rstCurrentBurn
    
    sim$burnMap[mod$spreadState$pixels] <- sim$burnMap[mod$spreadState$pixels] + 1
  }
  
  invisible(sim)
}

plot <- function(sim) 
{
  Plot(sim$rstCurrentBurn, sim$burnMap, title = c("Burn map", "Cumulative burn map"))
  
  invisible(sim)
}
