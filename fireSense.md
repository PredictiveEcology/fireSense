---
title: "fireSense"
author: "Jean Marchal"
date: "June 2022"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---



# Overview

A landscape fire model sensitive to environmental changes (e.g. weather and land-cover).

# Usage


```r
Require(c("SpaDES.core", "SpaDES.tools"))

set.seed(1)

nx <- ny <- 100L
n <- nx * ny
r <- raster(nrows = ny, ncols = nx, xmn = -nx/2, xmx = nx/2, ymn = -ny/2, ymx = ny/2)

# Create a map ignition probabilities
ignitionProbRaster <- gaussMap(r, scale = 10, var = .0001, speedup = nx/5e2, inMemory = TRUE)

# Create a map of escape probabilities
escapeProbRaster <- gaussMap(r, scale = 50, var = .01, speedup = nx/5e2, inMemory = TRUE)

# Create a map of spread probabilities
spreadProbRaster <- gaussMap(r, scale = 300, var = .05, speedup = nx/5e2, inMemory = TRUE)

#outputDir <- file.path(tempdir(), "outputs")
times <- list(start = 1, end = 100, timeunit = "year")

modules <- list("fireSense")

# Pass objects found in the global environment to the simList environment
objects <- list(
  ignitionProbRaster = ignitionProbRaster,
  escapeProbRaster = escapeProbRaster,
  spreadProbRaster = spreadProbRaster
)

paths <- list(
  # cachePath = file.path(outputDir, "cache"),
  modulePath = ".."
  # inputPath = inputDir,
  # outputPath = outputDir
)

mySim <- simInit(times = times, params = parameters, modules = modules, objects = objects, paths = paths)

spades(mySim)
```

# Parameters


```
## defineParameter: '.runInitialTime' is not of specified type 'numeric'.
```



|paramName             |paramClass |default      |min |max |paramDesc                                                                                                      |
|:---------------------|:----------|:------------|:---|:---|:--------------------------------------------------------------------------------------------------------------|
|.runInitialTime       |numeric    |start(sim)   |NA  |NA  |time to simulate initial fire                                                                                  |
|whichModulesToPrepare |character  |fireSens.... |NA  |NA  |Which fireSense fit modules to prep? defaults to all 3. Must include ignition                                  |
|.runInterval          |numeric    |1            |NA  |NA  |optional. Interval between two runs of this module, expressed in units of simulation time. By default, 1 year. |
|.saveInitialTime      |numeric    |NA           |NA  |NA  |optional. When to start saving output to a file.                                                               |
|.saveInterval         |numeric    |NA           |NA  |NA  |optional. Interval between save events.                                                                        |
|.plotInitialTime      |numeric    |NA           |NA  |NA  |optional. When to start plotting.                                                                              |
|.plotInterval         |numeric    |NA           |NA  |NA  |optional. Interval between plot events.                                                                        |

# Events

- ignite fires
- determine which fires escape
- spread escaped fires
- save
- plot

## Plotting

- **Burn map**: Pixels burned this timestep.
- **Cumulative burn map**: Number of times each pixel burned during the simulation..

## Saving

- **burnMap**: A RasterLayer describing how which pixels burned this timestep.
- **burnMapCumul**: A RasterLayer describing how many times each pixel burned over the course of the simulation.


# Data dependencies


```
## defineParameter: '.runInitialTime' is not of specified type 'numeric'.
```



|objectName                  |objectClass |desc                                     |sourceURL |
|:---------------------------|:-----------|:----------------------------------------|:---------|
|fireSense_IgnitionPredicted |data.frame  |A RasterLayer of ignition probabilities. |NA        |
|fireSense_EscapePredicted   |RasterLayer |A RasterLayer of escape probabilities.   |NA        |
|fireSense_SpreadPredicted   |RasterLayer |A RasterLayer of spread probabilities.   |NA        |

## Output data

Description of the module outputs.


```
## defineParameter: '.runInitialTime' is not of specified type 'numeric'.
```



|objectName      |objectClass |desc                                                      |
|:---------------|:-----------|:---------------------------------------------------------|
|rstCurrentBurn  |RasterLayer |A binary raster with 1 values representing burned pixels. |
|rstAnnualBurnID |RasterLayer |annual raster whose values distinguish individual fires   |
|burnMap         |RasterLayer |A raster of cumulative burns                              |
|burnDT          |data.table  |Data table with pixel IDs of most recent burn.            |
|burnSummary     |data.table  |Describes details of all burned pixels.                   |

# Links to other modules

This module should be coupled with a dynamic vegetation model.

