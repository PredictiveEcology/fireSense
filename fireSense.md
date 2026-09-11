---
title: "fireSense Manual"
subtitle: "v.0.0.0.9000"
date: "Last updated: 2025-04-08"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
  bibliography: citations/references_fireSense.bib
link-citations: true
always_allow_html: true
---

# fireSense Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:fireSense) *fireSense*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Jean Marchal <jean.d.marchal@gmail.com> [aut, cre], Alex M. Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

A landscape fire model sensitive to environmental changes (e.g., weather and land cover).

- [@Marchal:2017a; @Marchal:2017b; @Marchal:2019]

### Module summary

Provide a brief summary of what the module does / how to use the module.

Module documentation should be written so that others can use your module.
This is a template for module documentation, and should be changed to reflect your module.

### Module inputs and parameters

Describe input data required by the module and how to obtain it (e.g., directly from online sources or supplied by other modules)
If `sourceURL` is specified, `downloadData("fireSense", "..")` may be sufficient.
Table \@ref(tab:moduleInputs-fireSense) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense)(\#tab:moduleInputs-fireSense)List of (ref:fireSense) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireSense_IgnitionPredicted </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A SpatRaster of ignition probabilities. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_EscapePredicted </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A SpatRaster of escape probabilities. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_SpreadPredicted </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A SpatRaster of spread probabilities. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Provide a summary of user-visible parameters (Table \@ref(tab:moduleParams-fireSense))


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-fireSense)(\#tab:moduleParams-fireSense)List of (ref:fireSense) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> plotIgnitions </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> whether to plot ignitions, escapes, and burns </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. When to start plotting. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between plot events. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> time to simulate initial fire </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between two runs of this module, expressed in units of simulation time. By default, 1 year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. When to start saving output to a file. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between save events. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> whichModulesToPrepare </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> fireSens.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Which fireSense predict modules to prep? Defaults to all 3. Must include `fireSense_IgnitionPredict`. </td>
  </tr>
</tbody>
</table>

### Events

<!-- TODO update for latest version -->
- ignite fires
- determine which fires escape
- spread escaped fires
- save
- plot

### Plotting

<!-- TODO update for latest version -->
- **Burn map**: Pixels burned this timestep.
- **Cumulative burn map**: Number of times each pixel burned during the simulation..

### Saving

<!-- TODO update for latest version -->
- `burnMap`: A `RasterLayer` describing how which pixels burned this timestep.
- `burnMapCumul`: A `RasterLayer` describing how many times each pixel burned over the course of the simulation.

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-fireSense)(\#tab:moduleOutputs-fireSense)List of (ref:fireSense) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> burnDT </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Data table with pixel IDs of most recent burn. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> burnMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A raster of cumulative burns </td>
  </tr>
  <tr>
   <td style="text-align:left;"> burnSummary </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Describes details of all burned pixels. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstAnnualBurnID </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> annual raster whose values distinguish individual fires </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstCurrentBurn </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A binary raster with 1 values representing burned pixels. </td>
  </tr>
</tbody>
</table>

### Usage

<!-- TODO: update for latest version; use terra -->


``` r
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

### Links to other modules

<!-- TODO: describe links with the other fireSense modules -->
This module should be coupled with a dynamic vegetation model.

### Getting help

- <https://github.com/PredictiveEcology/fireSense/issues>
