# fireSense (development version)

- One fire per escaped ignition. `ignitionsAndEscapes` must have `escaped` (fireSense_IgnitionPredict >=
  1.0.0.9003); `escapes` is a coarse pixel's count repeated on each of its ignitions, and spreading `escapes`
  fires from every row gave, e.g., 8 fires where there were 2 escapes.
- New parameter `escapeSizeHa` (default 50), as in the spread fit: each escaped fire burns its first 50 ha
  whatever its spread probability (`SpaDES.tools::spreadCpp(minSize =)`), then spreads normally. `jumpTries`
  (default 0, off) and `jumpMeanDist` (default 3 cells) pass through to `spreadCpp()` for fires stuck under
  that size.
- New input `nonEscapedFireSizesHa` (from fireSense_dataPrepFit): each ignition that did not escape burns a
  small patch with a size drawn from the study area's observed fires below `escapeSizeHa`. Without it, as
  before, those ignitions burn nothing.
- Needs SpaDES.tools >= 2.1.3.9009 and fireSenseUtils >= 0.2.3.9044. Version 2.0.2.9002.

- Fires spread with `SpaDES.tools::spreadCpp()`, the spread the fit uses (fireSenseUtils' objective), instead
  of `spread2()`, so a forecast spreads fires as the fitted parameters assume. Several escapes on one pixel are one
  fire; `spread2()` stopped with "start has duplicates" at three or more. `burnDT` now holds `initialPixels`,
  `pixels` and `fire_id`.
- New input `fireSense_SpreadSD`, the per-year random effect fitted by fireSense_SpreadFit (`yearSpreadSD`): each
  year draws one z ~ N(0, 1) and every fire that year spreads with plogis(qlogis(p) + z * sd). A seasonal
  departure: a bad year makes every fire bigger. A raster sd (one per ELF) scales the shared z per ELF. `NULL`
  or 0: no effect.


## Breaking changes

- Removed parameters that nothing read: `plotIgnitions`, `.saveInitialTime`, `.saveInterval`. Stop setting them.
- Input `flammableRTM` is now declared `SpatRaster` (was `list`), which is how it was always used.

# fireSense 2.0.2

First release from `development` since `master` was last updated (2022-02-17). Full history: https://github.com/PredictiveEcology/fireSense/compare/504f990...v2.0.2

## Breaking changes

- Removed input `fireSense_EscapePredicted` (RasterLayer).
- Removed input `fireSense_IgnitionPredicted` (data.frame).
- Input `fireSense_SpreadPredicted` is now `SpatRaster` (was `RasterLayer`).
- Output `burnMap` is now `SpatRaster` (was `RasterLayer`).
- Output `rstAnnualBurnID` is now `SpatRaster` (was `RasterLayer`).
- Output `rstCurrentBurn` is now `SpatRaster` (was `RasterLayer`).
- Removed parameter: `.plotInitialTime`.

## New features

- New inputs: `flammableRTM`, `ignitionsAndEscapes`, `rasterToMatch`.
- New parameters: `.plots`, `plotIgnitions`.

## Dependencies

- No longer depends on `raster`.
- Now depends on `terra`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
