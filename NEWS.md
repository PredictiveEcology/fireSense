# fireSense (development version)

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
