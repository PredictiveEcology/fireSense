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
