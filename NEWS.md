# fireSense 1.1.0

This is the first release of fireSense as one family. `fireSense@main` loads the released version of every fireSense module, and `fireSense@development` their development versions, so a project no longer mixes released and in-progress modules by accident. To load exactly this set of releases, and nothing newer, use `PredictiveEcology/fireSense@fireSense-1.1.0`: every module carries that tag at its own release (ELFs 1.2.0, dataPrepFit 1.3.0, ignitionFit 1.2.0, spreadFit 1.2.0, dataPrepPredict 1.1.0, ignitionPredict 1.2.0, spreadPredict 1.2.0, burn 2.2.0).

The summary module is the exception for now: it still comes from its `modsForFireSense` branch. The parent module is licensed under GPL-3.

- License: GPL-3.
- Children are listed by name, so they follow the parent's account and branch (`fireSense@development` -> each child's `development`, `fireSense@main` -> each child's `main`); `fireSense_summary` stays on `modsForFireSense`.

# fireSense 1.0.0

- New parent module. Listing `PredictiveEcology/fireSense@development` in `setupProject(modules = )` fetches and runs the nine FireSense child modules.
- The name `fireSense` used to be the burn module; that module is now `fireSense_burn`.
- For the old-to-new module names and how to migrate, see the table in README.md.
