# stanSPsurv 0.0.1.9000

- Created the initial R package structure.
- Added candidate pooled, exchangeable, and spatial split-population survival models.
- Added Weibull and log-logistic baselines through a shared Stan backend.
- Added coefficient, prediction, and frailty helpers.
- Added initial input-validation tests and package branding.
- Compiled the shared backend with CmdStan 2.39.0 and ran all six candidate
  model combinations on simulated data.
- Fixed an improper auxiliary parameter block in the pooled specification.
- Added a reproducible validation script and recorded convergence diagnostics.
- Achieved `R CMD check --no-manual` status OK on R 4.2.1 for macOS.
