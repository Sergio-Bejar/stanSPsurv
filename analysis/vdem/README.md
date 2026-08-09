# V-Dem application: democratic breakdown

This application studies the duration of democratic spells and the risk of a
transition to autocracy. It uses V-Dem version 16 together with the Episodes of
Regime Transformation (ERT) version 16 release.

The unit of analysis used for estimation is a **democratic spell**, not a
country-year. A spell begins when `v2x_regime` enters electoral or liberal
democracy (categories 2 or 3) and ends in the last democratic year before a
transition to electoral or closed autocracy (categories 1 or 0). Spells still
democratic in the last observed year are right-censored. This avoids treating
the multiple annual rows of one spell as independent latent split-population
draws.

Covariates are measured at the start of each spell. The main specification
uses judicial constraints, legislative constraints, freedom of expression,
log GDP per capita, civil war, and prior democratic experience. The first five
are lagged one year whenever a pre-spell observation exists. ERT variables are
merged for descriptive validation; ERT democratization episodes are not
mistaken for democratic survival spells.

## Reproduce

From the repository root:

```sh
Rscript analysis/vdem/prepare_vdem.R
Rscript analysis/vdem/fit_models.R pooled weibull
Rscript analysis/vdem/fit_models.R exchangeable weibull
Rscript analysis/vdem/fit_models.R spatial weibull
```

`prepare_vdem.R` accepts an installed/cloned ERT source through the
`ERT_SOURCE` environment variable. Otherwise it downloads the official ERT
v16 source archive from GitHub. Raw inputs and fitted model objects are ignored
by Git; the compact diagnostics in `output/summary.csv` are reproducible.

The repository includes the compact output of the first full pooled Weibull
run (`output/pooled_weibull_summary.csv` and
`output/pooled_weibull_diagnostics.csv`). All four chains had R-hat at 1.00,
zero divergences, zero maximum-treedepth transitions, and adequate E-BFMI.

The spatial graph is a symmetric four-nearest-neighbour graph based on
country-polygon centroids from the `maps` package. It is a reproducible
geographic-proximity graph, not a historical border-contiguity matrix. That
choice keeps islands in the CAR graph and should be replaced by time-varying
historical borders for a definitive substantive analysis.

## Required packages

- preparation: `survival`, `countrycode`, `maps`
- estimation: `stanSPsurv`, `cmdstanr`

V-Dem and ERT have their own citation and data-use requirements. Consult the
official release files before distributing derived data.
