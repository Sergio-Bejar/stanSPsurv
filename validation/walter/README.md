/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
# Walter (2015) benchmark

This benchmark reproduces the data preparation and model formula documented by
`BayesSPsurv` for Barbara Walter's post-war peace data. The input has 1,237
country-years, 77 survival spells, 46 countries, and 24 recurrent-war events.
The published 800 km capital-distance graph has 27 undirected edges and 16
isolated countries.

The Stan runs use four chains, `adapt_delta = 0.97`, and
`max_treedepth = 13`. Pooled and exchangeable runs contain 500 warm-up and 500
sampling iterations per chain. The spatial diagnostic uses 250 warm-up and 250
sampling iterations per chain.

## Results

| Specification | Divergences | Max treedepth | R-hat | Assessment |
|---|---:|---:|---:|---|
| Pooled Weibull | 0/2,000 | 0 | <= 1.01 | Passed; 12 seconds |
| Exchangeable Weibull | 0/2,000 | 0 | <= 1.01 | Passed; 31 seconds |
| Spatial proper CAR, official 800 km | 0/1,000 | 0 | <= 1.01 | Passed; 841 seconds |

The original implementation produced 394/2,000 pooled divergences and
542/2,000 exchangeable divergences; its spatial smoke test stalled. The repair
uses structure-specific Stan programs, an episode-level mixture likelihood,
log-scale and non-centered frailty parameters, regularizing priors, centered
Walter predictors, and a component-aware proper-CAR mixture. These changes remove the
observed geometric failures. The spatial implementation is still computationally
expensive and its substantive estimates must be compared with the original
`BayesSPsurv` implementation before publication.
