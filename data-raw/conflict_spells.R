#!/usr/bin/env Rscript
# Generate the conflict_spells example dataset

set.seed(20261102)

N <- 300; J <- 20
countries <- paste0("Country_", LETTERS[1:J])
regions <- rep(paste0("Region_", 1:4), each = 5)

country_id <- rep(seq_len(J), length.out = N)

# Covariates
victory <- rbinom(N, 1, 0.3)
gdp <- rnorm(N)
peacekeeping <- rbinom(N, 1, 0.4)

# True parameters
beta <- c(-2.0, -0.5, 0.3, -0.8)   # intercept, victory, gdp, peacekeeping
gamma <- c(0.8, -0.6)               # intercept, gdp (split stage)
shape <- 1.2

# Frailties
u <- rnorm(J, 0, 0.3)
u <- u - mean(u)

# Split stage
X <- cbind(1, victory, gdp, peacekeeping)
Z <- cbind(1, gdp)
pi <- plogis(drop(Z %*% gamma))
susceptible <- rbinom(N, 1, pi)

# Survival stage
eta <- drop(X %*% beta) + u[country_id]
q <- runif(N)
latent <- (-log1p(-q) / exp(eta))^(1 / shape)
latent[susceptible == 0] <- Inf

censor <- rexp(N, rate = 0.08)
time <- pmin(latent, censor)
event <- as.integer(latent <= censor)

conflict_spells <- data.frame(
  time = round(time, 2),
  event = event,
  victory = victory,
  gdp = round(gdp, 3),
  peacekeeping = peacekeeping,
  country = factor(countries[country_id]),
  region = factor(regions[country_id]),
  country_id = country_id
)

usethis::use_data(conflict_spells, overwrite = TRUE)
