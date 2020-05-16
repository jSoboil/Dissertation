# ====================================================================================
# All Cause Mortality Model -----------------------------------------------
# ====================================================================================
library(R2jags)
library(rjags)
library(tidyverse)
library(bayesplot)
library(dampack)
library(readxl)
library(parallel)

options(mc.cores = detectCores())

# Mortality data for female population from ASSA model: 
mortFem_data <- read_excel("mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

N <- round(as.matrix(mortFem_data[, 1]), digits = 0)
Dead <- round(as.matrix(mortFem_data[, 2]), digits = 0)
N.pop <- as.numeric(unlist(N))
mort <- as.numeric(unlist(Dead))

model_string <-"
model {
# Sub-model 1: population probability for age-specific mortality
    # model parameters abbreviated by .pop

# Binomial Likelihood:
 for (i in 1:91) {
  mort[i] ~ dbin(pop.pi[i], N.pop[i])
  
  # Prior Sampling model:
  logit(pop.pi[i]) <-  mu.pop[i] + delta[i]
  
  # Priors on mort_pr:
  mu.pop[i] ~ dnorm(0, 1.0e-5)
  delta[i] ~ dnorm(rho.pop, 1 / sigma.pop ^ 2)
  
 }
 # Hyperpriors on delta:
 rho.pop ~ dnorm(0, 1.0e-5)
 sigma.pop ~ dunif(0, 10)
 
 # Odds ratio:
 OR <- exp(rho.pop)
 
}
"
writeLines(text = model_string, con = "mort_Prob.txt")

# Data format that is readable for JAGS:
data_JAGS <- list(mort = mort, N.pop = N.pop)
data_JAGS
# Initial starting values for JAGS sampler:
inits <- list(
 list(rho.pop = 0, delta = rep(.09, 91), 
      mu.pop = rep(.1, 91)),
 list(rho.pop = 0, delta = rep(1, 91), 
      mu.pop = rep(.5, 91))
)
# Parameters to monitor:
params <- c("OR", "pop.pi")

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "mort_Prob.txt", inits = inits, n.chains = 2, 
                 n.iter = 10000, n.burnin = 1000, n.thin = 18)
mod_JAGS
# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR", "pop.pi[63]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , 56:62], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("OR", "pop.pi[63]"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("OR", "pop.pi[63]"),
         lags = 250)

# End file ----------------------------------------------------------------