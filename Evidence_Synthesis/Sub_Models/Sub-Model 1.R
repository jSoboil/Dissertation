# ====================================================================================
# All Cause Mortality Model -----------------------------------------------
# ====================================================================================
library(R2jags)
library(rjags)
library(readxl)
library(parallel)

options(mc.cores = detectCores())

# Mortality data for female population from ASSA model: 
mort_data <- read_excel("mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

N <- round(as.matrix(mort_data[, 1]), digits = 0)
Dead <- round(as.matrix(mort_data[, 2]), digits = 0)
Dead / N
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
 
}
"
writeLines(text = model_string, con = "mortProb.txt")

# Data format that is readable for JAGS:
data_JAGS <- list(mort = mort, N.pop = N.pop)
data_JAGS
# Initial starting values for JAGS sampler:
inits <- list(
 list(rho.pop = 1, delta = rep(.09, 91), 
      mu.pop = rep(.1, 91)),
 list(rho.pop = 0, delta = rep(1, 91), 
      mu.pop = rep(.5, 91))
)
# Parameters to monitor:
params <- c("d", "OR", "pEfficacy")

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = c("pop.pi"), 
                 model.file = "mortProb.txt", inits = inits, n.chains = 2, 
                 n.iter = 10000, n.burnin = 1000)
mod_JAGS
