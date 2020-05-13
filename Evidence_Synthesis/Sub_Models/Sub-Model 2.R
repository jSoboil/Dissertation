library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(parallel)

options(mc.cores = detectCores())

# ==========================================================================================
# Age-specific HPV incidence ----------------------------------------------------
# ==========================================================================================

# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
group_Age <- c(rep(age_group[1], length(15:16)), rep(age_group[2], 1), 
           rep(age_group[3], 1), rep(age_group[4], 1), rep(age_group[5], 1), 
           rep(age_group[6], 1), rep(age_group[7], length(22:23)), 
           rep(age_group[8], length(24:29)), rep(age_group[9], length(30:49)), 
           rep(age_group[10], length(50:104)))
group_Age

incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))
r_Age <- c(rep(incidence[1], length(15:16)), rep(incidence[2], 1), 
           rep(incidence[3], 1), rep(incidence[4], 1), rep(incidence[5], 1), 
           rep(incidence[6], 1), rep(incidence[7], length(22:23)), 
           rep(incidence[8], length(24:29)), rep(incidence[9], length(30:49)), 
           rep(incidence[10], length(50:104)))
r_Age
barplot(r_Age, names.arg = "From Age 15 to 100", 
        main = "Rate of HPV+", sub =  "At constant t-cycle of 1-year")

## Copyright Gianluca Baio 2012
lognPar <- function(mu, sdev) {
	sdev.sq <- sdev ^ 2
	mu.log <- log(mu) - .5 * log(1 + sdev.sq / mu ^ 2)
	sdev.sq.log <- log(1 + (sdev.sq / mu ^ 2))
	sigma.log <- sqrt(sdev.sq.log)
	list(mu.log = mu.log, sigma.log= sigma.log)
}

r_mu.log <- lognPar(mu = r_Age, sdev = r_Age)$mu.log
r_mu.log
r_sigma.log <- 1 / lognPar(mu = r_Age, sdev = r_Age)$sigma.log ^ 2
r_sigma.log

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
 sigma.pop ~ dunif(0, 20)

# Sub-model 2: age-specific probability of infection:
  # model parameters abbreviated by .age
  
  # Probability of infection for 1 cycle:
    # note: equivalent of sampling directly from prior.
  for (i in 1:90) {
    omega.age[i] ~ dlnorm(r_mu.log[i], r_sigma.log[i])
 }
 
}
"
writeLines(text = model_string, con = "mortProb.txt")

data_JAGS <- list(mort = mort, N.pop = N.pop, r_mu.log = r_mu.log, 
                  r_sigma.log = r_sigma.log)
data_JAGS

jags(data = data_JAGS, parameters.to.save = c("omega.age"), 
     model.file = "mortProb.txt", n.chains = 2, n.iter = 10000, n.burnin = 1000)
