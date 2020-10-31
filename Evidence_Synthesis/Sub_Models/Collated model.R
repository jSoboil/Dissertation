library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)

# Vaccine efficacy data ----------------------------------------------------
# PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 11
rA.vac <- c(435, 41, 28, 38, 219, 15, 10, 21, 56, 20, 17)
nA.vac <- c(5375, 355, 277, 2239, 2924, 392, 175, 7838, 7312, 372, 2502)
rB.vac <- c(32, 12, 1, 3, 61, 0, 0, 2, 4, 0, 1)
nB.vac <- c(5406, 366, 310, 2190, 2910, 387, 193, 7788, 7344, 401, 2497)
# Check:
length(rA.vac) == length(nA.vac)
length(rB.vac) == length(nB.vac)

# ==========================================================================================
# Age-specific infection ----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# Informative prior sampling model:
# omega.age[i] ~ dlnorm(mu.log[i], prec.log[i])

# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤21", "22-23", 
               "24-29", "30-49", "≥55")

# Estimated Prevalence:
Prevalence <- c(.09516258, .1130796, .139292, .1563352, .139292, 
                .1130796, .09516258, .04877058, .009950166, .004987521)

mu.a.log <- log(Prevalence)

# 'Emprical bayes' method:
# mu.a.log <- lnorm_params(m = Prevalence, v = .01)$mu
# sigma.a.log <- lnorm_params(m = Prevalence, v = .01)$sigma
# prec.age <- 1 / (sigma.a.log * sigma.a.log)
# prec.age

# ==========================================================================================
# Sub-model 1 ------------------------------------------------------
# ==========================================================================================
model_String <- "
model {

# SUB-MODEL 1: AGE-SPECIFIC PREVALENCE - model parameters abbreviated by .age. Note: 
# this is equivalent to Monte Carlo PSA, as it is technically sampling directly from 
# a prior and is not propogated into a posterior using a likelihood model. However,
# a hyperprior is used for the population variance to account for a greater uncertainty.
  for (i in 1:10) {
  # Monte Carlo model:
    omega.age[i] ~ dlnorm(mu.a.log[i], prec.age[i])
    
    # Note use of pow() function, -2 is inverse method,
    # i.e. 1 / x^2
    log(prec.age[i]) <- pow(sigma.age[i], -2)
    sigma.age[i] ~ dt(0, eta.age, 1)T(0, )
  }
  
  # Hyper-priors on pop. variance:
  eta.age ~ dunif(0, 100)
 
# END OF SUB-MODEL 1.

# SUB-MODEL 2: VACCINE-EFFICACY - model parameters abbreviated by .vac.
  for (i in 1:Nstud.vac) {
    # Likelihood:
    rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
    rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])
    
    # Random Effect Logistic model:
    logit(pA.vac[i]) <- mu.vac[i]
    logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
    
    # Average effect prior for sub-model 2:
    mu.vac[i] ~ dnorm(0, 1e-6)
    # Prior for sub-model 2 (Random. pop. effect):
    delta.vac[i] ~ dnorm(psi.vac, prec.vac)
  }
  
   # Hyperpriors for Sub-model 2:
   psi.vac ~ dnorm(0, 1.0e-6)
   prec.vac <- pow(tau.vac, -2)
   tau.vac ~  dunif(0, 10)
  
  # Transformations for Sub-model 2
  
   # Convert LOR to OR
   OR.vac <- exp(psi.vac)
   # Convert OR to probability
   # for vaccine efficacy
   pEfficacy.vac <- 1 / (1 + OR.vac)

# END OF SUB-MODEL 2.

# SUB-MODEL 3: 

 }
"
writeLines(text = model_String, con = "Age_and_Efficacy.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
  Nstud.vac = Nstud.vac, 
  rA.vac = rA.vac, nA.vac = nA.vac, 
  rB.vac = rB.vac, nB.vac = nB.vac,
  # Population prevalence:
  mu.a.log = mu.a.log
)

# Parameters to monitor:
params <- c(
  # Vaccine efficacy parameters:
  "OR.vac", "pEfficacy.vac", 
  # Age prevalence:
  "omega.age"
  )

# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 25000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "Age_and_Efficacy.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS