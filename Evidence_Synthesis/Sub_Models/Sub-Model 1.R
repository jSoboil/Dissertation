# setwd("/Users/joshuamusson/Desktop/Analytics/R/Intergrated-CEA-thesis/Evidence_Synthesis/Sub_Models/")

# ====================================================================================
# All Cause Mortality Model -----------------------------------------------
# ====================================================================================
library(R2jags)
library(rjags)
library(readxl)
library(parallel)

options(mc.cores = detectCores())

# Mortality data for female population from ASSA model: 
mort_data <- read_excel("Evidence_Synthesis/mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

N <- round(as.matrix(mort_data[, 1]), digits = 0)
Dead <- round(as.matrix(mort_data[, 2]), digits = 0)
Dead / N
N <- as.numeric(unlist(N))
Dead <- as.numeric(unlist(Dead))
typeof(N) == typeof(Dead)

# p(d) ~ dbin(dead, N)

model_string <-"
model {
# Binomial Likelihood:
 for (i in 1:91) {
  Dead[i] ~ dbin(p[i], N[i])
  
  # Prior Sampling model:
  p[i] ~ dbeta(alpha, beta)
  
 }
 # Priors on mortality:
 alpha ~ dunif(0, 10)
 beta ~ dunif(0, 100)
}
"
writeLines(text = model_string, con = "mortProb.txt")

data_JAGS <- list(Dead = Dead, N = N)
data_JAGS

jags(data = data_JAGS, parameters.to.save = "p", model.file = "mortProb.txt", n.chains = 2,
    n.iter = 10000, n.burnin = 1000)
