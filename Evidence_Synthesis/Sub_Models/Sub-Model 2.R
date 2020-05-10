setwd("/Users/joshuamusson/Desktop/Analytics/R/Intergrated-CEA-thesis/Evidence_Synthesis/Sub_Models/")
# ==========================================================================================
# Age-specific incidence --------------------------------------------------
# ==========================================================================================
library(rjags)
library(R2jags)
library(parallel)

options(mc.cores = detectCores())

# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))

# Convert rate to probability:
p_Age <- 1 - exp(-incidence * 1)

Pr_Age <- cbind(age_group, round(p_Age, 4))
Pr_Age

p_Age <- c(rep(p_Age[1], length(15:16)), rep(p_Age[2], 1), rep(p_Age[3], 1), 
           rep(p_Age[4], 1), rep(p_Age[5], 1), rep(p_Age[6], 1), 
           rep(p_Age[7], length(22:23)), rep(p_Age[8], length(24:29)), 
           rep(p_Age[9], length(30:49)), rep(p_Age[10], length(50:100)))
length(p_Age)
# I THINK IT IS BEST TO RUN THE VIRTUAL COHORT FROM AGE OF VACCINATION TO 100.
barplot(p_Age, names.arg = "From Age 15 to 100", ylab = "Probability of HPV+")

# Model -------------------------------------------------------------------
# Must make this a lognormal, so it is not truncated...
model_string <- "model {
 for (i in 1:86) {
 # Likelihood
 p_Age[i] ~ dnorm(mu[i], prec[i])T(0, )
 
 # Prior on mu:
 mu[i] ~ dnorm(rho, psi)
 # Precision function:
 prec[i] <- 1 / (beta * beta)
 }
# Priors
rho ~ dunif(0, 1)
phi ~ dunif(0, 10)
psi <- 1  / (phi * phi)
beta ~ dunif(0, 10)
}"
writeLines(text = model_string, con = "ASInc.txt")

data_JAGS <- list(p_Age = p_Age)
data_JAGS
jags(data = data_JAGS, parameters.to.save = "p_Age", model.file = "ASInc.txt")








