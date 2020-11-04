library(tidyverse)
library(dampack)
library(readxl)
library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(MCMCpack)
library(Compositional)

# This script adds sub-model 5 to the overall model. Specifically, states LSIL to HPV,
# LSIL to HSIL, LSIL to LSIL, HSIL to LSIL, HSIL to Normal/Well, HSIL to Stage-I Cancer, 
# and HSIL to HSIL. Like the states from Infection, this can be modelled from a single
# distribution as the rest of the probabilities are conditional on moving out the state.

# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# ====================================================================================
# Sub-Model for progressions from LSIL --------------------
# ====================================================================================
# To Infection or Normal ----------------------------------------
# Ages 15-34:
LSIL_15to34 <- 1 - exp(-.65 * 6)
alpha.LSIL_15to34 <- beta_params(mean = LSIL_15to34, sigma = 0.05)$alpha
beta.LSIL_15to34 <- beta_params(mean = LSIL_15to34, sigma = 0.05)$beta

# Proportion reverting to Normal is .9:
LSILtoNormal_15to34 <- ((1 - exp(-.65 * 6)) * .9)
# Proportion reverting to HPV/Infection
LSILtoHPV_15to34 <- (1 - exp(-.65 * 6)) - ((1 - exp(-.65 * 6)) * .9)

# Ages ≥ 35:
LSIL_35up <- 1 - exp(-.4 * 6)
alpha.LSIL_35up <- beta_params(mean = LSIL_35up, sigma = 0.05)$alpha
beta.LSIL_35up <- beta_params(mean = LSIL_35up, sigma = 0.05)$beta
# Proportion reverting to Normal is .9:
LSILtoNormal_35up <- ((1 - exp(-.4 * 6)) * .9)
# Proportion reverting to HPV/Infection
LSILtoHPV_35up <- (1 - exp(-.4 * 6)) - ((1 - exp(-.4 * 6)) * .9)

# To HSIL -------------------------------------------------------
# Conditional on not regressing to Normal or Infection...
# Ages 15-34:
LSILtoHSIL_15to34 <- (1 - (1 - exp(-.65 * 6))) * (1 - exp(-.1 * 6))

# Ages ≥ 35:
LSILtoHSIL_35up <- (1 - (1 - exp(-.4 * 6))) * (1 - exp(-.35 * 6))

# To LSIL -----------------------------------------------------------------
# Setting 1 - \sum{p_i}
LSILtoLSIL_15to34 <- 1 - (LSILtoHPV_15to34 + LSILtoHSIL_15to34 + LSILtoNormal_15to34)

# Calibration for LSIL State ----------------------------------------------
LSILtoLSIL_15to34 + LSILtoHPV_15to34 + LSILtoHSIL_15to34 + LSILtoNormal_15to34

# ====================================================================================
# Sub-Model for progressions from HSIL --------------------
# ====================================================================================
HSIL <- (1 - exp(-.35 * 6))
alpha.HSIL <- beta_params(mean = HSIL, sigma = 0.05)$alpha
beta.HSIL <- beta_params(mean = HSIL, sigma = 0.05)$beta

# To Normal ----------------------------------------
HSILtoNormal <- ((1 - exp(-.35 * 6)) * .5)

# To LSIL ----------------------------------------
HSILtoLSIL <- (1 - exp(-.35 * 6)) - ((1 - exp(-.35 * 6)) * .5)

# To Stage-I Cancer -------------------------------------------------------
# Probability of progression every 10 years = .4
# Convert to yearly rate
- (1 / 10) * log(1 - .4)
# Convert to annual progression probability to Stage I Cancer:
HSILtoStageI <- 1 - exp(-0.05108256 * 1)

# To HSIL -----------------------------------------------------------------
HSILtoHSIL <- 1 - (HSILtoNormal + HSILtoLSIL + HSILtoStageI)
HSILtoHSIL

# Calibration for HSIL State ----------------------------------------------
HSILtoHSIL + HSILtoLSIL + HSILtoNormal + HSILtoStageI


# This script adds sub-model 4 to the overall model. Specifically, states Normal/Well to HPV,
# HPV to LSIL, HPV to HSIL.

# ==========================================================================================
# Normal/Well State Progression -------------------------------------------
# ==========================================================================================

# Normal to HPV -----------------------------------------------------------
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

# Optional 'emprical bayes' method:
# mu.a.log <- lnorm_params(m = Prevalence, v = .01)$mu
# sigma.a.log <- lnorm_params(m = Prevalence, v = .01)$sigma
# prec.age <- 1 / (sigma.a.log * sigma.a.log)
# prec.age

# Use rep(omega.age[, i], N.age_group[i])

# Normal to Death ---------------------------------------------------------
# Import ASSA mortality table:
mort_data <- read_excel(
"/Users/joshuamusson/Desktop/Analytics/R/Dissertation/Evidence_Synthesis/mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

# Total Pop. divided by total deaths:
v_r_mort_by_age <- mort_data[, 2] / mort_data[, 1]
v_r_mort_by_age[1:90, ]
plot(unlist(v_r_mort_by_age), type = "l", lwd = 4)

# Ages 1:14
annual.mortality_1to14 <- v_r_mort_by_age[1:14, ]
# Ages 15:24
annual.mortality_15to24 <- v_r_mort_by_age[15:24, ]
# Ages 25:29
annual.mortality_25to29 <- v_r_mort_by_age[25:29, ]
# Ages ≥30:
annual.mortality_30plus <- v_r_mort_by_age[30:91, ]

# Normal to Normal --------------------------------------------------------
# 1 - (Age mortality[, i] + omega.age[, i])

# ==========================================================================================
# HPV/Infection State Progression -----------------------------------------
# ==========================================================================================
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Below are the equations for relationships between all states from HPV/infection:

# HPV to Normal -----------------------------------------------------------
# Only way I can make sense of this is if it is conditional on not regressing to Normal.
# Conditional grouping abbreviated by starting age group.

# Ages 15-24:
HPV_Normal_15 <- (1 - exp(-0.7 * 1.5))
1 - HPV_Normal_15

# Ages 25-29:
HPV_Normal_25 <- (1 - exp(-.5 * 1.5))
1 - HPV_Normal_25

# Ages ≥30:
HPV_Normal_30plus <- (1 - exp(-.15 * 1.5))
1 - HPV_Normal_30up

# Parameters for this regression:
alpha.HPVtoNormal_15to24 <- beta_params(mean = HPV_Normal_15, sigma = 0.1)$alpha
beta.HPVtoNormal_15to24 <- beta_params(mean = HPV_Normal_15, sigma = 0.1)$beta

alpha.HPVtoNormal_25to29 <- beta_params(mean = HPV_Normal_25, sigma = 0.1)$alpha
beta.HPVtoNormal_25to29 <- beta_params(mean = HPV_Normal_25, sigma = 0.1)$beta

alpha.HPVtoNormal_30toPlus <- beta_params(mean = HPV_Normal_30plus, sigma = 0.1)$alpha
beta.HPVtoNormal_30toPlus <- beta_params(mean = HPV_Normal_30plus, sigma = 0.1)$beta

# HPV to LSIL -------------------------------------------------------------
# Ages 15-24:
HPV_LSIL_15 <- (1 - HPV_Normal_15) * ((1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))
HPV_LSIL_15

# Ages 25-29:
HPV_LSIL_25 <- (1 - HPV_Normal_25) * ((1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))
HPV_LSIL_25

# Ages ≥30:
HPV_LSIL_30plus <- (1 - HPV_Normal_30plus) * ((1 - exp(-.2 * 3)) - 
                                               ((1 - exp(-.2 * 3)) * .1))
HPV_LSIL_30plus


# HPV to HSIL -------------------------------------------------------------
# Ages 15-24:
HPV_HSIL_15 <- ((1 - HPV_Normal_15) * ((1 - exp(-.2 * 3)) * .1))
HPV_HSIL_15

# Ages 25-29:
HPV_HSIL_25 <- ((1 - HPV_Normal_25) * ((1 - exp(-.2 * 3)) * .1))
HPV_HSIL_25

# Ages ≥30:
HPV_HSIL_30plus <- ((1 - HPV_Normal_30plus) * ((1 - exp(-.2 * 3)) * .1))
HPV_HSIL_30plus

# HPV to Death ------------------------------------------------------------
# Death not dependent on regression to normal...

# Ages 15-24:
HPV_Death_15 <- v_r_mort_by_age[15:24, ]
HPV_Death_15

# Ages 25-29:
HPV_Death_25 <- v_r_mort_by_age[25:29, ]
HPV_Death_25

# Ages ≥30:
HPV_Death_30plus <- v_r_mort_by_age[30:91, ]
HPV_Death_30plus

# HPV to HPV --------------------------------------------------------------

# Ages 15-24:
HPV_HPV_15 <- (1 - HPV_Normal_15) - (HPV_Death_15 + HPV_LSIL_15 + 
                                HPV_HSIL_15)

# Ages 25-29:
HPV_HPV_25 <- (1 - HPV_Normal_25) - (HPV_Death_25 + HPV_LSIL_25 + 
                                HPV_HSIL_25)

# Ages ≥30:
HPV_HPV_30up <- (1 - HPV_Normal_30plus) - (HPV_Death_30plus + HPV_LSIL_30plus + 
                                HPV_HSIL_30plus)

# State Calibration -------------------------------------------------------------

# Calibration - Ages 15-34. Must equal 1 across all ages:
(HPV_HPV_15 + HPV_HSIL_15 + HPV_LSIL_15 + HPV_Death_15) + HPV_Normal_15

# Calibration - Ages 25-29. Must equal 1 across all ages:
HPV_HPV_25 + HPV_HSIL_25 + HPV_LSIL_25 + HPV_Death_25 + HPV_Normal_25

# Calibration - Ages ≥30. Must equal 1 across all ages:
HPV_HPV_30up + HPV_HSIL_30plus + HPV_LSIL_30plus + HPV_Death_30plus + HPV_Normal_30plus

# ==========================================================================================
# Cancer State Progression ------------------------------------------------
# ==========================================================================================

# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Distribution of Stage I progression:
# alpha = (0.4376587 + 0.15 + 0.4123413)
# Stage_I ~ ddirch(alpha)
alpha.Stage_I <- c(0.4376587, 0.15, 0.4123413)

# Distribution of Stage II progression:
# alpha = (0.5358411 + 0.225 + 0.2391589)
# Stage_II ~ ddirch(alpha)
alpha.Stage_II <- c(0.5358411, 0.225, 0.2391589)

# Distribution of Stage III progression:
# alpha = (0.6837724 + 0.1897366 + 0.126491)
# Stage_III ~ ddirch(alpha)
alpha.Stage_III <- c(0.6837724, 0.1897366, 0.126491)

# Distribution of Stage IV progression:
# alpha = (0.9 + 0.1)
# Stage_IV ~ ddirch(alpha)
alpha.Stage_IV <- beta_params(mean = .9, sigma = .05)$alpha
beta.Stage_IV <- beta_params(mean = .9, sigma = .05)$beta

# ==========================================================================================
# Vaccine efficacy data ----------------------------------------------------
# ==========================================================================================
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
# Sub-model 3 ------------------------------------------------------
# ==========================================================================================
model_String <- "
model {

# SUB-MODEL 1: AGE-SPECIFIC PREVALENCE
# Model parameters abbreviated by .age. Note: this is equivalent to Monte Carlo PSA, 
# as it is technically sampling directly from a prior and is not propogated into a 
# posterior using a likelihood model. However, a hyperprior is used for the population 
# variance to account for a greater uncertainty in each age-population.
  for (i in 1:10) {
    # Monte Carlo:
    omega.age[i] ~ dlnorm(mu.a.log[i], prec.age[i])
    
    # Note in use of pow() function, using -2 is a shorthand inverse
    # method equivalent to 1 / x^2.
    log(prec.age[i]) <- pow(sigma.age[i], -2)
    # Prior on variance for each age group. Note use of half Student-t to draw
    # variance away from 0. See Prior distribution for variance parameters in 
    # hierarchical models (Gelman, 2006):
    sigma.age[i] ~ dt(0, eta.age, 1)T(0, )
  }
  
   # Wide hyper-prior on prior variance parameter for SUB-MODEL 1:
   eta.age ~ dunif(0, 1000)
 
# END OF SUB-MODEL 1.

# SUB-MODEL 2: VACCINE-EFFICACY.
# Model parameters abbreviated by .vac.
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
  
   # Hyperpriors for SUB-MODEL 2:
   psi.vac ~ dnorm(0, 1.0e-6)
   prec.vac <- pow(tau.vac, -2)
   tau.vac ~  dunif(0, 10)
  
  # Transformations for SUB-MODEL 2:
   # Convert LOR to OR
   OR.vac <- exp(psi.vac)
   # Convert OR to probability
   # for vaccine efficacy
   pEfficacy.vac <- 1 / (1 + OR.vac)

# END OF SUB-MODEL 2.

# SUB-MODEL 3: CANCER PROGRESSION STAGES I-IV.
# Model parameters abbreviated by .canc. Note: this is equivalent to a standard 
# Monte Carlo PSA, as it is technically sampling directly from a prior and it is
# *not* propogated into a posterior using a likelihood model. 

   # Monte Carlo:
    Stage.I.canc ~ ddirch(alpha.Stage_I)
    Stage.II.canc ~ ddirch(alpha.Stage_II)
    Stage.III.canc ~ ddirch(alpha.Stage_III)
    Stage.IV.canc ~ dbeta(alpha.Stage_IV, beta.Stage_IV)

# END OF SUB-MODEL 3.

# SUB-MODEL 4: INFECTION PROGRESSION:
# Note: this is equivalent to a standard Monte Carlo PSA, as it is technically sampling
# directly from a prior and it is *not* propogated into a posterior using a likelihood 
# model. 
   
   # Monte Carlo:
   # From HPV to Normal across age groups:
   # Note: because all other states except Death are assumed to be dependent and disjoint for
   # regression to normal from state of HPV/Infection, one can calculate all other relevant 
   # states from the complement of the transitions that are obtained from the model below:
    HPV_Well_15to24 ~ dbeta(alpha.HPVtoNormal_15to24, beta.HPVtoNormal_15to24)
    HPV_Well_25to29 ~ dbeta(alpha.HPVtoNormal_25to29, beta.HPVtoNormal_25to29)
    HPV_Well_30toEnd ~ dbeta(alpha.HPVtoNormal_30toPlus, beta.HPVtoNormal_30toPlus)

# END OF SUB-MODEL 4.

# SUB-MODEL 5: LSIL & HSIL PROGRESSION:
# Note: this is equivalent to a standard Monte Carlo PSA, as it is technically sampling
# directly from a prior and it is *not* propogated into a posterior using a likelihood 
# model. 
   LSIL_15 ~ dbeta(alpha.LSIL_15to34, beta.LSIL_15to34)
   LSIL_35 ~ dbeta(alpha.LSIL_35up, beta.LSIL_35up)
   HSIL_n ~ dbeta(alpha.HSIL, beta.HSIL)


 }
"
writeLines(text = model_String, con = "SUBMOD5.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
  Nstud.vac = Nstud.vac, 
  rA.vac = rA.vac, nA.vac = nA.vac, 
  rB.vac = rB.vac, nB.vac = nB.vac,
  # Population prevalence:
  mu.a.log = mu.a.log,
  
  # Cancer Stage alpha's:
  alpha.Stage_I = alpha.Stage_I, alpha.Stage_II = alpha.Stage_II,
  alpha.Stage_III = alpha.Stage_III, 
  alpha.Stage_IV = alpha.Stage_IV, beta.Stage_IV = beta.Stage_IV,
  
  # HPV to Normal data:
  alpha.HPVtoNormal_15to24 = alpha.HPVtoNormal_15to24, 
  beta.HPVtoNormal_15to24 = beta.HPVtoNormal_15to24,
  alpha.HPVtoNormal_25to29 = alpha.HPVtoNormal_25to29, 
  beta.HPVtoNormal_25to29 = beta.HPVtoNormal_25to29,
  alpha.HPVtoNormal_30toPlus = alpha.HPVtoNormal_30toPlus,
  beta.HPVtoNormal_30toPlus = beta.HPVtoNormal_30toPlus,
  
  # LSIL and HSIL transition data:
  alpha.LSIL_15to34 = alpha.LSIL_15to34, beta.LSIL_15to34 = beta.LSIL_15to34,
  alpha.LSIL_35up = alpha.LSIL_35up, beta.LSIL_35up = beta.LSIL_35up,
  alpha.HSIL = alpha.HSIL, beta.HSIL = beta.HSIL
  
)

# Parameters to monitor:
params <- c(
  # Vaccine efficacy parameters:
  "OR.vac", "pEfficacy.vac", 
  # Age prevalence:
  "omega.age",
  # Cancer Stage Progression:
  "Stage.I.canc", "Stage.II.canc",
  "Stage.III.canc", "Stage.IV.canc",
  # HPV to Normal Progression:
  "HPV_Well_15to24", "HPV_Well_25to29",
  "HPV_Well_30toEnd",
  # LSIL and HSIL Progression:
  "LSIL_15", "LSIL_35",
  "HSIL_n"
  
  )

# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 25000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "SUBMOD5.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS

attach.jags(mod_JAGS)

# Equations for state transitions:
# Proportion reverting to Normal is .9:
LSILtoNormal_15 <- LSIL_15 * .9
LSILtoNormal_15
# Proportion reverting to HPV/Infection
LSILtoHPV_15 <- LSIL_15 - LSILtoNormal_15
LSILtoHPV_15
# To HSIL:
LSILtoHSIL_15 <- (1 - LSIL_15) * (1 - exp(-.1 * 6))
LSILtoHSIL_15
# To LSIL:
LSILtoLSIL_15 <- 1 - (LSILtoHPV_15 + LSILtoHSIL_15 + LSILtoNormal_15)

# Calibration for SUB-MODEL 5 ----------------------------------------------
LSILtoLSIL_15 + LSILtoHPV_15 + LSILtoHSIL_15 + LSILtoNormal_15
# Yay!