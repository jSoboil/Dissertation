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

# This script adds sub-model 6 to the overall model. Specifically, all progression from Cancer
# states to cancer survivor states.

source("parameter_Inputs.R")

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
  for (i in 1:85) {
    # Monte Carlo:
    # No vaccine:
    omega.age[i] ~ dlnorm(mu.a.log[i], prec.age[i])
    # With vaccine:
    omega.age.vac[i] <- omega.age[i] * pEfficacy.vac
    
    # Note in use of pow() function, using -2 is a shorthand inverse
    # method equivalent to 1 / x^2.
    log(prec.age[i]) <- pow(sigma.age[i], -2)
    # Prior on variance for each age group. Note use of half Student-t to draw
    # variance away from 0. See Prior distribution for variance parameters in 
    # hierarchical models (Gelman, 2006):
    sigma.age[i] ~ dt(0, eta.age, 1)T(0, )
  }
  
   # Wide hyper-prior on prior variance parameter for SUB-MODEL 1:
   eta.age ~ dunif(0, 100)
 
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

# SUB-MODEL 3: CANCER PROGRESSION AND 5-YEAR SURVIVAL STAGES I-IV.
# Model parameters abbreviated by .canc. Note: this is equivalent to a standard 
# Monte Carlo PSA, as it is technically sampling directly from a prior and it is
# *not* propogated into a posterior using a likelihood model. 

   # Distribution according to Stage:
    Stage.I.canc ~ dbeta(alpha.StageI, beta.StageI)
    Stage.II.canc ~ dbeta(alpha.StageII, beta.StageII)
    Stage.III.canc ~ dbeta(alpha.StageIII, beta.StageIII)
    StageIV.Detected ~ dbeta(alpha.StageIV, beta.StageIV)
    
    # Functional parameters for Cancer stages:
    StageI.Detected <- Stage.I.canc * 0.15
    StageII.Detected <- Stage.II.canc * 0.225
    
   # Stage I Cervical Cancer detected 5-year Survival:
   surv.StageI_year1 ~ dbeta(alpha.StageI_YearI, beta.StageI_YearI)
   surv.StageI_year2 ~ dbeta(alpha.StageI_YearII, beta.StageI_YearII)
   surv.StageI_year3 ~ dbeta(alpha.StageI_YearIII, beta.StageI_YearIII)
   surv.StageI_year4 ~ dbeta(alpha.StageI_YearIV, beta.StageI_YearIV)
   surv.StageI_year5 ~ dbeta(alpha.StageI_YearV, beta.StageI_YearV)
   # Stage II Cervical Cancer detected 5-year Survival:
   surv.StageII_year1 ~ dbeta(alpha.StageII_YearI, beta.StageII_YearI)
   surv.StageII_year2 ~ dbeta(alpha.StageII_YearII, beta.StageII_YearII)
   surv.StageII_year3 ~ dbeta(alpha.StageII_YearIII, beta.StageII_YearIII)
   surv.StageII_year4 ~ dbeta(alpha.StageII_YearIV, beta.StageII_YearIV)
   surv.StageII_year5 ~ dbeta(alpha.StageII_YearV, beta.StageII_YearV)
   # Stage III Cervical Cancer detected 5-year Survival:
   surv.StageIII_year1 ~ dbeta(alpha.StageIII_YearI, beta.StageIII_YearI)
   surv.StageIII_year2 ~ dbeta(alpha.StageIII_YearII, beta.StageIII_YearII)
   surv.StageIII_year3 ~ dbeta(alpha.StageIII_YearIII, beta.StageIII_YearIII)
   surv.StageIII_year4 ~ dbeta(alpha.StageIII_YearIV, beta.StageIII_YearIV)
   surv.StageIII_year5 ~ dbeta(alpha.StageIII_YearV, beta.StageIII_YearV)
   # Stage IV Cervical Cancer detected 5-year Survival:
   surv.StageIV_year1 ~ dbeta(alpha.StageIV_YearI, beta.StageIV_YearI)
   surv.StageIV_year2 ~ dbeta(alpha.StageIV_YearII, beta.StageIV_YearII)
   surv.StageIV_year3 ~ dbeta(alpha.StageIV_YearIII, beta.StageIV_YearIII)
   surv.StageIV_year4 ~ dbeta(alpha.StageIV_YearIV, beta.StageIV_YearIV)
   surv.StageIV_year5 ~ dbeta(alpha.StageIV_YearV, beta.StageIV_YearV)
   
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

# END OF SUB-MODEL 5.

 }
"
writeLines(text = model_String, con = "SUBMOD6.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
  Nstud.vac = Nstud.vac, 
  rA.vac = rA.vac, nA.vac = nA.vac, 
  rB.vac = rB.vac, nB.vac = nB.vac,
  
  # Population prevalence:
  mu.a.log = mu.a.log,
  
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
  alpha.HSIL = alpha.HSIL, beta.HSIL = beta.HSIL,
  
  # Cancer Stage alpha's:
  alpha.StageI = alpha.StageI, beta.StageI = beta.StageI,
  alpha.StageII = alpha.StageII, beta.StageII = beta.StageII,
  alpha.StageIII = alpha.StageIII, beta.StageIII = beta.StageIII,
  alpha.StageIV = alpha.StageIV, beta.StageIV = beta.StageIV,
  
  
  # Stage I 5-year Survival:
  alpha.StageI_YearI = alpha.StageI_YearI, beta.StageI_YearI = beta.StageI_YearI,
  alpha.StageI_YearII = alpha.StageI_YearII, beta.StageI_YearII = beta.StageI_YearII,
  alpha.StageI_YearIII = alpha.StageI_YearIII, beta.StageI_YearIII = beta.StageI_YearIII,
  alpha.StageI_YearIV = alpha.StageI_YearIV, beta.StageI_YearIV = beta.StageI_YearIV,
  alpha.StageI_YearV = alpha.StageI_YearV, beta.StageI_YearV = beta.StageI_YearV,
  # Stage II 5-year Survival:
  alpha.StageII_YearI = alpha.StageII_YearI, beta.StageII_YearI = beta.StageII_YearI,
  alpha.StageII_YearII = alpha.StageII_YearII, beta.StageII_YearII = beta.StageII_YearII,
  alpha.StageII_YearIII = alpha.StageII_YearIII, beta.StageII_YearIII = beta.StageII_YearIII,
  alpha.StageII_YearIV = alpha.StageII_YearIV, beta.StageII_YearIV = beta.StageII_YearIV,
  alpha.StageII_YearV = alpha.StageII_YearV, beta.StageII_YearV = beta.StageII_YearV,
  # Stage III 5-year Survival:
  alpha.StageIII_YearI = alpha.StageIII_YearI, beta.StageIII_YearI = beta.StageIII_YearI, 
  alpha.StageIII_YearII = alpha.StageIII_YearII, beta.StageIII_YearII = beta.StageIII_YearII,
  alpha.StageIII_YearIII = alpha.StageIII_YearIII, beta.StageIII_YearIII = beta.StageIII_YearIII,
  alpha.StageIII_YearIV = alpha.StageIII_YearIV, beta.StageIII_YearIV = beta.StageIII_YearIV,
  alpha.StageIII_YearV = alpha.StageIII_YearV, beta.StageIII_YearV = beta.StageIII_YearV,
  # Stage IV 5-year Survival:
  alpha.StageIV_YearI = alpha.StageIV_YearI, beta.StageIV_YearI = beta.StageIV_YearI,
  alpha.StageIV_YearII = alpha.StageIV_YearII, beta.StageIV_YearII = beta.StageIV_YearII,
  alpha.StageIV_YearIII = alpha.StageIV_YearIII, beta.StageIV_YearIII = beta.StageIV_YearIII,
  alpha.StageIV_YearIV = alpha.StageIV_YearIV, beta.StageIV_YearIV = beta.StageIV_YearIV,
  alpha.StageIV_YearV = alpha.StageIV_YearV, beta.StageIV_YearV = beta.StageIV_YearV
  
)

# Parameters to monitor:
params <- c(
  # Vaccine efficacy parameters:
  "OR.vac", "pEfficacy.vac",
  # Well to infection prevalence:
  "omega.age", "omega.age.vac",
  # HPV Progression:
  "HPV_Well_15to24", "HPV_Well_25to29",
  "HPV_Well_30toEnd",
  # LSIL and HSIL Progression:
  "LSIL_15", "LSIL_35",
  "HSIL_n",
   # Cancer Stage Progression:
  "Stage.I.canc", "Stage.II.canc",
  "Stage.III.canc", "StageIV.Detected",
  "StageI.Detected", "StageII.Detected",
  # Cancer survival progression:
  "surv.StageI_year1", "surv.StageI_year2", "surv.StageI_year3", 
  "surv.StageI_year4", "surv.StageI_year5",
  "surv.StageII_year1", "surv.StageII_year2", "surv.StageII_year3", 
  "surv.StageII_year4", "surv.StageII_year5",
  "surv.StageIII_year1", "surv.StageIII_year2", "surv.StageIII_year3",
  "surv.StageIII_year4", "surv.StageIII_year5",
  "surv.StageIV_year1", "surv.StageIV_year2", "surv.StageIV_year3",
  "surv.StageIV_year4", "surv.StageIV_year5"
  )

# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 30000
n.burnin <- 10000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "SUBMOD6.txt", n.chains = 6, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS

attach.jags(mod_JAGS)

# Death from cancer, for example:
1 - (surv.StageI_year1 - v_p_mort_lessHPV[1, 2])

# Stage I-IV equations:
StageItoStageII <- (Stage.I.canc - (Stage.I.canc * 0.15))
StageItoStageII

StageItoTreat <- (Stage.I.canc * 0.15)
StageItoTreat

StageItoDeath <- v_p_mort_lessHPV[22, 2]

StageItoStageI <- 1 - (StageItoDeath + StageItoTreat + StageItoStageII)

# LOTP Check:
(StageItoDeath + StageItoTreat + StageItoStageII + StageItoStageI)

# ==========================================================================================
# Markov Model setup: -----------------------------------------------------------
# ==========================================================================================
n_age_init <- 0 # age at baseline
n_age_max <- 85 # maximum age of follow up
n_t <- n_age_max - n_age_init # time horizon, number of cycles

# the 30 health states of the model:
v_n <- c("Well", "Infection", "LSIL", "HSIL", "Stage-I Cancer", "Stage-II Cancer",
         "Stage-III Cancer", "Stage-IV Cancer", "Detected.Stage-I Year 1",
         "Detected.Stage-I Year 2", "Detected.Stage-I Year 3", "Detected.Stage-I Year 4",
         "Detected.Stage-I Year 5", "Detected.Stage-II Year 1", "Detected.Stage-II Year 2",
         "Detected.Stage-II Year 3", "Detected.Stage-II Year 4", "Detected.Stage-II Year 5",
         "Detected.Stage-III Year 1", "Detected.Stage-III Year 2", "Detected.Stage-III Year 3",
         "Detected.Stage-III Year 4", "Detected.Stage-III Year 5", "Detected.Stage-IV Year 1",
         "Detected.Stage-IV Year 2", "Detected.Stage-IV Year 3", "Detected.Stage-IV Year 4",
         "Detected.Stage-IV Year 5", "Cancer Survivor", "Death")
n_states <- length(v_n) # number of health states 

d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs
v_names_str <- c("Status quo", "Bivalent HPV Vaccine") # strategy names

a_P <- array(0, dim = c(n_states, n_states, n_t, n.sims),
             dimnames = list(v_n, v_n, 0:(n_t - 1), 1:n.sims))
str(a_P)

# ==========================================================================================
# Fill Transition Array ------------------------------------------------
# ==========================================================================================
# Status Quo  ----------------------------------------------------------


# New Treatment  ----------------------------------------------------------
str(a_P)
str(omega.age)
for (i in 1:85) {
 for (j in 1:n.sims) {
  a_P["Well", "Infection", i, ] <- omega.age[j, i] * (1 - pEfficacy.vac[j, ])
 }
}
a_P[, , 19, ]

str(pEfficacy.vac)

for (j in 1:n.sims) {
 treatment <- omega.age[j, ] * (1- pEfficacy.vac[j, ])
}
treatment



a_P["Stable", "Progression", , ] <- pi_tox[, 2]
a_P["Stable", "Response", , ] <- (1 - pi_tox[, 2]) * beta_tr[, 2]
a_P["Stable", "Death", , ] <- 1 - ((1 - pi_tox[, 2]) * (1 - beta_tr[, 2]) + pi_tox[, 2] + 
                                    (1 - pi_tox[, 2]) * beta_tr[, 2])
# All transitions from the state Response
a_P["Response", "Response", , ] <- ((1 - pi_tox[, 2] ) * pi_res[, 2])
a_P["Response", "Progression", , ] <-  (pi_tox[, 2]) + (1 - pi_tox[, 2]) * (1 - pi_res[, 2])

# All transitions from the state Progression
a_P["Progression", "Progression", , ] <- 1 - beta_dth[, 2]
a_P["Progression", "Death", , ] <- beta_dth[, 2]

# All transitions from the state Death
a_P["Death", "Death", , ] <- 1
a_P





