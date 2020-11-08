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
library(rbenchmark)

start_time <- Sys.time()

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
  for (i in 1:86) {
    # Monte Carlo:
    # No vaccine:
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
   eta.age ~ dunif(0, 10000)
 
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
   LSIL_15_34 ~ dbeta(alpha.LSIL_15to34, beta.LSIL_15to34)T(0, )
   LSIL_35_End ~ dbeta(alpha.LSIL_35up, beta.LSIL_35up)T(0, )
   HSIL_n ~ dbeta(alpha.HSIL, beta.HSIL)T(0, )

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
  "omega.age",
  # HPV Progression:
  "HPV_Well_15to24", "HPV_Well_25to29",
  "HPV_Well_30toEnd",
  # LSIL and HSIL Progression:
  "LSIL_15_34", "LSIL_35_End",
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
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "SUBMOD6.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS

attach.jags(mod_JAGS)

# ==========================================================================================
# Markov Model setup: -----------------------------------------------------------
# ==========================================================================================
n_age_init <- 0 # age at baseline
n_age_max <- 86 # maximum age of follow up - note: the maximum age is actually 85, however, 
# due to starting at 0, it has to be compensated for by adding n + 1.
n_t <- n_age_max - n_age_init # time horizon, number of cycles
n_t
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

# Transition array for each alternative:
a_P_1 <- array(0, dim = c(n_states, n_states, n_t, n.sims),
             dimnames = list(v_n, v_n, 0:(n_t - 1), 1:n.sims))
a_P_2 <- array(0, dim = c(n_states, n_states, n_t, n.sims),
             dimnames = list(v_n, v_n, 0:(n_t - 1), 1:n.sims))
str(a_P_2)

# ==========================================================================================
# Transition Array ------------------------------------------------
# ==========================================================================================
# The following sections fill in the transition array probabilities across both the time
# horizon and all probabilistic simulations.

# First is the Status-quo treatment, which is non-vaccine probabilities. Second is the the 
# New treatment, which is vaccine associated probabilities. 

# ==========================================================================================
# Status-Quo  ----------------------------------------------------------
# ==========================================================================================
# These are all the state transitions with no vaccine treatment. The following creates a
# transition array that is then used to simulate the Markov cohort without treatment.

# Transitions from Well State ---------------------------------------------

# The following enters all transition probabilities for ages 0-85 for each transition from
# the state Well, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death-less-HPV
# within each individual health state.
for (i in 1:86) {
 for (j in 1:n.sims) {
  a_P_1["Well", "Death", i, j] <-  v_p_mort_lessHPV[i]
  
  a_P_1["Well", "Infection", i, j] <- (omega.age[j, i])
  
  a_P_1["Well", "Well", i, j] <- 1 - (v_p_mort_lessHPV[i] + omega.age[j, i])
  
 }
}

# Note: all transition probabilities for ages 0-14 are already filled by the natural structure
# of the transition array, i.e., all initial transitions are = 0 before being filled in by the
# proceeding code. This is why one can omit these transitions where applicable from the for 
# loop. Also see that this is why the mortality data for deaths from well state are the only
# data filled in for age groups 0-14.

# Transitions from Infection State ----------------------------------------

# The following enters all transition probabilities for ages 15-24 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
 for (i in 16:25) {
     for (j in 1:n.sims) {
      a_P_1["Infection", "Well", i, j] <- HPV_Well_15to24[j]
      
      a_P_1["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
      a_P_1["Infection", "LSIL", i, j] <- ((1 - HPV_Well_15to24[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1)))
      
      a_P_1["Infection", "HSIL", i, j] <- ((1 - HPV_Well_15to24[j]) * ((1 - exp(-.2 * 3)) * .1))
      
      a_P_1["Infection", "Infection", i, j] <- 1 - (HPV_Well_15to24[j] + v_p_mort_lessHPV[i] + 
                                                     ((1 - HPV_Well_15to24[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))) + (
        (1 - HPV_Well_15to24[j]) * ((1 - exp(-.2 * 3)) * .1)))

     }
 }

# The following enters all transition probabilities for ages 25-29 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
 for (i in 26:30) {
     for (j in 1:n.sims) {
      a_P_1["Infection", "Well", i, j] <- HPV_Well_25to29[j]
      
      a_P_1["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
      a_P_1["Infection", "LSIL", i, j] <- ((1 - HPV_Well_25to29[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1)))
      
      a_P_1["Infection", "HSIL", i, j] <- ((1 - HPV_Well_25to29[j]) * ((1 - exp(-.2 * 3)) * .1))
      
      a_P_1["Infection", "Infection", i, j] <- 1 - (HPV_Well_25to29[j] + v_p_mort_lessHPV[i] + 
                                                     ((1 - HPV_Well_25to29[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))) + (
        (1 - HPV_Well_25to29[j]) * ((1 - exp(-.2 * 3)) * .1)))

     }
 }

# The following enters all transition probabilities for ages 30-85 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
 for (i in 31:86) {
     for (j in 1:n.sims) {
      a_P_1["Infection", "Well", i, j] <- HPV_Well_30toEnd[j]
      
      a_P_1["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
      a_P_1["Infection", "LSIL", i, j] <- ((1 - HPV_Well_30toEnd[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1)))
      
      a_P_1["Infection", "HSIL", i, j] <- ((1 - HPV_Well_30toEnd[j]) * ((1 - exp(-.2 * 3)) * .1))
      
      a_P_1["Infection", "Infection", i, j] <- 1 - (HPV_Well_30toEnd[j] + v_p_mort_lessHPV[i] + 
                                                     ((1 - HPV_Well_30toEnd[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))) + (
        (1 - HPV_Well_30toEnd[j]) * ((1 - exp(-.2 * 3)) * .1)))

     }
 }

# ==========================================================================================
# New Treatment  ----------------------------------------------------------
# ==========================================================================================
# These are all the state transitions with added vaccine treatment. The following creates a
# transition array that is then used to simulate the Markov cohort with treatment.

# Transitions from Well State ---------------------------------------------

# The following enters all transition probabilities for ages 0-85 for each transition from
# the state Well, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death-less-HPV
# within each individual health state.
for (i in 1:86) {
 for (j in 1:n.sims) {
  a_P_2["Well", "Death", i, j] <-  v_p_mort_lessHPV[i]
  
  a_P_2["Well", "Infection", i, j] <- (omega.age[j, i] * (1 - pEfficacy.vac[j]))
  
  a_P_2["Well", "Well", i, j] <- 1 - (v_p_mort_lessHPV[i] + (
   omega.age[j, i] * (1 - pEfficacy.vac[j])))
  
 }
}

# Note: all transition probabilities for ages 0-14 are already filled by the natural structure
# of the transition array, i.e. all initial transitions are = 0 before being filled in by the
# proceeding code. This is why you can omit these transitions where applicable from the for 
# loop. Also note that this is why the mortality data for deaths from well state are the only
# data filled in for age groups 0-14.

# Transitions from Infection State ----------------------------------------

# The following enters all transition probabilities for ages 15-24 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
 for (i in 16:25) {
     for (j in 1:n.sims) {
      a_P_2["Infection", "Well", i, j] <- HPV_Well_15to24[j]
      
      a_P_2["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
      a_P_2["Infection", "LSIL", i, j] <- ((1 - HPV_Well_15to24[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1)))
      
      a_P_2["Infection", "HSIL", i, j] <- ((1 - HPV_Well_15to24[j]) * ((1 - exp(-.2 * 3)) * .1))
      
      a_P_2["Infection", "Infection", i, j] <- 1 - (HPV_Well_15to24[j] + v_p_mort_lessHPV[i] + 
                                                     ((1 - HPV_Well_15to24[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))) + (
        (1 - HPV_Well_15to24[j]) * ((1 - exp(-.2 * 3)) * .1)))

     }
 }

# The following enters all transition probabilities for ages 25-29 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
 for (i in 26:30) {
     for (j in 1:n.sims) {
      a_P_2["Infection", "Well", i, j] <- HPV_Well_25to29[j]
      
      a_P_2["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
      a_P_2["Infection", "LSIL", i, j] <- ((1 - HPV_Well_25to29[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1)))
      
      a_P_2["Infection", "HSIL", i, j] <- ((1 - HPV_Well_25to29[j]) * ((1 - exp(-.2 * 3)) * .1))
      
      a_P_2["Infection", "Infection", i, j] <- 1 - (HPV_Well_25to29[j] + v_p_mort_lessHPV[i] + 
                                                     ((1 - HPV_Well_25to29[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))) + (
        (1 - HPV_Well_25to29[j]) * ((1 - exp(-.2 * 3)) * .1)))

     }
 }

# The following enters all transition probabilities for ages 30-85 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
 for (i in 31:86) {
     for (j in 1:n.sims) {
      a_P_2["Infection", "Well", i, j] <- HPV_Well_30toEnd[j]
      
      a_P_2["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
      a_P_2["Infection", "LSIL", i, j] <- ((1 - HPV_Well_30toEnd[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1)))
      
      a_P_2["Infection", "HSIL", i, j] <- ((1 - HPV_Well_30toEnd[j]) * ((1 - exp(-.2 * 3)) * .1))
      
      a_P_2["Infection", "Infection", i, j] <- 1 - (HPV_Well_30toEnd[j] + v_p_mort_lessHPV[i] + 
                                                     ((1 - HPV_Well_30toEnd[j]) * (
       (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))) + (
        (1 - HPV_Well_30toEnd[j]) * ((1 - exp(-.2 * 3)) * .1)))

     }
 }
a_P_2[, , 22, ]
# Transitions from LSIL State ---------------------------------------------

# The following enters all transition probabilities for ages 15-34 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
for (i in 16:35) {
    for (j in 1:n.sims) {
     a_P_2["LSIL", "Well", i, j] <- (LSIL_15_34[i] * 0.9)
     a_P_2["LSIL", "Infection", i, j] <- (LSIL_15_34[i] - (LSIL_15_34[i] * 0.9))
     a_P_2["LSIL", "HSIL", i, j] <- (1 - LSIL_15_34[i]) * (1 - exp(-0.1 * 6))
     a_P_2["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     a_P_2["LSIL", "LSIL", i, j] <- (1 - LSIL_15_34[i]) - ((1 - LSIL_15_34[i]) * 
                                                         (1 - exp(-0.1 * 6)) + 
                                                         v_p_mort_lessHPV[i])
    }
}
# Check:
# a_P_2["LSIL", "LSIL", 22, ] + a_P_2["LSIL", "Death", 22, ] + a_P_2["LSIL", "HSIL", 22, ] + 
# a_P_2["LSIL", "Infection", 22, ] + a_P_2["LSIL", "Well", 22, ]

# The following enters all transition probabilities for ages ≥35 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.
for (i in 36:86) {
    for (j in 1:n.sims) {
     a_P_2["LSIL", "Well", i, j] <- (LSIL_35_End[i] * 0.9)
     a_P_2["LSIL", "Infection", i, j] <- (LSIL_35_End[i] - (LSIL_35_End[i] * 0.9))
     a_P_2["LSIL", "HSIL", i, j] <- (1 - LSIL_35_End[i]) * (1 - exp(-0.35 * 6))
     a_P_2["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     a_P_2["LSIL", "LSIL", i, j] <- (1 - LSIL_35_End[i]) - ((1 - LSIL_35_End[i]) * 
                                                         (1 - exp(-0.1 * 6)) + 
                                                         v_p_mort_lessHPV[i])
    }
}


# Transitions from HSIL State ---------------------------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state HSIl, across the appropriate time horizon i and all probabilistic 
# simulations j. Note: all state transitions are dependent on any regression or death less HPV
# from that occurs within the state.

for (i in 16:86) {
    for (j in 1:n.sims) {

    }
}


# HSIL to Normal:
HSIL_n * 0.5
# HSIL to LSIL:
(HSIL_n - (HSIL_n * 0.5))
# HSIL to Stage I Cancer:
HSILCANC <- ((1 - (HSIL_n + v_p_mort_lessHPV[22])) * (1 - exp(-0.4 * 10)))
# HSIL to Death:
HSILDEAD <- v_p_mort_lessHPV[22]
# HSIL to HSIL:
HSILHSIL <- 1 - (HSILCANC + HSILDEAD + HSIL_n)





# Code run time:
end_time <- Sys.time()

end_time - start_time
