# ==========================================================================================
# Evidence Synthesis and Transition Probability Matrix  ---------------------
# ==========================================================================================
# Note: this code below sits on top of the parameter Inputs.R file, which is the source code for 
# all parameter data (04_parameter_inputs.R).

source("R/04_parameter_inputs.R")

# Several distributions have had to be truncated at an upper bound of 1 - (max mortality prob). 
# We found that the mortality rate when it included HIV/AIDS related mortality was extremely high 
# given the context of South Africa in 2009. In order to account for this, it would have be 
# desirable to incorporate another parameter to provide a conditional annual risk of dying 
# from HIV/AIDS. However, because this was purely a replication exercise, we did not 
# find it necessary and focused on the outcomes of the comparative models.

# ==========================================================================================
# Evidence Synthesis Model ------------------------------------------------
# ==========================================================================================
model_String <- "
model {

### SUB-MODEL 1: AGE-SPECIFIC PREVALENCE
# Model parameters abbreviated by .age. Note: this is equivalent to Monte Carlo PSA, 
# as it is technically sampling directly from a prior and is not propogated into a 
# posterior using a likelihood model. However, a hyperprior is used for the population 
# variance to account for a greater uncertainty in each age-population.
  for (i in 1:86) {
    # Monte Carlo:
     omega.age[i] ~ dlnorm(mu.a.log[i], prec.age[i])
    
    # Note use of pow() function, using -2 is a shorthand inverse method equivalent 
    # to (1 / x^2):
     log(prec.age[i]) <- pow(sigma.age[i], -2)
     
    # Relatively wide prior on variance for each age group.
     sigma.age[i] ~ dunif(0, 10)
  }

### END OF SUB-MODEL 1.

### SUB-MODEL 2: POP. LEVEL VACCINE-EFFICACY. 
# Note: this is a fully integrated Bayesian model, as it combines evidence directly via the 
# likelihood and a prior. Model parameters abbreviated by .vac.
  for (i in 1:Nstud.vac) {
    # Likelihood:
     rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
     rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])

    # Logistic link function:
     logit(pA.vac[i]) <- mu.vac[i]
     logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
    
    # Average effect prior for SUB-MODEL 2:
     mu.vac[i] ~ dnorm(0, 1e-6)
    # Prior for sub-model 2 (Random. pop. effect):
     delta.vac[i] ~ dt(psi.vac, prec.vac, 1) # if desired can be ~ dnorm(psi.vac, prec.vac)
    
     ### Mixed predictive check for SUB-MODEL 2:
       # Predictive likelihood:
        rA.mxd[i] ~ dbin(pA.new[i], nA.vac[i])
        
       # Predictive logit link function:
        logit(pA.new[i]) <- mu.vac[i] + delta.new
               
       # Mixed predictve p-value:
        pA.mxd[i] <- step(rA.mxd[i] - rA.vac[i]) - 0.5 * equals(rA.mxd[i], rA.vac[i])

  }
  
   # Hyperpriors for SUB-MODEL 2:
    psi.vac ~ dnorm(0, 1.0e-6)
    prec.vac <- pow(tau.vac, -2)
    tau.vac ~  dunif(0, 10)
  
  # Transformations for SUB-MODEL 2:
   # Convert LOR to OR
    OR.vac <- exp(psi.vac)
   # Convert OR to probability for vaccine efficacy
    pEfficacy.vac <- 1 / (1 + OR.vac)
   
     # Predicted average 
     # treatment effect:
      delta.new ~ dnorm(psi.vac, prec.vac)

### END OF SUB-MODEL 2.

### SUB-MODEL 3: CANCER PROGRESSION AND 5-YEAR SURVIVAL STAGES I-IV.
# Model parameters abbreviated by .canc. Note: this is equivalent to a standard 
# Monte Carlo PSA, as it is technically sampling directly from a prior and it is
# *not* propogated into a posterior using a likelihood model. We have had to truncate these 
# distributions in order for the ASSA mortality data to be properly combined and the state
# transition probabilities to be proper.

   # Distribution according to Stage:
    Stage.I.canc ~ dbeta(alpha.StageI, beta.StageI)T(0, 0.88)
    Stage.II.canc ~ dbeta(alpha.StageII, beta.StageII)T(0, 0.88)
    Stage.III.canc ~ dbeta(alpha.StageIII, beta.StageIII)T(0, 0.88)
    StageIV.Detected ~ dbeta(alpha.StageIV, beta.StageIV)T(0, 0.88)
    
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
   
### END OF SUB-MODEL 3.

### SUB-MODEL 4: INFECTION PROGRESSION:
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

### END OF SUB-MODEL 4.

### SUB-MODEL 5: LSIL & HSIL PROGRESSION:
# Note: this is equivalent to a standard Monte Carlo PSA, as it is technically sampling
# directly from a prior and it is *not* propogated into a posterior using a likelihood 
# model. I have had to truncate these distributions in order for the ASSA mortality data 
# to be properly combined and for the state transition probabilities to be proper.
   LSIL_15_34 ~ dbeta(alpha.LSIL_15to34, beta.LSIL_15to34)T(0, 0.88)
   LSIL_35_85 ~ dbeta(alpha.LSIL_35to85, beta.LSIL_35to85)T(0, 0.88)
   HSIL_n ~ dbeta(alpha.HSIL, beta.HSIL)T(0, 0.88)

### END OF SUB-MODEL 5.

 }
"
writeLines(text = model_String, con = "data/jags_model.txt")

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
  alpha.LSIL_35to85 = alpha.LSIL_35to85, beta.LSIL_35to85 = beta.LSIL_35to85,
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
  "OR.vac", "pEfficacy.vac", "delta.vac", "psi.vac",
  # Well to infection prevalence:
  "omega.age",
  # HPV Progression:
  "HPV_Well_15to24", "HPV_Well_25to29",
  "HPV_Well_30toEnd",
  # LSIL and HSIL Progression:
  "LSIL_15_34", "LSIL_35_85",
  "HSIL_n",
   # Cancer Stage Progression:
  "Stage.I.canc", "Stage.II.canc",
  "Stage.III.canc", "StageIV.Detected",
  # Cancer Stage Survival progression:
  "surv.StageI_year1", "surv.StageI_year2", "surv.StageI_year3", 
  "surv.StageI_year4", "surv.StageI_year5",
  "surv.StageII_year1", "surv.StageII_year2", "surv.StageII_year3", 
  "surv.StageII_year4", "surv.StageII_year5",
  "surv.StageIII_year1", "surv.StageIII_year2", "surv.StageIII_year3",
  "surv.StageIII_year4", "surv.StageIII_year5",
  "surv.StageIV_year1", "surv.StageIV_year2", "surv.StageIV_year3",
  "surv.StageIV_year4", "surv.StageIV_year5",
  
  # Mixed predictive check. Provides a conservative p-value measure for 
  # inconsistency of evidence across each study comparative to the rest:
  "pA.mxd"
  )

# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 45000
n.burnin <- 10000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "data/jags_model.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS
# One can automate convergence. However, the improvement is negligible, and it increased the run
# time of the model from < 60secs to > 2 mins. However, if desired, uncomment the line of code 
# below to add auto chain convergence:
# mod_JAGS <- autojags(mod_JAGS, Rhat = 1.01)

# Posterior inspection:
# posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
# dimnames(x = posterior)

# color_scheme_set("mix-gray-brightblue")
# mcmc_combo(posterior, pars = c("pEfficacy.vac", "LSIL_15_34", "HSIL_n", "Stage.III.canc"))

# mcmc_dens_chains_data(posterior, pars = c("pA.mxd[1]", "pA.mxd[2]", "pA.mxd[3]", "pA.mxd[4]",
#                                "pA.mxd[5]", "pA.mxd[6]", "pA.mxd[7]", "pA.mxd[8]",
#                                "pA.mxd[9]", "pA.mxd[10]", "pA.mxd[11]"))

# Attach JAGS model to envir.
attach.jags(mod_JAGS)

# ==========================================================================================
# Probability Matrix --------------------------------------------------
# ==========================================================================================
 # The following creates a transition array that is then used to simulate the Markov cohort for both
# treatments, Status Quo and New Treatment.

# The following sets up the structure of the transitions between Health States:
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

# Transition array for each alternative:
a_P_1 <- array(0, dim = c(n_states, n_states, n_t + 1, n.sims),
             dimnames = list(v_n, v_n, 0:(n_t), 1:n.sims))
a_P_2 <- array(0, dim = c(n_states, n_states, n_t + 1, n.sims),
             dimnames = list(v_n, v_n, 0:(n_t), 1:n.sims))

# Note: all progression probabilities are assumed dependent on the probability of regression. Moreover, 
# due to the natural structure of the transition array, there is no need to fill in probabilities = 0, 
# such as the transitions from Infection to LSIL for ages ≤ 14.

# ==========================================================================================
# Status-Quo  ----------------------------------------------------------
# ==========================================================================================
# These are for all the state transitions with no vaccine treatment.

# Transitions from Well State ---------------------------------------------

# The following enters all transition probabilities for ages 0-85 for each transition from
# the state Well, across the appropriate time horizon i and all probabilistic 
# simulations j
for (i in 1:n_t) {
 for (j in 1:n.sims) {
  a_P_1["Well", "Death", i, j] <-  v_p_mort_lessHPV[i]
  
  a_P_1["Well", "Infection", i, j] <- omega.age[j, i]
  
  a_P_1["Well", "Well", i, j] <- 1 - (v_p_mort_lessHPV[i] + omega.age[j, i])
 }
}

# Transitions from Infection State ----------------------------------------

# The following enters all transition probabilities for ages 15-24 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 16:25) {
    for (j in 1:n.sims) {
     a_P_1["Infection", "Well", i, j] <- HPV_Well_15to24[j]
      
     a_P_1["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
     a_P_1["Infection", "LSIL", i, j] <- (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.9)
                                           
     a_P_1["Infection", "HSIL", i, j] <- (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.1)
      
     a_P_1["Infection", "Infection", i, j] <- 1 - (HPV_Well_15to24[j] + v_p_mort_lessHPV[i] + (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.9) + (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.1))
    }
}
# The following enters all transition probabilities for ages 25-29 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 26:30) {
    for (j in 1:n.sims) {
     a_P_1["Infection", "Well", i, j] <- HPV_Well_25to29[j]
      
     a_P_1["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
     a_P_1["Infection", "LSIL", i, j] <- (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.9)
                                           
     a_P_1["Infection", "HSIL", i, j] <- (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.1)
      
     a_P_1["Infection", "Infection", i, j] <- 1 - (HPV_Well_25to29[j] + v_p_mort_lessHPV[i] + (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.9) + (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.1))
    }
}
# The following enters all transition probabilities for ages 30-85 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 31:86) {
    for (j in 1:n.sims) {
     a_P_1["Infection", "Well", i, j] <- HPV_Well_30toEnd[j]
      
     a_P_1["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
     a_P_1["Infection", "LSIL", i, j] <- (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 3))) * 0.9)
                                           
     a_P_1["Infection", "HSIL", i, j] <- (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 36))) * 0.1)
      
     a_P_1["Infection", "Infection", i, j] <- 1 - (HPV_Well_30toEnd[j] + v_p_mort_lessHPV[i] + (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 3))) * 0.9) + (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 36))) * 0.1))
    }
}

# Transitions from LSIL State ---------------------------------------------

# The following enters all transition probabilities for ages 15-34 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 16:35) {
    for (j in 1:n.sims) {
     a_P_1["LSIL", "Well", i, j] <- (LSIL_15_34[j] * 0.9)
     
     a_P_1["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["LSIL", "Infection", i, j] <- (LSIL_15_34[j] * 0.1)
     
     a_P_1["LSIL", "HSIL", i, j] <- ((1 - (LSIL_15_34[j] + v_p_mort_lessHPV[i])) * (1 - exp(-0.1 * 6)))
     
     a_P_1["LSIL", "LSIL", i, j] <- 1 - ((LSIL_15_34[j] * 0.9) + (LSIL_15_34[j] * 0.1) + v_p_mort_lessHPV[i] + ((1 - (LSIL_15_34[j] + v_p_mort_lessHPV[i])) * (1 - exp(-0.1 * 6))))
    }
}

# The following enters all transition probabilities for ages ≥25 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 36:86) {
    for (j in 1:n.sims) {
     a_P_1["LSIL", "Well", i, j] <- (LSIL_35_85[j] * 0.9)
     
     a_P_1["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["LSIL", "Infection", i, j] <- (LSIL_35_85[j] * 0.1)
     
     a_P_1["LSIL", "HSIL", i, j] <- ((1 - LSIL_35_85[j]) * (1 - exp(-0.35 * 6)))
     
     a_P_1["LSIL", "LSIL", i, j] <- 1 - ((LSIL_35_85[j] * 0.9) + (LSIL_35_85[j] * 0.1) + v_p_mort_lessHPV[i] + ((1 - LSIL_35_85[j]) * (1 - exp(-0.35 * 6))))
    }
}

# Transitions from HSIL State ---------------------------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state HSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["HSIL", "Well", i, j] <- (HSIL_n[j] * 0.5)
     
     a_P_1["HSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["HSIL", "LSIL", i, j] <- (HSIL_n[j] * 0.5)
     
     a_P_1["HSIL", "Stage-I Cancer", i, j] <- ((1 - ((HSIL_n[j] * 0.5) + (HSIL_n[j] * 0.5))) * (1 - exp(-0.4 * 10)))
     
     a_P_1["HSIL", "HSIL", i, j] <- 1 - ((HSIL_n[j] * 0.5) + (HSIL_n[j] * 0.5) + v_p_mort_lessHPV[i] + ((1 - ((HSIL_n[j] * 0.5) + (HSIL_n[j] * 0.5))) * (1 - exp(-0.4 * 10))))
    }
}

# Transitions from Stage-I Cervix Cancer State ----------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-I Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Stage-I Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["Stage-I Cancer", "Stage-II Cancer", i, j] <- Stage.I.canc[j]
     
     a_P_1["Stage-I Cancer", "Detected.Stage-I Year 1", i, j] <- ((1 - Stage.I.canc[j]) * 0.15)
     
     a_P_1["Stage-I Cancer", "Stage-I Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + Stage.I.canc[j] + ((1 - Stage.I.canc[j]) * 0.15))
    }
}

# Transitions from Stage-II Cervix Cancer State ---------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-II Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Stage-II Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["Stage-II Cancer", "Stage-III Cancer", i, j] <- Stage.II.canc[j]
     
     a_P_1["Stage-II Cancer", "Detected.Stage-II Year 1", i, j] <- ((1 - Stage.II.canc[j]) * 0.225)
     
     a_P_1["Stage-II Cancer", "Stage-II Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + Stage.II.canc[j] + ((1 - Stage.II.canc[j]) * 0.225))
    }
}

# Transitions from Stage-III Cervix Cancer State --------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-III Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Stage-III Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["Stage-III Cancer", "Stage-IV Cancer", i, j] <- Stage.III.canc[j]
     
     a_P_1["Stage-III Cancer", "Detected.Stage-III Year 1", i, j] <- ((1 - Stage.III.canc[j]) * 0.6)
     
     a_P_1["Stage-III Cancer", "Stage-III Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + Stage.III.canc[j] + ((1 - Stage.III.canc[j]) * 0.6))
    }
}

# Transitions from Stage-IV Cervix Cancer State ---------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-IV Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Stage-IV Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_1["Stage-IV Cancer", "Detected.Stage-IV Year 1", i, j] <- StageIV.Detected[j]
     
     a_P_1["Stage-IV Cancer", "Stage-IV Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + StageIV.Detected[j])
    }
}

# Transitions from All Stage-I Cancer Survivor/Detected States ---------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-I Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Detected.Stage-I Year 1", "Detected.Stage-I Year 2", i, j] <- surv.StageI_year1[j]
     
     a_P_1["Detected.Stage-I Year 1", "Death", i, j] <- (1 - surv.StageI_year1[j])
     
     a_P_1["Detected.Stage-I Year 2", "Detected.Stage-I Year 3", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j])
     
     a_P_1["Detected.Stage-I Year 2", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j]))
     
     a_P_1["Detected.Stage-I Year 3", "Detected.Stage-I Year 4", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j])
     
     a_P_1["Detected.Stage-I Year 3", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j]))

     a_P_1["Detected.Stage-I Year 4", "Detected.Stage-I Year 5", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j])
     
     a_P_1["Detected.Stage-I Year 4", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j]))
     
     a_P_1["Detected.Stage-I Year 5", "Cancer Survivor", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j] * surv.StageI_year5[j])
     
     a_P_1["Detected.Stage-I Year 5", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j] * surv.StageI_year5[j]))
     }
}

# Transitions from All Stage-II Cancer Survivor/Detected States -----------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-II Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Detected.Stage-II Year 1", "Detected.Stage-II Year 2", i, j] <- surv.StageII_year1[j]
     
     a_P_1["Detected.Stage-II Year 1", "Death", i, j] <- (1 - surv.StageII_year1[j])
     
     a_P_1["Detected.Stage-II Year 2", "Detected.Stage-II Year 3", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j])
     
     a_P_1["Detected.Stage-II Year 2", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j]))
     
     a_P_1["Detected.Stage-II Year 3", "Detected.Stage-II Year 4", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j])
     
     a_P_1["Detected.Stage-II Year 3", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j]))

     a_P_1["Detected.Stage-II Year 4", "Detected.Stage-II Year 5", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j])
     
     a_P_1["Detected.Stage-II Year 4", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j]))
     
     a_P_1["Detected.Stage-II Year 5", "Cancer Survivor", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j] * surv.StageII_year5[j])
     
     a_P_1["Detected.Stage-II Year 5", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j] * surv.StageII_year5[j]))
    }
}

# Transitions from All Stage-III Cancer Survivor/Detected States ----------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-III Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Detected.Stage-III Year 1", "Detected.Stage-III Year 2", i, j] <- surv.StageIII_year1[j]
     
     a_P_1["Detected.Stage-III Year 1", "Death", i, j] <- (1 - surv.StageIII_year1[j])
     
     a_P_1["Detected.Stage-III Year 2", "Detected.Stage-III Year 3", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j])
     
     a_P_1["Detected.Stage-III Year 2", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j]))
     
     a_P_1["Detected.Stage-III Year 3", "Detected.Stage-III Year 4", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j])
     
     a_P_1["Detected.Stage-III Year 3", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j]))

     a_P_1["Detected.Stage-III Year 4", "Detected.Stage-III Year 5", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j])
     
     a_P_1["Detected.Stage-III Year 4", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j]))
     
     a_P_1["Detected.Stage-III Year 5", "Cancer Survivor", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j] * surv.StageIII_year5[j])
     
     a_P_1["Detected.Stage-III Year 5", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j] * surv.StageIII_year5[j]))
    }
}

# Transitions from All Stage-IV Cancer Survivor/Detected States -----------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-IV Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Detected.Stage-IV Year 1", "Detected.Stage-IV Year 2", i, j] <- surv.StageIV_year1[j]
     
     a_P_1["Detected.Stage-IV Year 1", "Death", i, j] <- (1 - surv.StageIV_year1[j])
     
     a_P_1["Detected.Stage-IV Year 2", "Detected.Stage-IV Year 3", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j])
     
     a_P_1["Detected.Stage-IV Year 2", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j]))
     
     a_P_1["Detected.Stage-IV Year 3", "Detected.Stage-IV Year 4", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j])
     
     a_P_1["Detected.Stage-IV Year 3", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j]))

     a_P_1["Detected.Stage-IV Year 4", "Detected.Stage-IV Year 5", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j])
     
     a_P_1["Detected.Stage-IV Year 4", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j]))
     
     a_P_1["Detected.Stage-IV Year 5", "Cancer Survivor", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j] * surv.StageIV_year5[j])
     
     a_P_1["Detected.Stage-IV Year 5", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j] * surv.StageIV_year5[j]))
    }
}

# Transitions from Cancer Survivor State ----------------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state of Cancer Survivor, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_1["Cancer Survivor", "Cancer Survivor", i, j] <- (1 - v_p_mort_lessHPV[i])
     
     a_P_1["Cancer Survivor", "Death", i, j] <- v_p_mort_lessHPV[i]
    }
}

# Transitions from Death State --------------------------------------------

# The following enters all transition probabilities for ages 0-85 for each transition from
# the state of Death, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 0:n_t) {
 for (j in 1:n.sims) {
  a_P_1["Death", "Death", i, j] <- 1
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
# simulations j
for (i in 0:n_t) {
 for (j in 1:n.sims) {
  a_P_2["Well", "Death", i, j] <-  v_p_mort_lessHPV[i]
  
  a_P_2["Well", "Infection", i, j] <- (1 - pEfficacy.vac[j]) * omega.age[j, i]
  
  a_P_2["Well", "Well", i, j] <- 1 - (v_p_mort_lessHPV[i] + ((1 - pEfficacy.vac[j]) * omega.age[j, i]))
 } 
}

# Transitions from Infection State ----------------------------------------

# The following enters all transition probabilities for ages 15-24 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 16:25) {
    for (j in 1:n.sims) {
     a_P_2["Infection", "Well", i, j] <- HPV_Well_15to24[j]
 
     a_P_2["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
     a_P_2["Infection", "LSIL", i, j] <- (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.9)
                                           
     a_P_2["Infection", "HSIL", i, j] <- (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.1)
      
     a_P_2["Infection", "Infection", i, j] <- 1 - (HPV_Well_15to24[j] + v_p_mort_lessHPV[i] + (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.9) + (((1 - HPV_Well_15to24[j]) * (1 - exp(-.2 * 3))) * 0.1))
    }
}
# The following enters all transition probabilities for ages 25-29 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 26:30) {
    for (j in 1:n.sims) {
     a_P_2["Infection", "Well", i, j] <- HPV_Well_25to29[j]
      
     a_P_2["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
     a_P_2["Infection", "LSIL", i, j] <- (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.9)
                                           
     a_P_2["Infection", "HSIL", i, j] <- (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.1)
      
     a_P_2["Infection", "Infection", i, j] <- 1 - (HPV_Well_25to29[j] + v_p_mort_lessHPV[i] + (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.9) + (((1 - HPV_Well_25to29[j]) * (1 - exp(-.2 * 3))) * 0.1))
    }
}
# The following enters all transition probabilities for ages 30-85 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 31:86) {
    for (j in 1:n.sims) {
     a_P_2["Infection", "Well", i, j] <- HPV_Well_30toEnd[j]
      
     a_P_2["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
      
     a_P_2["Infection", "LSIL", i, j] <- (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 3))) * 0.9)
                                           
     a_P_2["Infection", "HSIL", i, j] <- (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 36))) * 0.1)
      
     a_P_2["Infection", "Infection", i, j] <- 1 - (HPV_Well_30toEnd[j] + v_p_mort_lessHPV[i] + (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 3))) * 0.9) + (((1 - HPV_Well_30toEnd[j]) * (1 - exp(-.2 * 36))) * 0.1))
    }
}

# Transitions from LSIL State ---------------------------------------------

# The following enters all transition probabilities for ages 15-34 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 16:35) {
    for (j in 1:n.sims) {
     a_P_2["LSIL", "Well", i, j] <- (LSIL_15_34[j] * 0.9)
     
     a_P_2["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["LSIL", "Infection", i, j] <- (LSIL_15_34[j] * 0.1)
     
     a_P_2["LSIL", "HSIL", i, j] <- ((1 - (LSIL_15_34[j] + v_p_mort_lessHPV[i])) * (1 - exp(-0.1 * 6)))
     
     a_P_2["LSIL", "LSIL", i, j] <- 1 - ((LSIL_15_34[j] * 0.9) + (LSIL_15_34[j] * 0.1) + v_p_mort_lessHPV[i] + ((1 - (LSIL_15_34[j] + v_p_mort_lessHPV[i])) * (1 - exp(-0.1 * 6))))
    }
}

# The following enters all transition probabilities for ages ≥25 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 36:86) {
    for (j in 1:n.sims) {
     a_P_2["LSIL", "Well", i, j] <- (LSIL_35_85[j] * 0.9)
     
     a_P_2["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["LSIL", "Infection", i, j] <- (LSIL_35_85[j] * 0.1)
     
     a_P_2["LSIL", "HSIL", i, j] <- ((1 - LSIL_35_85[j]) * (1 - exp(-0.35 * 6)))
     
     a_P_2["LSIL", "LSIL", i, j] <- 1 - ((LSIL_35_85[j] * 0.9) + (LSIL_35_85[j] * 0.1) + v_p_mort_lessHPV[i] + ((1 - LSIL_35_85[j]) * (1 - exp(-0.35 * 6))))
    }
}

# Transitions from HSIL State ---------------------------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state HSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["HSIL", "Well", i, j] <- (HSIL_n[j] * 0.5)
     
     a_P_2["HSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["HSIL", "LSIL", i, j] <- (HSIL_n[j] * 0.5)
     
     a_P_2["HSIL", "Stage-I Cancer", i, j] <- ((1 - ((HSIL_n[j] * 0.5) + (HSIL_n[j] * 0.5))) * (1 - exp(-0.4 * 10)))
     
     a_P_2["HSIL", "HSIL", i, j] <- 1 - ((HSIL_n[j] * 0.5) + (HSIL_n[j] * 0.5) + v_p_mort_lessHPV[i] + ((1 - ((HSIL_n[j] * 0.5) + (HSIL_n[j] * 0.5))) * (1 - exp(-0.4 * 10))))
    }
}

# Transitions from Stage-I Cervix Cancer State ----------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-I Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Stage-I Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["Stage-I Cancer", "Stage-II Cancer", i, j] <- Stage.I.canc[j]
     
     a_P_2["Stage-I Cancer", "Detected.Stage-I Year 1", i, j] <- ((1 - Stage.I.canc[j]) * 0.15)
     
     a_P_2["Stage-I Cancer", "Stage-I Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + Stage.I.canc[j] + ((1 - Stage.I.canc[j]) * 0.15))
    }
}

# Transitions from Stage-II Cervix Cancer State ---------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-II Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Stage-II Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["Stage-II Cancer", "Stage-III Cancer", i, j] <- Stage.II.canc[j]
     
     a_P_2["Stage-II Cancer", "Detected.Stage-II Year 1", i, j] <- ((1 - Stage.II.canc[j]) * 0.225)
     
     a_P_2["Stage-II Cancer", "Stage-II Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + Stage.II.canc[j] + ((1 - Stage.II.canc[j]) * 0.225))
    }
}

# Transitions from Stage-III Cervix Cancer State --------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-III Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Stage-III Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["Stage-III Cancer", "Stage-IV Cancer", i, j] <- Stage.III.canc[j]
     
     a_P_2["Stage-III Cancer", "Detected.Stage-III Year 1", i, j] <- ((1 - Stage.III.canc[j]) * 0.6)
     
     a_P_2["Stage-III Cancer", "Stage-III Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + Stage.III.canc[j] + ((1 - Stage.III.canc[j]) * 0.6))
    }
}

# Transitions from Stage-IV Cervix Cancer State ---------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-IV Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Stage-IV Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_2["Stage-IV Cancer", "Detected.Stage-IV Year 1", i, j] <- StageIV.Detected[j]
     
     a_P_2["Stage-IV Cancer", "Stage-IV Cancer", i, j] <- 1 - (v_p_mort_lessHPV[i] + StageIV.Detected[j])
    }
}

# Transitions from All Stage-I Cancer Survivor/Detected States ---------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-I Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Detected.Stage-I Year 1", "Detected.Stage-I Year 2", i, j] <- surv.StageI_year1[j]
     
     a_P_2["Detected.Stage-I Year 1", "Death", i, j] <- (1 - surv.StageI_year1[j])
     
     a_P_2["Detected.Stage-I Year 2", "Detected.Stage-I Year 3", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j])
     
     a_P_2["Detected.Stage-I Year 2", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j]))
     
     a_P_2["Detected.Stage-I Year 3", "Detected.Stage-I Year 4", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j])
     
     a_P_2["Detected.Stage-I Year 3", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j]))

     a_P_2["Detected.Stage-I Year 4", "Detected.Stage-I Year 5", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j])
     
     a_P_2["Detected.Stage-I Year 4", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j]))
     
     a_P_2["Detected.Stage-I Year 5", "Cancer Survivor", i, j] <- (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j] * surv.StageI_year5[j])
     
     a_P_2["Detected.Stage-I Year 5", "Death", i, j] <- (1 - (surv.StageI_year1[j] * surv.StageI_year2[j] * surv.StageI_year3[j] * surv.StageI_year4[j] * surv.StageI_year5[j]))
     }
}

# Transitions from All Stage-II Cancer Survivor/Detected States -----------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-II Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Detected.Stage-II Year 1", "Detected.Stage-II Year 2", i, j] <- surv.StageII_year1[j]
     
     a_P_2["Detected.Stage-II Year 1", "Death", i, j] <- (1 - surv.StageII_year1[j])
     
     a_P_2["Detected.Stage-II Year 2", "Detected.Stage-II Year 3", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j])
     
     a_P_2["Detected.Stage-II Year 2", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j]))
     
     a_P_2["Detected.Stage-II Year 3", "Detected.Stage-II Year 4", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j])
     
     a_P_2["Detected.Stage-II Year 3", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j]))

     a_P_2["Detected.Stage-II Year 4", "Detected.Stage-II Year 5", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j])
     
     a_P_2["Detected.Stage-II Year 4", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j]))
     
     a_P_2["Detected.Stage-II Year 5", "Cancer Survivor", i, j] <- (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j] * surv.StageII_year5[j])
     
     a_P_2["Detected.Stage-II Year 5", "Death", i, j] <- (1 - (surv.StageII_year1[j] * surv.StageII_year2[j] * surv.StageII_year3[j] * surv.StageII_year4[j] * surv.StageII_year5[j]))
    }
}

# Transitions from All Stage-III Cancer Survivor/Detected States ----------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-III Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Detected.Stage-III Year 1", "Detected.Stage-III Year 2", i, j] <- surv.StageIII_year1[j]
     
     a_P_2["Detected.Stage-III Year 1", "Death", i, j] <- (1 - surv.StageIII_year1[j])
     
     a_P_2["Detected.Stage-III Year 2", "Detected.Stage-III Year 3", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j])
     
     a_P_2["Detected.Stage-III Year 2", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j]))
     
     a_P_2["Detected.Stage-III Year 3", "Detected.Stage-III Year 4", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j])
     
     a_P_2["Detected.Stage-III Year 3", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j]))

     a_P_2["Detected.Stage-III Year 4", "Detected.Stage-III Year 5", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j])
     
     a_P_2["Detected.Stage-III Year 4", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j]))
     
     a_P_2["Detected.Stage-III Year 5", "Cancer Survivor", i, j] <- (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j] * surv.StageIII_year5[j])
     
     a_P_2["Detected.Stage-III Year 5", "Death", i, j] <- (1 - (surv.StageIII_year1[j] * surv.StageIII_year2[j] * surv.StageIII_year3[j] * surv.StageIII_year4[j] * surv.StageIII_year5[j]))
    }
}

# Transitions from All Stage-IV Cancer Survivor/Detected States -----------

# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-IV Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Detected.Stage-IV Year 1", "Detected.Stage-IV Year 2", i, j] <- surv.StageIV_year1[j]
     
     a_P_2["Detected.Stage-IV Year 1", "Death", i, j] <- (1 - surv.StageIV_year1[j])
     
     a_P_2["Detected.Stage-IV Year 2", "Detected.Stage-IV Year 3", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j])
     
     a_P_2["Detected.Stage-IV Year 2", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j]))
     
     a_P_2["Detected.Stage-IV Year 3", "Detected.Stage-IV Year 4", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j])
     
     a_P_2["Detected.Stage-IV Year 3", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j]))

     a_P_2["Detected.Stage-IV Year 4", "Detected.Stage-IV Year 5", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j])
     
     a_P_2["Detected.Stage-IV Year 4", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j]))
     
     a_P_2["Detected.Stage-IV Year 5", "Cancer Survivor", i, j] <- (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j] * surv.StageIV_year5[j])
     
     a_P_2["Detected.Stage-IV Year 5", "Death", i, j] <- (1 - (surv.StageIV_year1[j] * surv.StageIV_year2[j] * surv.StageIV_year3[j] * surv.StageIV_year4[j] * surv.StageIV_year5[j]))
    }
}

# Transitions from Cancer Survivor State ----------------------------------

# The following enters all transition probabilities for ages 15-85 for each transition from
# the state of Cancer Survivor, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 16:86) {
    for (j in 1:n.sims) {
     a_P_2["Cancer Survivor", "Cancer Survivor", i, j] <- (1 - v_p_mort_lessHPV[i])
     
     a_P_2["Cancer Survivor", "Death", i, j] <- v_p_mort_lessHPV[i]
    }
}

# Transitions from Death State --------------------------------------------

# The following enters all transition probabilities for ages 0-85 for each transition from
# the state of Death, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 0:n_t) {
 for (j in 1:n.sims) {
  a_P_2["Death", "Death", i, j] <- 1
  }
}

# End file ----------------------------------------------------------------