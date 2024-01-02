# ==========================================================================================
# Evidence Synthesis and Transition Probability Matrix  ---------------------
# ==========================================================================================
# Note: this code below sits on top of the parameter Inputs.R file, which is the source code
# for all parameter data (04_parameter_inputs.R).

source("R/04_parameter_inputs.R")

# We found that the mortality rate when it included HIV/AIDS related mortality was extremely
# high given the context of South Africa in 2009. In order to account for this, it would be
# have been desirable to incorporate a parameter which accounts for the conditional annual 
# risk of dying from HIV/AIDS. However, because this was purely a replication exercise, we
# did not and focused rather on the outcomes of the comparative models.

# ==========================================================================================
# Evidence Synthesis Model ------------------------------------------------
# ==========================================================================================
model_String <- "
model {
### SUB-MODEL 1: AGE-SPECIFIC PREVALENCE
# Model parameters abbreviated by .age.:
 # Incidence for ages 0-14:
  for (i in 1:15) {
    # Assumption of 0 incidence from ages 0-14
     p.age[i] <- 0
     
  }
 # Incidence for ages 15-85:
  for (i in 16:86) {
    # Likelihood for incidence rate:
     rate_Incidence[i] ~ dexp(lambda.age[i])

    # Prior on rate with separate variance within each age-group:
     lambda.age[i] ~ dlnorm(0, prec.age[i]) # the uncertain prior  pulls likelihood
                                              # due to uncertainty of US data
                                              # applied to SA context.
     # Independent variance to precision within each age-group:
     prec.age[i] <- pow(tau.age[i], -2)
     # Stdev for each age-group:
     tau.age[i] ~  dunif(0, 10) # with uniform prior bounded at U(a = 0, b = 10)
     
    # Expected rate for each age group per cycle:
     mu.age[i] <- 1 / lambda.age[i] # reparametrise as scale parameter, i.e., 
                                    # 1/lambda = E[lambda]

    # Transformation of mean rate to annual probability:
     p.age[i] <- 1 - exp(-mu.age[i] * 1)
  
  }

    # Hyperprior for SUB-MODEL 1:
     # prec.age <- pow(tau.age, -2) # inverse power transformation

### END OF SUB-MODEL 1.

### SUB-MODEL 2: POP. LEVEL VACCINE-EFFICACY. 
# Note: this is a Bayesian Posterior model, as it combines evidence directly via the 
# likelihood and prior. Model parameters are abbreviated by .vac.
  for (i in 1:Nstud.vac) {
    # Binomial Likelihood:
     rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
     rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])

    # Logistic link function:
     logit(pA.vac[i]) <- mu.vac[i]
     logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
    
    # Average effect prior for SUB-MODEL 2:
     mu.vac[i] ~ dnorm(0, 1e-4) # average effect across groups.
    # Prior for sub-model 2 (Random. pop. effect):
     delta.vac[i] ~ dt(psi.vac, prec.vac, 1)  # with separate random effect between 
                                              # groups.
    
     ### Mixed predictive check for SUB-MODEL 2:
       # Predictive likelihood:
        rA.mxd[i] ~ dbin(pA.new[i], nA.vac[i])
        
       # Predictive logit link function:
        logit(pA.new[i]) <- mu.vac[i] + delta.new
               
       # Mixed predictve p-value:
        pA.mxd[i] <- step(rA.mxd[i] - rA.vac[i]) - 0.5 * equals(rA.mxd[i], rA.vac[i])

  }
  
  # Hyperpriors for SUB-MODEL 2:
     psi.vac ~ dnorm(0, 1.0e-4)
     prec.vac <- pow(tau.vac, -2)
     tau.vac ~  dt(0, 1, 1)T(0, )
  
  # Transformations for SUB-MODEL 2:
   # Convert LOR to OR
    OR.vac <- exp(psi.vac)
   # Convert OR to probability for vaccine efficacy
    pEfficacy.vac <- 1 / (1 + OR.vac)
   
     # Predicted average treatment effect:
      delta.new ~ dnorm(psi.vac, prec.vac)

### END OF SUB-MODEL 2.

### SUB-MODEL 3: CANCER PROGRESSION AND 5-YEAR SURVIVAL STAGES I-IV.
# Model parameters abbreviated by .canc. Note: this is equivalent to standard PSA, as it 
# is technically sampling directly from an integrated probability distribution and is 
# not explicitly propogated into a posterior via Bayes' theorem:

   # Distribution according to Stage:
    Stage.I.canc ~ dbeta(alpha.StageI, beta.StageI)
    Stage.II.canc ~ dbeta(alpha.StageII, beta.StageII)
    Stage.III.canc ~ dbeta(alpha.StageIII, beta.StageIII)
    StageIV.Detected ~ dbeta(alpha.StageIV, beta.StageIV)
    
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
# Note: this is equivalent to standard PSA, as it is technically sampling directly
# from an integrated probability distribution and is not explicitly propogated into a 
# posterior via Bayes' theorem:

   # From HPV to Normal across age groups:
   # Note: because all other states except Death are assumed to be dependent and disjoint 
   # for regression to states normal or Infection from states LSIL & HSIL, one may calculate 
   # all other relevant states from the complement of the transitions:
    HPV_Well_15to24 ~ dbeta(alpha.HPVtoNormal_15to24, beta.HPVtoNormal_15to24)
    HPV_Well_25to29 ~ dbeta(alpha.HPVtoNormal_25to29, beta.HPVtoNormal_25to29)
    HPV_Well_30toEnd ~ dbeta(alpha.HPVtoNormal_30toPlus, beta.HPVtoNormal_30toPlus)

### END OF SUB-MODEL 4.

### SUB-MODEL 5: LSIL & HSIL PROGRESSION:
# Note: this is equivalent to standard PSA, as it is technically sampling directly
# from an integrated probability distribution and not explicitly propogated into a 
# posterior via Bayes' theorem:

   # From LSIL & HSIL to Normal or Infection across age groups:
   # Note: because all other states except Death are assumed to be dependent and disjoint 
   # for regression to states normal or Infection from states LSIL & HSIL, one may calculate 
   # all other relevant states from the complement of the transitions:
    LSIL_15_34 ~ dbeta(alpha.LSIL_15to34, beta.LSIL_15to34)
    LSIL_35_85 ~ dbeta(alpha.LSIL_35to85, beta.LSIL_35to85)
    HSIL_n ~ dbeta(alpha.HSIL, beta.HSIL)

### END OF SUB-MODEL 5.

 }
"
writeLines(text = model_String, con = "data/jags_model.txt")

# Transform data into list so that it can be parsed by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
  Nstud.vac = Nstud.vac, 
  rA.vac = rA.vac, nA.vac = nA.vac, 
  rB.vac = rB.vac, nB.vac = nB.vac,
  
  # Population prevalence:
  rate_Incidence = rate_Incidence,
  
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
  "OR.vac", "pEfficacy.vac",
  # Well to infection prevalence:
  "p.age",
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
n.iter <- 20000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "data/jags_model.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS
# Note: mixed predictive checks show that variation in data is largely 
# consistent between studies and that there is only weak evidence that 
# each study is inconsistent with the results of the other studies.

# Attach JAGS model to local envir.:
attach.jags(mod_JAGS)
# It is possible automate convergence. However, the improvement is negligible, 
# and it increased the run time from < 1.25 min to > 2.5 mins. Uncomment 
# the line of code below to add auto chain convergence: 
# mod_JAGS <- autojags(mod_JAGS, Rhat = 1.01)

# Visual posterior inspection: uncomment to run visual inspection
# posterior <- mod_JAGS$BUGSoutput$sims.array
# mcmc_dens_chains(x = posterior, pars = c("pEfficacy.vac"))
## Ensuring p.age chains have converged to similar distribution:
# color_scheme_set("mix-gray-brightblue")
# Trace and density:
# mcmc_combo(x = posterior, pars = "p.age[18]")
# Trace and density:
# mcmc_combo(x = posterior, pars = "p.age[21]")
# Trace and density:
# mcmc_combo(x = posterior, pars = "p.age[25]")
# Trace and density:
# mcmc_combo(x = posterior, pars = "p.age[50]")
# Trace and density:
# mcmc_combo(x = posterior, pars = "p.age[80]")

# ==========================================================================================
# Probability Matrix --------------------------------------------------
# ==========================================================================================
 # This section creates a transition array that is used to simulate the Markov cohort
# for treatments, over time. Note that all transitions are conditional on surviving.

# The following initialises the structure of the transitions between health states:
n_age_min <- 0 # age at baseline
n_age_max <- 85 # maximum age of follow up
n_cycles <- n_age_max - n_age_min # time horizon, number of cycles

# Define the names of the health states of the model:
v_names_states <- c("Well", "Infection", "LSIL", "HSIL", "Stage-I Cancer", 
                    "Stage-II Cancer", "Stage-III Cancer", "Stage-IV Cancer",
                    "Detected.Stage-I Year 1", "Detected.Stage-I Year 2",
                    "Detected.Stage-I Year 3", "Detected.Stage-I Year 4",
                    "Detected.Stage-I Year 5", "Detected.Stage-II Year 1",
                    "Detected.Stage-II Year 2", "Detected.Stage-II Year 3",
                    "Detected.Stage-II Year 4", "Detected.Stage-II Year 5",
                    "Detected.Stage-III Year 1", "Detected.Stage-III Year 2",
                    "Detected.Stage-III Year 3", "Detected.Stage-III Year 4",
                    "Detected.Stage-III Year 5", "Detected.Stage-IV Year 1", 
                    "Detected.Stage-IV Year 2", "Detected.Stage-IV Year 3",
                    "Detected.Stage-IV Year 4", "Detected.Stage-IV Year 5", 
                    "Cancer Survivor", "Death")
n_states <- length(v_names_states) # record the number of health states

# Transition array for each comparative intervention:
a_P_SoC <- array(0, 
               dim = c(n_states, n_states, n_cycles, n.sims),
               dimnames = list(
                v_names_states, v_names_states, 0:(n_cycles - 1), 1:n.sims)
               )
a_P_NT <- a_P_SoC
# Note: due to the natural structure of the transition array, there is no need to fill 
# in probabilities = 0, such as the transitions from Infection to LSIL for ages ≤ 14.

# ==========================================================================================
# Standard of Care --------------------------------------------------------
# ==========================================================================================
# Transitions from Well State ---------------------------------------------
# The following enters all transition probabilities for ages 0-85 for each transition from
# the state Well, across the appropriate time horizon i and all probabilistic 
# simulations j. All transitions, across all states, are conditional on surviving.
for (i in 1:n_cycles) {
 for (j in 1:n.sims) {
  a_P_SoC["Well", "Death", i, j] <-  v_p_mort_lessHPV[i] 
  
  a_P_SoC["Well", "Infection", i, j] <- p.age[j, i] * (1 - v_p_mort_lessHPV[i])
  
  a_P_SoC["Well", "Well", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - p.age[j, i])
 }
}

# Transitions from Infection State ----------------------------------------
# The following enters all transition probabilities for ages 15-24 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 15:(n_cycles - 61)) {
    for (j in 1:n.sims) {
     a_P_SoC["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Infection", "Well", i, j] <- HPV_Well_15to24[j] * (1 - v_p_mort_lessHPV[i])
      
     a_P_SoC["Infection", "LSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_15to24[j]) * (1 - v_p_mort_lessHPV[i])
                                           
     a_P_SoC["Infection", "HSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_15to24[j]) * (1 - v_p_mort_lessHPV[i])
      
     a_P_SoC["Infection", "Infection", i, j] <- (1 - (1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_15to24[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# The following enters all transition probabilities for ages 25-29 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 25:(n_cycles - 56)) {
    for (j in 1:n.sims) {
     a_P_SoC["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Infection", "Well", i, j] <- HPV_Well_25to29[j] * (1 - v_p_mort_lessHPV[i])
      
     a_P_SoC["Infection", "LSIL", i, j] <- ((1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_25to29[j]) * (1 - v_p_mort_lessHPV[i])
                                           
     a_P_SoC["Infection", "HSIL", i, j] <- ((1 - exp(-.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_25to29[j]) * (1 - v_p_mort_lessHPV[i])
      
     a_P_SoC["Infection", "Infection", i, j] <- (1 - (1 - exp(-.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_25to29[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# The following enters all transition probabilities for ages 30-85 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 30:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Infection", "Well", i, j] <- HPV_Well_30toEnd[j] * (1 - v_p_mort_lessHPV[i])
      
     a_P_SoC["Infection", "LSIL", i, j] <- ((1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_30toEnd[j]) * (1 - v_p_mort_lessHPV[i])
                                           
     a_P_SoC["Infection", "HSIL", i, j] <- ((1 - exp(-.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_30toEnd[j]) * (1 - v_p_mort_lessHPV[i])
      
     a_P_SoC["Infection", "Infection", i, j] <- (1 - (1 - exp(-.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_30toEnd[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from LSIL State ---------------------------------------------
# The following enters all transition probabilities for ages 15-34 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 15:(n_cycles - 51)) {
    for (j in 1:n.sims) {
     a_P_SoC["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
          
     a_P_SoC["LSIL", "Well", i, j] <- (LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i]) * 0.9 # split proportion between LSIL to Well
     
     a_P_SoC["LSIL", "Infection", i, j] <- (LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i]) * 0.1 # split proportion between LSIL to Infection 
     
     a_P_SoC["LSIL", "HSIL", i, j] <- (1 - exp(-0.1 * 6)) * (1 - LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["LSIL", "LSIL", i, j] <- (1 - (1 - exp(-0.1 * 6))) * (1 - LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# The following enters all transition probabilities for ages ≥25 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 35:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
          
     a_P_SoC["LSIL", "Well", i, j] <- (LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i]) * 0.9 # split proportion between LSIL to Well
     
     a_P_SoC["LSIL", "Infection", i, j] <- (LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i]) * 0.1 # split proportion between LSIL to Infection 
     
     a_P_SoC["LSIL", "HSIL", i, j] <- (1 - exp(-0.1 * 6)) * (1 - LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["LSIL", "LSIL", i, j] <- (1 - (1 - exp(-0.1 * 6))) * (1 - LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from HSIL State ---------------------------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state HSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["HSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["HSIL", "Well", i, j] <- HSIL_n[j] * (1 - v_p_mort_lessHPV[i]) * 0.5 # Proportion reverting to Well
     
     a_P_SoC["HSIL", "LSIL", i, j] <- HSIL_n[j] * (1 - v_p_mort_lessHPV[i]) * 0.5 # Proportion reverting to LSIL
     
     a_P_SoC["HSIL", "Stage-I Cancer", i, j] <- (1 - exp(-0.4 * 10)) * (1 - v_p_mort_lessHPV[i]) * (1 - HSIL_n[j]) # 1 - e^(-0.4 * 10) is progression rate from HSIL to cancer
     
     a_P_SoC["HSIL", "HSIL", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - HSIL_n[j]) * (1 - (1 - exp(-0.4 * 10)))
    }
}

# Transitions from Stage-I Cervix Cancer State ----------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-I Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Stage-I Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Stage-I Cancer", "Detected.Stage-I Year 1", i, j] <- StageI.Detect.mu * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Stage-I Cancer", "Stage-II Cancer", i, j] <-  Stage.I.canc[j] * (1 - v_p_mort_lessHPV[i]) * (1 - StageI.Detect.mu)
     
     a_P_SoC["Stage-I Cancer", "Stage-I Cancer", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - Stage.I.canc[j]) * (1 - StageI.Detect.mu)
    }
}

# Transitions from Stage-II Cervix Cancer State ---------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-II Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Stage-II Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Stage-II Cancer", "Detected.Stage-II Year 1", i, j] <- StageII.Detect.mu * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Stage-II Cancer", "Stage-III Cancer", i, j] <- Stage.II.canc[j] * (1 - v_p_mort_lessHPV[i]) * (1 - StageII.Detect.mu)
     
     a_P_SoC["Stage-II Cancer", "Stage-II Cancer", i, j] <- (1 - Stage.II.canc[j]) * (1 - v_p_mort_lessHPV[i]) * (1 - StageII.Detect.mu)
    }
}

# Transitions from Stage-III Cervix Cancer State --------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-III Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Stage-III Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Stage-III Cancer", "Detected.Stage-III Year 1", i, j] <- StageIII.Detect.mu * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Stage-III Cancer", "Stage-IV Cancer", i, j] <- Stage.III.canc[j] * (1 - StageIII.Detect.mu) * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Stage-III Cancer", "Stage-III Cancer", i, j] <- (1 - Stage.III.canc[j]) * (1 - StageIII.Detect.mu) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from Stage-IV Cervix Cancer State ---------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-IV Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Stage-IV Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_SoC["Stage-IV Cancer", "Detected.Stage-IV Year 1", i, j] <- StageIV.Detected[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Stage-IV Cancer", "Stage-IV Cancer", i, j] <- (1 - StageIV.Detected[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from All Stage-I Cancer Survivor/Detected States ---------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-I Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Detected.Stage-I Year 1", "Detected.Stage-I Year 2", i, j] <- surv.StageI_year1[j] * (1 - v_p_mort_lessHPV[i]) # conditional on surviving all cause mortality
     
     a_P_SoC["Detected.Stage-I Year 1", "Death", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - surv.StageI_year1[j]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-I Year 2", "Detected.Stage-I Year 3", i, j] <- surv.StageI_year2[j] * (1 - v_p_mort_lessHPV[i]) 
     
     a_P_SoC["Detected.Stage-I Year 2", "Death", i, j] <- (1 - surv.StageI_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-I Year 3", "Detected.Stage-I Year 4", i, j] <- surv.StageI_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-I Year 3", "Death", i, j] <- (1 - surv.StageI_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_SoC["Detected.Stage-I Year 4", "Detected.Stage-I Year 5", i, j] <- surv.StageI_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-I Year 4", "Death", i, j] <- (1 - surv.StageI_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-I Year 5", "Cancer Survivor", i, j] <- surv.StageI_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-I Year 5", "Death", i, j] <- (1 - surv.StageI_year5[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     }
}

# Transitions from All Stage-II Cancer Survivor/Detected States -----------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-II Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Detected.Stage-II Year 1", "Detected.Stage-II Year 2", i, j] <- surv.StageII_year1[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-II Year 1", "Death", i, j] <- (1 - surv.StageII_year1[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-II Year 2", "Detected.Stage-II Year 3", i, j] <- surv.StageII_year2[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-II Year 2", "Death", i, j] <- (1 - surv.StageII_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-II Year 3", "Detected.Stage-II Year 4", i, j] <- surv.StageII_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-II Year 3", "Death", i, j] <- (1 - surv.StageII_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_SoC["Detected.Stage-II Year 4", "Detected.Stage-II Year 5", i, j] <- surv.StageII_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-II Year 4", "Death", i, j] <- (1 - surv.StageII_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-II Year 5", "Cancer Survivor", i, j] <- surv.StageII_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-II Year 5", "Death", i, j] <- (1 - surv.StageII_year5[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
    }
}

# Transitions from All Stage-III Cancer Survivor/Detected States ----------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-III Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Detected.Stage-III Year 1", "Detected.Stage-III Year 2", i, j] <- surv.StageIII_year1[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-III Year 1", "Death", i, j] <- (1 - surv.StageIII_year1[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-III Year 2", "Detected.Stage-III Year 3", i, j] <- surv.StageIII_year2[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-III Year 2", "Death", i, j] <- (1 - surv.StageIII_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-III Year 3", "Detected.Stage-III Year 4", i, j] <- surv.StageIII_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-III Year 3", "Death", i, j] <- (1 - surv.StageIII_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_SoC["Detected.Stage-III Year 4", "Detected.Stage-III Year 5", i, j] <- surv.StageIII_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-III Year 4", "Death", i, j] <- (1 - surv.StageIII_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-III Year 5", "Cancer Survivor", i, j] <- surv.StageIII_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-III Year 5", "Death", i, j] <- (1 - surv.StageIII_year5[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
    }
}

# Transitions from All Stage-IV Cancer Survivor/Detected States -----------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-IV Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Detected.Stage-IV Year 1", "Detected.Stage-IV Year 2", i, j] <- surv.StageIV_year1[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-IV Year 1", "Death", i, j] <- (1 - surv.StageIV_year1[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-IV Year 2", "Detected.Stage-IV Year 3", i, j] <- surv.StageIV_year2[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-IV Year 2", "Death", i, j] <- (1 - surv.StageIV_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-IV Year 3", "Detected.Stage-IV Year 4", i, j] <- surv.StageIV_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-IV Year 3", "Death", i, j] <- (1 - surv.StageIV_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_SoC["Detected.Stage-IV Year 4", "Detected.Stage-IV Year 5", i, j] <- surv.StageIV_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-IV Year 4", "Death", i, j] <- (1 - surv.StageIV_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_SoC["Detected.Stage-IV Year 5", "Cancer Survivor", i, j] <- surv.StageIV_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Detected.Stage-IV Year 5", "Death", i, j] <- (1 - surv.StageIV_year5[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from Cancer Survivor State ----------------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state of Cancer Survivor, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_SoC["Cancer Survivor", "Cancer Survivor", i, j] <- (1 - v_p_mort_lessHPV[i])
     
     a_P_SoC["Cancer Survivor", "Death", i, j] <- v_p_mort_lessHPV[i]
    }
}

# Transitions from Death State --------------------------------------------
# The following enters all transition probabilities for ages 0-85 for each transition from
# the state of Death, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 1:n_cycles) {
 for (j in 1:n.sims) {
  a_P_SoC["Death", "Death", i, j] <- 1
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
# simulations j. All transitions are conditional on surviving.
for (i in 1:n_cycles) {
 for (j in 1:n.sims) {
  a_P_NT["Well", "Death", i, j] <-  v_p_mort_lessHPV[i] 
  
  a_P_NT["Well", "Infection", i, j] <-  (1 - (1 - p.age[j, i]) ^ (1 - pEfficacy.vac[j])) * (1 - v_p_mort_lessHPV[i])
  
  a_P_NT["Well", "Well", i, j] <-  ((1 - p.age[j, i]) ^ (1 - pEfficacy.vac[j])) * (1 - v_p_mort_lessHPV[i])
 }
}

# Transitions from Infection State ----------------------------------------
# The following enters all transition probabilities for ages 15-24 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 15:(n_cycles - 61)) {
    for (j in 1:n.sims) {
     a_P_NT["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Infection", "Well", i, j] <- HPV_Well_15to24[j] * (1 - v_p_mort_lessHPV[i])
      
     a_P_NT["Infection", "LSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_15to24[j]) * (1 - v_p_mort_lessHPV[i])
                                           
     a_P_NT["Infection", "HSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_15to24[j]) * (1 - v_p_mort_lessHPV[i])
      
     a_P_NT["Infection", "Infection", i, j] <- (1 - (1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_15to24[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# The following enters all transition probabilities for ages 25-29 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 25:(n_cycles - 56)) {
    for (j in 1:n.sims) {
     a_P_NT["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Infection", "Well", i, j] <- HPV_Well_25to29[j] * (1 - v_p_mort_lessHPV[i])
      
     a_P_NT["Infection", "LSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_25to29[j]) * (1 - v_p_mort_lessHPV[i])
                                           
     a_P_NT["Infection", "HSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-.2 * 3)) * 0.9) * (1 - HPV_Well_25to29[j]) * (1 - v_p_mort_lessHPV[i])
      
     a_P_NT["Infection", "Infection", i, j] <- (1 - (1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_25to29[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# The following enters all transition probabilities for ages 30-85 for each transition from
# the state Infection, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 30:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Infection", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Infection", "Well", i, j] <- HPV_Well_30toEnd[j] * (1 - v_p_mort_lessHPV[i])
      
     a_P_NT["Infection", "LSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_30toEnd[j]) * (1 - v_p_mort_lessHPV[i])
                                           
     a_P_NT["Infection", "HSIL", i, j] <- ((1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_30toEnd[j]) * (1 - v_p_mort_lessHPV[i])
      
     a_P_NT["Infection", "Infection", i, j] <- (1 - (1 - exp(-0.2 * 3)) * 0.1) * (1 - (1 - exp(-0.2 * 3)) * 0.9) * (1 - HPV_Well_30toEnd[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from LSIL State ---------------------------------------------
# The following enters all transition probabilities for ages 15-34 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 15:(n_cycles - 51)) {
    for (j in 1:n.sims) {
     a_P_NT["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
          
     a_P_NT["LSIL", "Well", i, j] <- (LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i]) * 0.9 # split proportion between LSIL to Well
     
     a_P_NT["LSIL", "Infection", i, j] <- (LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i]) * 0.1 # split proportion between LSIL to Infection 
     
     a_P_NT["LSIL", "HSIL", i, j] <- (1 - exp(-0.1 * 6)) * (1 - LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["LSIL", "LSIL", i, j] <- (1 - (1 - exp(-0.1 * 6))) * (1 - LSIL_15_34[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# The following enters all transition probabilities for ages ≥25 for each transition from
# the state LSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 35:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["LSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
          
     a_P_NT["LSIL", "Well", i, j] <- (LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i]) * 0.9 # split proportion between LSIL to Well
     
     a_P_NT["LSIL", "Infection", i, j] <- (LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i]) * 0.1 # split proportion between LSIL to Infection 
     
     a_P_NT["LSIL", "HSIL", i, j] <- (1 - exp(-0.1 * 6)) * (1 - LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["LSIL", "LSIL", i, j] <- (1 - (1 - exp(-0.1 * 6))) * (1 - LSIL_35_85[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from HSIL State ---------------------------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state HSIL, across the appropriate time horizon i and all probabilistic 
# simulations j.
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["HSIL", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["HSIL", "Well", i, j] <- HSIL_n[j] * (1 - v_p_mort_lessHPV[i]) * 0.5 # Proportion reverting to Well
     
     a_P_NT["HSIL", "LSIL", i, j] <- HSIL_n[j] * (1 - v_p_mort_lessHPV[i]) * 0.5 # Proportion reverting to LSIL
     
     a_P_NT["HSIL", "Stage-I Cancer", i, j] <- (1 - exp(-0.4 * 10)) * (1 - v_p_mort_lessHPV[i]) * (1 - HSIL_n[j]) # 1 - e^(-0.4 * 10) is progression rate from HSIL to cancer
     
     a_P_NT["HSIL", "HSIL", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - HSIL_n[j]) * (1 - (1 - exp(-0.4 * 10)))
    }
}

# Transitions from Stage-I Cervix Cancer State ----------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-I Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Stage-I Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Stage-I Cancer", "Detected.Stage-I Year 1", i, j] <- StageI.Detect.mu * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Stage-I Cancer", "Stage-II Cancer", i, j] <-  Stage.I.canc[j] * (1 - v_p_mort_lessHPV[i]) * (1 - StageI.Detect.mu)
     
     a_P_NT["Stage-I Cancer", "Stage-I Cancer", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - Stage.I.canc[j]) * (1 - StageI.Detect.mu)
    }
}

# Transitions from Stage-II Cervix Cancer State ---------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-II Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Stage-II Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Stage-II Cancer", "Detected.Stage-II Year 1", i, j] <- StageII.Detect.mu * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Stage-II Cancer", "Stage-III Cancer", i, j] <- Stage.II.canc[j] * (1 - v_p_mort_lessHPV[i]) * (1 - StageII.Detect.mu)
     
     a_P_NT["Stage-II Cancer", "Stage-II Cancer", i, j] <- (1 - Stage.II.canc[j]) * (1 - v_p_mort_lessHPV[i]) * (1 - StageII.Detect.mu)
    }
}

# Transitions from Stage-III Cervix Cancer State --------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-III Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Stage-III Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Stage-III Cancer", "Detected.Stage-III Year 1", i, j] <- StageIII.Detect.mu * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Stage-III Cancer", "Stage-IV Cancer", i, j] <- Stage.III.canc[j] * (1 - StageIII.Detect.mu) * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Stage-III Cancer", "Stage-III Cancer", i, j] <- (1 - Stage.III.canc[j]) * (1 - StageIII.Detect.mu) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from Stage-IV Cervix Cancer State ---------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state Stage-IV Cancer, across the appropriate time horizon i and all probabilistic 
# simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Stage-IV Cancer", "Death", i, j] <- v_p_mort_lessHPV[i]
     
     a_P_NT["Stage-IV Cancer", "Detected.Stage-IV Year 1", i, j] <- StageIV.Detected[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Stage-IV Cancer", "Stage-IV Cancer", i, j] <- (1 - StageIV.Detected[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from All Stage-I Cancer Survivor/Detected States ---------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-I Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Detected.Stage-I Year 1", "Detected.Stage-I Year 2", i, j] <- surv.StageI_year1[j] * (1 - v_p_mort_lessHPV[i]) # conditional on surviving all cause mortality
     
     a_P_NT["Detected.Stage-I Year 1", "Death", i, j] <- (1 - v_p_mort_lessHPV[i]) * (1 - surv.StageI_year1[j]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-I Year 2", "Detected.Stage-I Year 3", i, j] <- surv.StageI_year2[j] * (1 - v_p_mort_lessHPV[i]) 
     
     a_P_NT["Detected.Stage-I Year 2", "Death", i, j] <- (1 - surv.StageI_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-I Year 3", "Detected.Stage-I Year 4", i, j] <- surv.StageI_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-I Year 3", "Death", i, j] <- (1 - surv.StageI_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_NT["Detected.Stage-I Year 4", "Detected.Stage-I Year 5", i, j] <- surv.StageI_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-I Year 4", "Death", i, j] <- (1 - surv.StageI_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-I Year 5", "Cancer Survivor", i, j] <- surv.StageI_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-I Year 5", "Death", i, j] <- (1 - surv.StageI_year5[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     }
}

# Transitions from All Stage-II Cancer Survivor/Detected States -----------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-II Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Detected.Stage-II Year 1", "Detected.Stage-II Year 2", i, j] <- surv.StageII_year1[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-II Year 1", "Death", i, j] <- (1 - surv.StageII_year1[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-II Year 2", "Detected.Stage-II Year 3", i, j] <- surv.StageII_year2[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-II Year 2", "Death", i, j] <- (1 - surv.StageII_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-II Year 3", "Detected.Stage-II Year 4", i, j] <- surv.StageII_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-II Year 3", "Death", i, j] <- (1 - surv.StageII_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_NT["Detected.Stage-II Year 4", "Detected.Stage-II Year 5", i, j] <- surv.StageII_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-II Year 4", "Death", i, j] <- (1 - surv.StageII_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-II Year 5", "Cancer Survivor", i, j] <- surv.StageII_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-II Year 5", "Death", i, j] <- (1 - surv.StageII_year5[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
    }
}

# Transitions from All Stage-III Cancer Survivor/Detected States ----------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-III Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Detected.Stage-III Year 1", "Detected.Stage-III Year 2", i, j] <- surv.StageIII_year1[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-III Year 1", "Death", i, j] <- (1 - surv.StageIII_year1[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-III Year 2", "Detected.Stage-III Year 3", i, j] <- surv.StageIII_year2[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-III Year 2", "Death", i, j] <- (1 - surv.StageIII_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-III Year 3", "Detected.Stage-III Year 4", i, j] <- surv.StageIII_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-III Year 3", "Death", i, j] <- (1 - surv.StageIII_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_NT["Detected.Stage-III Year 4", "Detected.Stage-III Year 5", i, j] <- surv.StageIII_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-III Year 4", "Death", i, j] <- (1 - surv.StageIII_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-III Year 5", "Cancer Survivor", i, j] <- surv.StageIII_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-III Year 5", "Death", i, j] <- (1 - surv.StageIII_year5[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
    }
}

# Transitions from All Stage-IV Cancer Survivor/Detected States -----------
# The following enters all transition probabilities for ages 15-85 for each transition from
# all states of Detected Stage-IV Cancer, across the appropriate time horizon i and all 
# probabilistic simulations j. 
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Detected.Stage-IV Year 1", "Detected.Stage-IV Year 2", i, j] <- surv.StageIV_year1[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-IV Year 1", "Death", i, j] <- (1 - surv.StageIV_year1[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-IV Year 2", "Detected.Stage-IV Year 3", i, j] <- surv.StageIV_year2[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-IV Year 2", "Death", i, j] <- (1 - surv.StageIV_year2[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-IV Year 3", "Detected.Stage-IV Year 4", i, j] <- surv.StageIV_year3[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-IV Year 3", "Death", i, j] <- (1 - surv.StageIV_year3[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]

     a_P_NT["Detected.Stage-IV Year 4", "Detected.Stage-IV Year 5", i, j] <- surv.StageIV_year4[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-IV Year 4", "Death", i, j] <- (1 - surv.StageIV_year4[j]) * (1 - v_p_mort_lessHPV[i]) + v_p_mort_lessHPV[i]
     
     a_P_NT["Detected.Stage-IV Year 5", "Cancer Survivor", i, j] <- surv.StageIV_year5[j] * (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Detected.Stage-IV Year 5", "Death", i, j] <- (1 - surv.StageIV_year5[j]) * (1 - v_p_mort_lessHPV[i])
    }
}

# Transitions from Cancer Survivor State ----------------------------------
# The following enters all transition probabilities for ages 15-85 for each transition from
# the state of Cancer Survivor, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 15:n_cycles) {
    for (j in 1:n.sims) {
     a_P_NT["Cancer Survivor", "Cancer Survivor", i, j] <- (1 - v_p_mort_lessHPV[i])
     
     a_P_NT["Cancer Survivor", "Death", i, j] <- v_p_mort_lessHPV[i]
    }
}

# Transitions from Death State --------------------------------------------
# The following enters all transition probabilities for ages 0-85 for each transition from
# the state of Death, across the appropriate time horizon i and all 
# probabilistic simulations j.
for (i in 1:n_cycles) {
 for (j in 1:n.sims) {
  a_P_NT["Death", "Death", i, j] <- 1
  }
}

# End file ----------------------------------------------------------------