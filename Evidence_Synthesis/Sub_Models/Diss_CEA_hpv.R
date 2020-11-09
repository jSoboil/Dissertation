library(tidyverse)
library(dampack)
library(readxl)
library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(reshape2)
library(BCEA)

source("HPV_Markov_Model.R")

start_time <- Sys.time()

# ==========================================================================================
# Cost-Effectiveness Analysis ---------------------------------------------
# ==========================================================================================
# This is code that runs the final Cost-Effectiveness Analysis script for the HPV model developed 
# by Sinanovic, E., et al. 2009): "The potential cost-effectiveness of adding a human papillomavirus
# vaccine to the cervical cancer screening programme in South Africa." Note that this code sits on
# top of the code for the Markov Models for either treatment.

# Note: m_M_ad_1 is the Status Quo treatment model; m_M_ad_2 is the vaccine treatment model.

# Health State utilities --------------------------------------------------
u_Well <- 1 # utility of being in Well for one cycle
u_Infection <- 1 # utility of being in Infected for one cycle
u_LSIL <- 0.91 # utility of being in LSIL for one cycle
u_HSIL <- 0.87 # utility of being in HSIL for one cycle
u_StageI <- 0.65 # utility of being in Stage I Cancer for one cycle
u_StageII <- 0.56 # utility of being in Stage II Cancer for one cycle
u_StageIII <- 0.56 # utility of being in Stage III Cancer for one cycle
u_StageIV <- 0.48 # utility of being in Stage IV Cancer for one cycle
u_CancerSurvivor <- 0.84 # utility of being in Stage II Cancer for one cycle

# Vector of state utilities for both treatments:
v_Utilities <- array(c("Well" = u_Well, "Infection" = u_Infection, "LSIL" = u_LSIL, "HSIL" = u_HSIL,
                 "Stage-I Cancer" = u_StageI, "Stage-II Cancer" = u_StageII,
                 "Stage-III Cancer" = u_StageIII, "Stage-IV Cancer" = u_StageIV,
                 "Detected.Stage-I Year 1" = u_StageI,
                 "Detected.Stage-I Year 2" = u_StageI,
                 "Detected.Stage-I Year 3" = u_StageI,
                 "Detected.Stage-I Year 4" = u_StageI,
                 "Detected.Stage-I Year 5" = u_StageI,
                 "Detected.Stage-II Year 1" = u_StageII,
                 "Detected.Stage-II Year 2" = u_StageII,
                 "Detected.Stage-II Year 3" = u_StageII,
                 "Detected.Stage-II Year 4" = u_StageII,
                 "Detected.Stage-II Year 5" = u_StageII,
                 "Detected.Stage-III Year 1" = u_StageIII,
                 "Detected.Stage-III Year 2" = u_StageIII,
                 "Detected.Stage-III Year 3" = u_StageIII,
                 "Detected.Stage-III Year 4" = u_StageIII,
                 "Detected.Stage-III Year 5" = u_StageIII,
                 "Detected.Stage-IV Year 1" = u_StageIV,
                 "Detected.Stage-IV Year 2" = u_StageIV,
                 "Detected.Stage-IV Year 3" = u_StageIV,
                 "Detected.Stage-IV Year 4" = u_StageIV,
                 "Detected.Stage-IV Year 5" = u_StageIV,
                 "Cancer Survivor" = u_CancerSurvivor,
                 "Death" = 0), dim = c(c(n_t + 1, n_states), n.sims), 
                 dimnames = list(0:n_t, v_n, 1:n.sims))

a_QALYs_SQ <- array(0, dim = c(c(n_t + 1, n_states), n.sims), 
                 dimnames = list(0:n_t, v_n, 1:n.sims))

 for (t in 1:n_t) {
  for (i in 1:n_states) {
   for (j in 1:n.sims) {
    a_QALYs_SQ[t, i, j] <- m_M_ad_1[t, i, j] %*% v_Utilities[t, i, j]
   }
  }
 }

m_M_ad_1
str(a_QALYs_SQ)
str(v_Utilities)
# Health State costs ------------------------------------------------------
# All costs in $US:
c_Vaccine <- 559 # once off cost of vaccine from age 12.

# All costs for LSIL, HSIL, Infection, and Cancer Stages are at screening intervals of 
# 30, 40, and 50:
c_Infection <- 276 # cost of Infection screening
c_LSIL <- 39 # cost of LSIL screening
c_HSIL <- 689 # cost of HSIL screening
c_Cancer <- 60 # cost of cervical cytology screening

# Cost for treating all cancer stages:
c_StageI <- 3923 # cost of treatment of Stage I Cancer for one cycle
c_StageII <- 5361 # cost of treatment of Stage II Cancer for one cycle
c_StageIII <- 5361 # cost of treatment of Stage III Cancer for one cycle
c_StageIV <- 7323 # cost of treatment of Stage IV Cancer for one cycle

# Vector of state costs under status quo treatment:
v_c_SQ <- c("Well" = 0, "Infection" = 0, "LSIL" = u_LSIL, "HSIL" = u_HSIL,
                 "Stage-I Cancer" = u_StageI, "Stage-II Cancer" = u_StageII,
                 "Stage-III Cancer" = u_StageIII, "Stage-IV Cancer" = u_StageIV,
                 "Detected.Stage-I Year 1" = u_StageI,
                 "Detected.Stage-I Year 2" = u_StageI,
                 "Detected.Stage-I Year 3" = u_StageI,
                 "Detected.Stage-I Year 4" = u_StageI,
                 "Detected.Stage-I Year 5" = u_StageI,
                 "Detected.Stage-II Year 1" = u_StageII,
                 "Detected.Stage-II Year 2" = u_StageII,
                 "Detected.Stage-II Year 3" = u_StageII,
                 "Detected.Stage-II Year 4" = u_StageII,
                 "Detected.Stage-II Year 5" = u_StageII,
                 "Detected.Stage-III Year 1" = u_StageIII,
                 "Detected.Stage-III Year 2" = u_StageIII,
                 "Detected.Stage-III Year 3" = u_StageIII,
                 "Detected.Stage-III Year 4" = u_StageIII,
                 "Detected.Stage-III Year 5" = u_StageIII,
                 "Detected.Stage-IV Year 1" = u_StageIV,
                 "Detected.Stage-IV Year 2" = u_StageIV,
                 "Detected.Stage-IV Year 3" = u_StageIV,
                 "Detected.Stage-IV Year 4" = u_StageIV,
                 "Detected.Stage-IV Year 5" = u_StageIV,
                 "Cancer Survivor" = u_CancerSurvivor,
                 "Death" = 0)


str(m_M_ad_1)
