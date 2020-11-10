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
# This is the code that runs the final Cost-Effectiveness Analysis script for the HPV model 
# developed by Sinanovic, E., et al. 2009): "The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa." Note that 
# this code sits on top of the code for the Markov Models for either treatment.

d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs

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
u_Death <- 0 # utility of being in state of Death for one cycle

# Vector of state utilities for both treatments:
v_utilities <- c("Well" = u_Well, "Infection" = u_Infection, "LSIL" = u_LSIL, "HSIL" = u_HSIL,
                 "Stage-I Cancer" = u_StageI, "Stage-II Cancer" = u_StageII,
                 "Stage-III Cancer" = u_StageIII, "Stage-IV Cancer" = u_StageIV,
                 "Detected.Stage-I Year 1" = 0,
                 "Detected.Stage-I Year 2" = 0,
                 "Detected.Stage-I Year 3" = 0,
                 "Detected.Stage-I Year 4" = 0,
                 "Detected.Stage-I Year 5" = 0,
                 "Detected.Stage-II Year 1" = 0,
                 "Detected.Stage-II Year 2" = 0,
                 "Detected.Stage-II Year 3" = 0,
                 "Detected.Stage-II Year 4" = 0,
                 "Detected.Stage-II Year 5" = 0,
                 "Detected.Stage-III Year 1" = 0,
                 "Detected.Stage-III Year 2" = 0,
                 "Detected.Stage-III Year 3" = 0,
                 "Detected.Stage-III Year 4" = 0,
                 "Detected.Stage-III Year 5" = 0,
                 "Detected.Stage-IV Year 1" = 0,
                 "Detected.Stage-IV Year 2" = 0,
                 "Detected.Stage-IV Year 3" = 0,
                 "Detected.Stage-IV Year 4" = 0,
                 "Detected.Stage-IV Year 5" = 0,
                 "Cancer Survivor" = u_CancerSurvivor,
                 "Death" = u_Death)

# Health State costs ------------------------------------------------------
# All costs in $US:
c_Vaccine <- 559 # once off cost of vaccine from age 12.

# All costs for LSIL, HSIL, Infection, and Cancer Stages are at screening intervals of 
# 30, 40, and 50:
c_Infection <- 276 # cost of Infection screening
c_LSIL <- 39 # cost of LSIL screening
c_HSIL <- 689 # cost of HSIL screening
c_Cancer <- 93 # cost of cervical cytology screening

# Cost for treating all cancer stages:
c_StageI <- 4615 # cost of treatment of Stage I Cancer for one cycle
c_StageII <- 6307 # cost of treatment of Stage II Cancer for one cycle
c_StageIII <- 6307 # cost of treatment of Stage III Cancer for one cycle
c_StageIV <- 8615 # cost of treatment of Stage IV Cancer for one cycle

# Vector of state costs under new treatment:
v_costs_SQ <- matrix(c("Well" = 0, "Infection" = 0, "LSIL" = 0, "HSIL" = 0,
                 "Stage-I Cancer" = 0, "Stage-II Cancer" = 0,
                 "Stage-III Cancer" = 0, "Stage-IV Cancer" = 0,
                 "Detected.Stage-I Year 1" = c_StageI,
                 "Detected.Stage-I Year 2" = c_StageI,
                 "Detected.Stage-I Year 3" = c_StageI,
                 "Detected.Stage-I Year 4" = c_StageI,
                 "Detected.Stage-I Year 5" = c_StageI,
                 "Detected.Stage-II Year 1" = c_StageII,
                 "Detected.Stage-II Year 2" = c_StageII,
                 "Detected.Stage-II Year 3" = c_StageII,
                 "Detected.Stage-II Year 4" = c_StageII,
                 "Detected.Stage-II Year 5" = c_StageII,
                 "Detected.Stage-III Year 1" = c_StageIII,
                 "Detected.Stage-III Year 2" = c_StageIII,
                 "Detected.Stage-III Year 3" = c_StageIII,
                 "Detected.Stage-III Year 4" = c_StageIII,
                 "Detected.Stage-III Year 5" = c_StageIII,
                 "Detected.Stage-IV Year 1" = c_StageIV,
                 "Detected.Stage-IV Year 2" = c_StageIV,
                 "Detected.Stage-IV Year 3" = c_StageIV,
                 "Detected.Stage-IV Year 4" = c_StageIV,
                 "Detected.Stage-IV Year 5" = c_StageIV,
                 "Cancer Survivor" = 0,
                 "Death" = 0), nrow = n_states, ncol = n_t, dimnames = list(v_n))

# Costs for HPV Infection screening for Ages 30, 40, and 50.
v_costs_SQ["Infection", 30] <- 309
v_costs_SQ["Infection", 40] <- 309
v_costs_SQ["Infection", 50] <- 309
# Costs for LSIL screening for Ages 30, 40, and 50.
v_costs_SQ["Infection", 30] <- 61
v_costs_SQ["Infection", 40] <- 61
v_costs_SQ["Infection", 50] <- 61
# Costs for HSIL screening for Ages 30, 40, and 50.
v_costs_SQ["Infection", 30] <- 764
v_costs_SQ["Infection", 40] <- 764
v_costs_SQ["Infection", 50] <- 764
# Costs for cervical cytology screening for Ages 30, 40, and 50.
v_costs_SQ["Stage-I Cancer", 30] <- 93
v_costs_SQ["Stage-I Cancer", 40] <- 93
v_costs_SQ["Stage-I Cancer", 50] <- 93
v_costs_SQ["Stage-II Cancer", 30] <- 93
v_costs_SQ["Stage-II Cancer", 40] <- 93
v_costs_SQ["Stage-II Cancer", 50] <- 93
v_costs_SQ["Stage-III Cancer", 30] <- 93
v_costs_SQ["Stage-III Cancer", 40] <- 93
v_costs_SQ["Stage-III Cancer", 50] <- 93
v_costs_SQ["Stage-IV Cancer", 30] <- 93
v_costs_SQ["Stage-IV Cancer", 40] <- 93
v_costs_SQ["Stage-IV Cancer", 50] <- 93

# Vector of state costs under new treatment:
v_costs_Vac <- matrix(c("Well" = 0, "Infection" = 0, "LSIL" = 0, "HSIL" = 0,
                 "Stage-I Cancer" = 0, "Stage-II Cancer" = 0,
                 "Stage-III Cancer" = 0, "Stage-IV Cancer" = 0,
                 "Detected.Stage-I Year 1" = c_StageI,
                 "Detected.Stage-I Year 2" = c_StageI,
                 "Detected.Stage-I Year 3" = c_StageI,
                 "Detected.Stage-I Year 4" = c_StageI,
                 "Detected.Stage-I Year 5" = c_StageI,
                 "Detected.Stage-II Year 1" = c_StageII,
                 "Detected.Stage-II Year 2" = c_StageII,
                 "Detected.Stage-II Year 3" = c_StageII,
                 "Detected.Stage-II Year 4" = c_StageII,
                 "Detected.Stage-II Year 5" = c_StageII,
                 "Detected.Stage-III Year 1" = c_StageIII,
                 "Detected.Stage-III Year 2" = c_StageIII,
                 "Detected.Stage-III Year 3" = c_StageIII,
                 "Detected.Stage-III Year 4" = c_StageIII,
                 "Detected.Stage-III Year 5" = c_StageIII,
                 "Detected.Stage-IV Year 1" = c_StageIV,
                 "Detected.Stage-IV Year 2" = c_StageIV,
                 "Detected.Stage-IV Year 3" = c_StageIV,
                 "Detected.Stage-IV Year 4" = c_StageIV,
                 "Detected.Stage-IV Year 5" = c_StageIV,
                 "Cancer Survivor" = 0,
                 "Death" = 0), nrow = n_states, ncol = n_t, dimnames = list(v_n))

v_costs_Vac["Well", 12] <- 570
# Costs for HPV Infection screening for Ages 30, 40, and 50.
v_costs_Vac["Infection", 30] <- 309
v_costs_Vac["Infection", 40] <- 309
v_costs_Vac["Infection", 50] <- 309
# Costs for LSIL screening for Ages 30, 40, and 50.
v_costs_Vac["Infection", 30] <- 61
v_costs_Vac["Infection", 40] <- 61
v_costs_Vac["Infection", 50] <- 61
# Costs for HSIL screening for Ages 30, 40, and 50.
v_costs_Vac["Infection", 30] <- 764
v_costs_Vac["Infection", 40] <- 764
v_costs_Vac["Infection", 50] <- 764
# Costs for cervical cytology screening for Ages 30, 40, and 50.
v_costs_Vac["Stage-I Cancer", 30] <- 93
v_costs_Vac["Stage-I Cancer", 40] <- 93
v_costs_Vac["Stage-I Cancer", 50] <- 93
v_costs_Vac["Stage-II Cancer", 30] <- 93
v_costs_Vac["Stage-II Cancer", 40] <- 93
v_costs_Vac["Stage-II Cancer", 50] <- 93
v_costs_Vac["Stage-III Cancer", 30] <- 93
v_costs_Vac["Stage-III Cancer", 40] <- 93
v_costs_Vac["Stage-III Cancer", 50] <- 93
v_costs_Vac["Stage-IV Cancer", 30] <- 93
v_costs_Vac["Stage-IV Cancer", 40] <- 93
v_costs_Vac["Stage-IV Cancer", 50] <- 93


m_u_SQ <- m_u_SQ <- matrix(NA, n.sims, n_t)
for (i in 1:n.sims) {
 for (j in 1:n_t) {
  m_u_SQ[i, j] <- v_utilities %*% m_M_ad_1[j, , i]
 }
}
m_c_SQ <- matrix(NA, n.sims, n_t)
for (i in 1:n.sims) {
 for (j in 1:n_t) {
  m_c_SQ[i, j] <- v_costs_SQ[, j] %*% m_M_ad_1[j, , i]
 }
}
m_u_NT <- matrix(NA, n.sims, n_t)
for (i in 1:n.sims) {
 for (j in 1:n_t) {
  m_u_NT[i, j] <- v_utilities %*% m_M_ad_2[j, , i]
 }
}
m_c_NT <- matrix(NA, n.sims, n_t)
for (i in 1:n.sims) {
 for (j in 1:n_t) {
  m_c_NT[i, j] <- v_costs_Vac[, j] %*% m_M_ad_2[j, , i]
 }
}



c <- matrix(NA, n.sims, 2)
c[, 1] <- apply(m_c_SQ, 1, sum)
c[, 2] <- apply(m_c_NT, 1, sum)

e <- matrix(NA, n.sims, 2)
e[, 1] <- apply(m_u_SQ, 1, sum)
e[, 2] <- apply(m_u_NT, 1, sum)

v_names_str <- c("Status quo", "Bivalent HPV Vaccine") # strategy names
m <- bcea(e, c, ref = 2, interventions = v_names_str)
BCEA::summary.bcea(m)
BCEA::plot.bcea(m)
BCEA::contour.bcea(m)

# Model run time ----------------------------------------------------------
end_time <- Sys.time()
end_time - start_time