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

start_time <- Sys.time()

source("HPV_Markov_Model.R")

# ==========================================================================================
# Cost-Effectiveness Analysis ---------------------------------------------
# ==========================================================================================
# This is the code that runs the final Cost-Effectiveness Analysis script for the HPV model 
# developed by Sinanovic, E., et al. 2009): "The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa." Note that 
# this code sits on top of the code for the Markov Models for either treatment.

# Note: m_M_ad_1 is the Status Quo treatment model
m_M_ad_1
# m_M_ad_2 is the New Treatment model.
m_M_ad_2

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
v_u_All <- c("Well" = u_Well,
         "Infection" = u_Infection,
         "LSIL" = u_LSIL,
         "HSIL" = u_HSIL,
         "Stage-I Cancer" = u_StageI,
         "Stage-II Cancer" = u_StageII,
         "Stage-III Cancer" = u_StageIII,
         "Stage-IV Cancer" = u_StageIV,
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
         "Death" = u_Death)

# Health State costs ------------------------------------------------------
# All costs in $US:
c_Vaccine <- 570 # once off cost of vaccine from age 12.

# All costs for LSIL, HSIL, Infection, and Cancer Stages are at screening intervals of 
# 30, 40, and 50:
c_Infection <- 309 + 75 # cost of Infection screening
c_LSIL <- 61 # cost of LSIL screening
c_HSIL <- 764 # cost of HSIL screening
c_Cancer <- 93 # cost of cervical cytology screening

# Cost for treating all cancer stages:
c_StageI <- 4615 # cost of treatment of Stage I Cancer for one cycle
c_StageII <- 6307 # cost of treatment of Stage II Cancer for one cycle
c_StageIII <- 6307 # cost of treatment of Stage III Cancer for one cycle
c_StageIV <- 8615 # cost of treatment of Stage IV Cancer for one cycle

# Vector of state costs under Status Quo without screening interval:
v_c_SQ <- c("Well" = 0, 
                        "Infection" = 0, 
                        "LSIL" = 0, 
                        "HSIL" = 0,
                        "Stage-I Cancer" = 0, 
                        "Stage-II Cancer" = 0,
                        "Stage-III Cancer" = 0, 
                        "Stage-IV Cancer" = 0,
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
                        "Death" = 0)

# Vector of state costs under New Treatment without screening interval:
v_c_NT <- c("Well" = 0,
            "Infection" = 0,
            "LSIL" = 0,
            "HSIL" = 0,
            "Stage-I Cancer" = 0,
            "Stage-II Cancer" = 0,
            "Stage-III Cancer" = 0,
            "Stage-IV Cancer" = 0,
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
            "Death" = 0)


# Initialize transition array
a_A <- array(0, 
             dim = c(n_states, n_states, (n_t + 1), n.sims), 
             dimnames = list(v_n, v_n, 0:n_t, 1:n.sims))


# Set first slice of A equal to the initial state vector in its diagonal 
for (i in 1:n.sims) {
 diag(a_A[, , 1, i]) <- v_s_init
}

# Iterative solution to produce the transition array
for (i in 1:n.sims) {
 for (t in 1:n_t){
  a_A[, , t + 1, i] <- m_M_ad_1[t, , i] * a_P_1[, , t, i]
  }
 }

# Array of state and transition utilities under Usual Care
a_R_u <- aperm(array(v_u_All,
                     dim = c(n_states, n_states, n_t + 1, n.sims),
                     dimnames = list(v_n, v_n, 0:n_t, 1:n.sims)), perm = c(2, 1, 3, 4))
# Array of state and transition costs under Usual Care
a_R_c_SQ <- aperm(array(v_c_SQ,
                        dim = c(n_states, n_states, n_t + 1, n.sims),
                        dimnames = list(v_n, v_n, 0:n_t, 1:n.sims)), perm = c(2, 1, 3, 4))

a_R_c_SQ["Infection", "Infection", c(31, 41,51), ] + c_Infection
a_R_c_SQ["Infection", "LSIL", c(31, 41,51), ] + c_Infection
a_R_c_SQ["Infection", "HSIL", c(31, 41,51), ] + c_Infection

a_R_c_SQ["Stage-I Cancer", "Stage-I Cancer", c(31, 41,51), ] + c_Infection
a_R_c_SQ["Stage-II Cancer", "Stage-II Cancer", c(31, 41,51), ] + c_Infection
a_R_c_SQ["Stage-III Cancer", "Stage-III Cancer", c(31, 41,51), ] + c_Infection
a_R_c_SQ["Stage-IV Cancer", "Stage-IV Cancer", c(31, 41,51), ] + c_Infection


for (i in 1:n.sims) {
 x <- (a_A[, , , i] * 0.5) * (a_R_c_SQ["Well", "Infection", c(31, 41,51), i] + c_Infection) 
}






# Create cost matrices:
m_utilities_NT <- m_costs_NT <- m_utilities_SQ <- m_costs_SQ <- matrix(0, n.sims, n_t)
for (i in 1:n.sims) {
 for (j in 0:n_t) {
  ### For Status Quo (screening only):
  ## Costs
  m_costs_SQ[i, j] <- v_c_SQ %*% m_M_ad_1[j, , i]
  ## QALYs
  m_utilities_SQ[i, j] <- v_u_All %*% m_M_ad_1[j, , i]
  ### For New Treatment (screening and vaccine):
  ## Costs
  m_costs_NT[i, j] <- v_c_NT %*% m_M_ad_2[j, , i]
  ## QALYs
  m_utilities_NT[i, j] <- v_u_All %*% m_M_ad_2[j, , i]
 }
}

mean(apply(m_costs_SQ, 1, sum))
mean(apply(m_utilities_SQ, 1, sum))
mean(apply(m_costs_NT, 1, sum))
mean(apply(m_utilities_NT, 1, sum))

# Discount rates:
d_e <- 0.03
d_c <- 0.03
# Discount weight for costs
v_dwc <- 1 / ((1 + d_c) ^ (0:(n_t)))
# Discount weight for effects
v_dwe <- 1 / ((1 + d_e) ^ (0:(n_t)))

# Discounted costs and utilities matrices:
m_utilities_NTdisc <- m_costs_NTdisc <- m_utilities_SQdisc <- m_costs_SQdisc <- matrix(0, n.sims, n_t)

for (i in 1:n.sims) {
 for (j in 0:n_t) {
  ### For Status Quo (screening only):
  ## Costs
  m_costs_SQdisc[i, j] <- m_costs_SQ[i, j] / v_dwc[j]
  ## QALYs
  m_utilities_SQdisc[i, j] <- m_utilities_SQ[i, j] / v_dwe[j]
  ### For New Treatment (screening and vaccine):
  ## Costs
  m_costs_NTdisc[i, j] <- m_costs_NT[i, j] / v_dwc[j]
  ## QALYs
  m_utilities_NTdisc[i, j] <- m_utilities_NT[i, j] / v_dwe[j]
 }
}

# Sum costs and effectiveness across all time points:
Effects <- Costs <- matrix(0, n.sims, 2)
### For Status Quo (screening only):
## Costs
Costs[, 1] <- apply(m_costs_SQdisc, 1, sum)
## QALYs
Effects[, 1] <- apply(m_utilities_SQdisc, 1, sum)
### For New Treatment (screening and vaccine):
## Costs
Costs[, 2] <- apply(m_costs_NTdisc, 1, sum)
## QALYs
Effects[, 2] <- apply(m_utilities_NTdisc, 1, sum)

# Strategy namesL
v_names_str <- c("Status Quo: screening only", "New Treatment: screening & vaccine")

df_cea <- bcea(Effects, Costs, ref = 2, interventions = v_names_str, Kmax = 1500)
summary.bcea(df_cea)
contour2(df_cea, wtp = 500, graph = "ggplot")

## Expected Costs and Utility for both treatments:
E_c <- apply(Costs, 2, mean)
E_u <- apply(Effects, 2, mean)
(E_c[2] - E_c[1]) / (E_u[2] - E_u[1])

cea_summary <- dampack::calculate_icers(cost = E_c, effect = E_u, strategies = v_names_str)
cea_summary

# Model run time ----------------------------------------------------------
end_time <- Sys.time()
end_time - start_time