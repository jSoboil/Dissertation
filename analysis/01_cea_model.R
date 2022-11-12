################################################################################
## BEFORE RUNNING THE MODEL, PLEASE ENSURE THAT YOUR WORKING DIRECTORY IS SET ##
## TO THE DIRECTORY OF THE Rproj FOLDER SAVED ON YOUR COMPUTER. IN RSTUDIO    ##
## THE EASIEST WAY TO DO THIS IS BY PRESSING Ctrl + Shift + H.                ##
################################################################################

# Libraries and Misc ------------------------------------------------------
pkgs <- c("bayesplot", "BCEA", "dampack", "readxl",  "reshape2", "R2jags", 
          "tidyverse", "rjags", "darthtools")
sapply(pkgs, require, character.only = TRUE)

# Detect cores:
options(mc.cores = parallel::detectCores())

# Initialise start time counter:
start_time <- Sys.time()

# Set seed for replication:
set.seed(100)

# Note that the source code for this section sits on top of the code for the Markov Models
# found in the R folder (R/02_markov_model.R):
source("R/02_markov_model.R")

# The advantages of using parallel processing for this model are negligible, which runs 
# consistently < 60 secs on a relatively low-powered Macbook pro 2018 13".

# ==========================================================================================
# Cost-Effectiveness Analysis ---------------------------------------------
# ==========================================================================================
# This is the code that runs the final Cost-Effectiveness Analysis script for the HPV model 
# developed by Sinanovic, E., et al. 2009): "The potential cost-effectiveness of adding a 
# human papillomavirus vaccine to the cervical cancer screening programme in South Africa" 

# m_M_ad_1 is the Status Quo (Screening only) treatment Markov model:
m_M_ad_1
# m_M_ad_2 is the New Treatment (Screening plus Vaccine) Markov model:
m_M_ad_2

# Health State utilities --------------------------------------------------
u_Well <- 1 # utility of being in Well for one cycle
u_Infection <- 1 # utility of being in Infected for one cycle
u_LSIL <- 0.91 # utility of being in LSIL for one cycle
u_HSIL <- 0.87 # utility of being in HSIL for one cycle
u_StageI <- 0.65 # utility of being in Stage I Cancer for one cycle
u_StageII <- 0.56 # utility of being in Stage II Cancer for one cycle
u_StageIII <- 0.56 # utility of being in Stage III Cancer for one cycle
u_StageIV <- 0.48 # utility of being in state Stage IV Cancer for one cycle
u_CancerSurvivor <- 0.84 # utility of being in state Cancer Survivor for one cycle
u_Death <- 0 # utility of being in state of Death for one cycle

# Vector of state utilities for both treatments:
m_u_SQ <- m_u_NT <- matrix(c("Well" = u_Well,
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
                             "Death" = u_Death), 
                           n_states, n_t + 1, dimnames = list(v_n, 0:(n_t)))

# Health costs from a society perspective --------------------------------
## All costs in $US.
c_Vaccine <- 570 # once off cost of vaccine from age 12.
c_Screening <- 93 + 309 + 75 # cost of screening using HPV DNA, VIA, cancer cytology tests.
c_LSIL <- 61 # cost of LSIL treatment*.
c_HSIL <- 764 # cost of HSIL treatment*.
### *Note: HPV 16/18 assumed asymptomatic in original mode and therefore not treated.

## Cost for treatment at each cancer stage:
c_StageI <- 4615 # cost of treatment of Stage I Cancer for one cycle.
c_StageII <- 6307 # cost of treatment of Stage II Cancer for one cycle.
c_StageIII <- 6307 # cost of treatment of Stage III Cancer for one cycle.
c_StageIV <- 8615 # cost of treatment of Stage IV Cancer for one cycle.

# Matrix of state costs based on time interval t under Status Quo treatment:
m_c_SQ <- matrix(c("Well" = 0,
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
            "Death" = 0), 
            n_states, n_t + 1, dimnames = list(v_n, 0:(n_t)))

# Screening costs for ages 30, 40, and 50:
m_c_SQ["Well", c(31, 41, 51)] <-(m_c_SQ["Well", c(31, 41, 51)] + c_Screening)
# LSIL costs for ages 30, 40, and 50:
m_c_SQ["LSIL", c(31, 41, 51)] <-(m_c_SQ["LSIL", c(31, 41, 51)] + c_LSIL)
# HSIL costs for ages 30, 40, and 50:
m_c_SQ["HSIL", c(31, 41, 51)] <-(m_c_SQ["HSIL", c(31, 41, 51)] + c_HSIL)
# Matrix of state costs based on time interval t under New Treatment:
m_c_NT <- m_c_SQ
# Vaccine cost at age 12:
m_c_NT["Well", c(13)] <-(m_c_NT["Well", c(13)] + (c_Vaccine * 0.8)) # Note: assumed 80% vaccine coverage
# Create cost and effects matrices:
m_utilities_NT <- m_costs_NT <- m_utilities_SQ <- m_costs_SQ <- matrix(
 0, n.sims, n_t + 1, dimnames = list(
  1:n.sims, 0:(n_t))
 )

# Loop costs over time interval t and n.sims i:
for (i in 1:n.sims) {
 for (t in 0:n_t) {
  ### For Status-Quo (screening only)
  ## Costs
  m_costs_SQ[i, t]  <- (m_M_ad_1[t, , i] * 0.5) %*% m_c_SQ[, t] # Note: 50% of cohort screened
  ## QALYs
  m_utilities_SQ[i, t]  <- m_M_ad_1[t, , i] %*% m_u_SQ[, t]
  ### For New Treatment (screening and vaccine):
  ## Costs
  m_costs_NT[i, t] <- (m_M_ad_2[t, , i] * 0.5) %*% m_c_NT[, t] # Note: 50% of cohort screened
  ## QALYs
  m_utilities_NT[i, t] <- m_M_ad_2[t, , i] %*% m_u_NT[, t]
 }
}

# Evaluation --------------------------------------------------------------
#### Non-discounted costs and effects:
### For Status Quo (screening only):
## Costs
mean(apply(m_costs_SQ, 1, sum))
## QALYs
mean(apply(m_utilities_SQ, 1, sum))
### For New Treatment (screening and vaccine):
## Costs
mean(apply(m_costs_NT, 1, sum))
## QALYs
mean(apply(m_utilities_NT, 1, sum))

#### Discounted costs and effects
# Discount rates:
d_e <- 0.03
d_c <- 0.03
# Discount weights for costs
v_dwc <- 1 / ((1 + d_c) ^ (0:(n_t)))
# Discount weights for effects
v_dwe <- 1 / ((1 + d_e) ^ (0:(n_t)))

# Discounted costs and utilities matrices:
m_costs_SQdisc <- matrix(0, n.sims, n_t + 1, 
                         dimnames = list(1:n.sims, 0:n_t))
m_utilities_NTdisc <- m_costs_NTdisc <- m_utilities_SQdisc <- m_costs_SQdisc

for (i in 1:n.sims) {
 for (t in 0:n_t) {
  ### For Status Quo (screening only):
  ## Costs
  m_costs_SQdisc[i, t] <- m_costs_SQ[i, t] * v_dwc[t]
  ## QALYs
  m_utilities_SQdisc[i, t] <- m_utilities_SQ[i, t] * v_dwe[t]
  ### For New Treatment (screening and vaccine):
  ## Costs
  m_costs_NTdisc[i, t] <- m_costs_NT[i, t] * v_dwc[t]
  ## QALYs
  m_utilities_NTdisc[i, t] <- m_utilities_NT[i, t] * v_dwe[t]
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

# Strategy names:
v_names_str <- c("Status Quo: screening only", "New Treatment: screening & vaccine")

## BCEA package summary:
df_cea <- bcea(Effects, Costs, ref = 2, interventions = v_names_str,
               Kmax = 5724, plot = TRUE)
BCEA::contour(he = df_cea)
BCEA::ceplane.plot(df_cea, wtp = 5724, graph = "ggplot")
BCEA::ceac.plot(df_cea, graph = "ggplot2")
BCEA::eib.plot(df_cea, graph = "ggplot")
BCEA::ib.plot(df_cea, wtp = 5724, graph = "ggplot")

## Expected Costs and Utility for both treatments across all simulations:
E_c <- apply(Costs, 2, mean)
E_u <- apply(Effects, 2, mean)
# Manual ICER:
(E_c[2] - E_c[1]) / (E_u[2] - E_u[1])
## dampack package summary:
cea_summary <- dampack::calculate_icers(cost = E_c, effect = E_u, strategies = v_names_str)
cea_summary

# Model run time ----------------------------------------------------------
# End time counter:
end_time <- Sys.time()
# Total run time:
end_time - start_time

# End file ----------------------------------------------------------------