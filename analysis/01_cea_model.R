################################################################################
## BEFORE RUNNING THE MODEL, PLEASE ENSURE THAT YOUR WORKING DIRECTORY IS SET ##
## TO THE DIRECTORY OF THE Rproj FOLDER SAVED ON YOUR COMPUTER. IN RSTUDIO    ##
## THE EASIEST WAY TO DO THIS IS BY PRESSING Ctrl + Shift + H.                ##
################################################################################

# Libraries and Misc ------------------------------------------------------
## Names of required packages for analysis:
pkgs <- c("bayesplot", "BCEA", "dampack", "readxl",  "reshape2", "R2jags", 
          "tidyverse", "rjags", "darthtools")

# Install packages not yet installed:
installed_packages <- pkgs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
 install.packages(pkgs[!installed_packages])
}

# Load required packages:
invisible(lapply(pkgs, library, character.only = TRUE))

# Detect and set parallel processing to number of cores:
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

# vector of state names:
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

# Initialize transition-dynamics array under SoC and New Treatment
a_A_SoC <- array(0,
             dim = c(n_states, n_states, (n_t + 1), n.sims),
dimnames = list(v_names_states, v_names_states, 0:n_t)) # Set first slice to the initial state vector in its diagonal diag(a_A_SoC[, , 1]) <- v_m_init
# New Treatment:
a_A_NT <- a_A_SoC

# Set first slice to the initial state vector in its diagonal
for (i in 1:n.sims) {
 diag(a_A_SoC[, , 1, i]) <- v_s_init
 diag(a_A_NT[, , 1, i]) <- v_s_init
}

# Iterative solution to produce the transition-dynamics array under SoC:
for (i in 1:n.sims) {
 for (t in 1:n_t){
  a_A_SoC[, , t + 1, i] <- diag(m_M_SoC[t, , i]) %*% a_P_SoC[ , , t, i] # using Markov trace found in 02_markov_model.R script
  }
}
# Iterative solution to produce the transition-dynamics array under new Treatment:
for (i in 1:n.sims) {
 for (t in 1:n_t){
  a_A_NT[, , t + 1, i] <- diag(m_M_NT[t, , i]) %*% a_P_NT[ , , t, i]
  }
}

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

# vector of costs:
v_c_SoC <- c(0, 0, 0, 0, 0, 0, 0, 0, c_StageI, c_StageI, c_StageI, c_StageI, 
             c_StageI, c_StageII, c_StageII, c_StageII, c_StageII, c_StageII,
             c_StageIII, c_StageIII, c_StageIII, c_StageIII, c_StageIII, 
             c_StageIV, c_StageIV, c_StageIV, c_StageIV, c_StageIV, 0, 0)

# Array of state costs based on time interval t under Status Quo treatment:
a_c_SoC <- array(matrix(v_c_SoC, nrow = n_states, ncol = n_states, byrow = T),
                  dim = c(n_states, n_states, n_t + 1, n.sims),
                  dimnames = list(v_names_states, v_names_states, 0:n_t))

# Screening costs for ages 30, 40, and 50:
a_c_SoC[c("Well", "Infection", "LSIL", "HSIL"), c("Well", "Infection", "LSIL", "HSIL"), c(30, 40, 50), ] <- a_c_SoC[c("Well", "Infection", "LSIL", "HSIL"), c("Well", "Infection", "LSIL", "HSIL"), c(30, 40, 50), ] + c_Screening # everyone not dead or with diagnosed cancer
# LSIL costs for ages 30, 40, and 50:
a_c_SoC["LSIL", c("LSIL", "Well", "Infection", "HSIL"), c(30, 40, 50), ] <- a_c_SoC["LSIL", c("LSIL", "Well", "Infection", "HSIL"), c(30, 40, 50), ] + c_LSIL
# HSIL costs for ages 30, 40, and 50:
a_c_SoC["HSIL", c("HSIL", "Well", "Infection", "LSIL"), c(30, 40, 50), ] <- a_c_SoC["HSIL", c("HSIL", "Well", "Infection", "LSIL"), c(30, 40, 50), ] + c_HSIL

# Array of state costs based on time interval t under New Treatment:
a_c_NT <- a_c_SoC
# Vaccine cost at age 12:
a_c_NT["Well", "Well", 12, ] <- a_c_NT["Well", "Well", 12, ] + c_Vaccine # only applied to Well state since patients can only move to Death from Well till 15

## Total costs Arrays ------------------------------------------------------
# Standard of Care
a_Y_c_SoC <- a_A_SoC * a_c_SoC
a_Y_c_SoC[c("Well", "Infection", "LSIL", "HSIL"), c("Well", "Infection", "LSIL", "HSIL"), c(30, 40, 50), ] <- a_Y_c_SoC[c("Well", "Infection", "LSIL", "HSIL"), c("Well", "Infection", "LSIL", "HSIL"), c(30, 40, 50), ] * 0.5 # 50% of cohort screened

# New Treatment
a_Y_c_NT <- a_A_NT * a_c_NT
a_Y_c_NT[, , 12, ] <- a_Y_c_NT[, , 12, ] * 0.8 # 80% of cohort vaccinated
a_Y_c_NT[c("Well", "Infection", "LSIL", "HSIL"), c("Well", "Infection", "LSIL", "HSIL"), c(30, 40, 50), ] <- a_Y_c_NT[c("Well", "Infection", "LSIL", "HSIL"), c("Well", "Infection", "LSIL", "HSIL"), c(30, 40, 50), ] * 0.5 # 50% of cohort screened

# Health State Utilities --------------------------------------------------
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
v_u_SoC <- c(u_Well, u_Infection, u_LSIL, u_HSIL, u_StageI, u_StageII, u_StageIII,
             u_StageIV, u_StageI, u_StageI, u_StageI, u_StageI, u_StageI, 
             u_StageII, u_StageII, u_StageII, u_StageII, u_StageII, 
             u_StageIII, u_StageIII, u_StageIII, u_StageIII, 
             u_StageIII, u_StageIV, u_StageIV, u_StageIV,
             u_StageIV, u_StageIV, u_CancerSurvivor,
             u_Death)

# Array of state utilities based on time interval t under Status Quo treatment:
a_u_SoC <- array(matrix(v_u_SoC, nrow = n_states, ncol = n_states, byrow = T),
                  dim = c(n_states, n_states, n_t + 1, n.sims),
                  dimnames = list(v_names_states, v_names_states, 0:n_t))

# Array of state utilities based on time interval t under New Treatment:
a_u_NT <- a_u_SoC

## Total QALYs Array -------------------------------------------------------
# Standard of Care
a_Y_u_SoC <- a_A_SoC * a_u_SoC
# New Treatment
a_Y_u_NT <- a_A_NT * a_u_NT

# Vectors of rewards under each Treatment ---------------------------------
m_costs_SoC <- matrix(NA, nrow = n_t + 1, ncol = n.sims, byrow = TRUE) # create matrix with n.sims slices
m_qaly_NT <- m_costs_NT <- m_qaly_SoC <- m_costs_SoC

for (i in 1:n.sims) {
 # Standard of Care:
 m_costs_SoC[, i] <- rowSums(t(colSums(a_Y_c_SoC[, , , i]))) # sum vector of costs
 m_qaly_SoC[, i] <- rowSums(t(colSums(a_Y_u_SoC[, , , i]))) # sum vector of utilities
  
 # New Treatment:
 m_costs_NT[, i] <- rowSums(t(colSums(a_Y_c_NT[, , , i]))) # sum vector of costs
 m_qaly_NT[, i] <- rowSums(t(colSums(a_Y_u_NT[, , , i]))) # sum vector of utilities
}

# Evaluation --------------------------------------------------------------
### For Standard of Care (screening only):
sum(apply(X = m_costs_SoC, 1, mean)) # Costs
sum(apply(X = m_qaly_SoC, 1, mean)) # QALYs

### For New Treatment (screening and vaccine):
sum(apply(X = m_costs_NT, 1, mean)) # Costs
sum(apply(X = m_qaly_NT, 1, mean)) # QALYs

## Discounted costs and effects --------------------------------------------
# Discount rates:
d_e <- 0.03
d_c <- 0.03
# Discount weights for costs
v_dwc <- 1 / ((1 + d_c) ^ (0:(n_t)))
# Discount weights for effects
v_dwe <- 1 / ((1 + d_e) ^ (0:(n_t)))

# Apply discount:
v_qaly_disc_SoC <- t(m_qaly_SoC) %*% v_dwe # SoC QALYs
v_costs_disc_SoC <- t(m_costs_SoC) %*% v_dwc # SoC costs

v_qaly_disc_NT <- t(m_qaly_NT) %*% v_dwe # NT QALYs
v_costs_disc_NT <- t(m_costs_NT) %*% v_dwc # NT costs

# Summary Results ---------------------------------------------------------
### BCEA package summary ###
n_comps <- 2 # number of comparators
# Comparator names:
v_names_comps <- c("Standard of Care: screening only", "New Treatment: screening & vaccine")

# Create n_t 
m_effects <- m_costs <- matrix(NA, nrow = n.sims, ncol = n_comp, 
                           byrow = TRUE) # create an n.sims x n_comp matrix

### For Standard of Care (screening only):
m_costs[, 1] <- t(m_costs_SoC) %*% v_dwc # Costs
m_effects[, 1] <- t(m_qaly_SoC) %*% v_dwe # QALYs

### For New Treatment (screening and vaccine):
m_costs[, 2] <- t(m_costs_NT) %*% v_dwc # Costs
m_effects[, 2] <- t(m_qaly_NT) %*% v_dwe # QALYs

df_cea <- bcea(eff = m_effects, cost = m_costs, ref = 2, 
               interventions = v_names_comps, Kmax = 5724, 
               plot = FALSE)
BCEA::contour(he = df_cea)
BCEA::ceplane.plot(df_cea, wtp = 5724, graph = "ggplot")
BCEA::ceac.plot(df_cea, graph = "ggplot2")
BCEA::eib.plot(df_cea, graph = "ggplot")
BCEA::ib.plot(df_cea, wtp = 5724, graph = "ggplot")
BCEA::ce_table(he = df_cea)

# Model run time ----------------------------------------------------------
# End time counter:
end_time <- Sys.time()
# Total run time:
end_time - start_time

# End file ----------------------------------------------------------------