# ==========================================================================================
# Parameter Input Values --------------------------------------------------
# ==========================================================================================
# Note: the following script is the source code for all model input parameter values.

# Most data are derived from (Sinanovic, E., et al. 2009): "The potential cost-effectiveness of 
# adding a human papillomavirus vaccine to the cervical cancer screening programme in South 
# Africa." Most parameters are simulated using an empirical Bayes' approach, bar vaccine
# efficacy.

# Note: across all empirical bayes estimates, we assume a 5% standard error.

# ==========================================================================================
# Probability Dying less HPV for all States -------------------------------
# ==========================================================================================
# Import ASSA mortality table:
# mort_data <- read_excel("data-raw/mortality tables.xls", sheet = "final tables", range = "D1:D87")
# v_p_mort_lessHPV <- as.matrix(mort_data[1])

# Like the original study, the replication study is hugely sensitive to the HIV caused mortality 
# in South Africa in 2009. The mortality table below excludes HIV/AIDS related mortality.  
# Uncomment the code above and comment the code below to include HIV/AIDS related mortality.
# However, as in the original model, we chose to exclude HIV-related mortality in order to be 
# able to run the model for the full time-horizon.
 mort_data <- read_excel("data-raw/mortality tables.xls", sheet = "final tables", range = "F1:F87")
 v_p_mort_lessHPV <- as.matrix(mort_data[1])
 
# ==========================================================================================
# Normal/Well State Progression -------------------------------------------
# ==========================================================================================
# Normal to HPV -----------------------------------------------------------:
# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤21", "22-23", 
               "24-29", "30-49", "≥55")

# Estimated Prevalence:
Prevalence <- c(rep(0, length(0:14)), rep(0.09516258, length(15:16)), 
                rep(0.1130796, length(17:17)), 
                rep(0.139292, length(18:18)), rep(0.1563352, length(19:19)), 
                rep(0.13929220, length(20:20)), 
                rep(0.1130796, length(21:21)), rep(0.09516258, length(22:23)), 
                rep(0.04877058, length(24:29)), 
                rep(0.009950166, length(30:49)), rep(0.004987521, length(50:85)))
# log-transformation for input to JAGS:
mu.a.log <- log(Prevalence)

# ==========================================================================================
# Vaccine efficacy data ----------------------------------------------------
# ==========================================================================================
# PLACEBO (A) AND VACCINE(B).
# Number of studies:
Nstud.vac <- 11
# Positive participants in Control group:
rA.vac <- c(435, 41, 28, 38, 219, 15, 10, 21, 56, 20, 17)
# Total number of participants in control group:
nA.vac <- c(5375, 355, 277, 2239, 2924, 392, 175, 7838, 7312, 372, 2502)
# Positive participants in Trial group:
rB.vac <- c(32, 12, 1, 3, 61, 0, 0, 2, 4, 0, 1)
nB.vac <- c(5406, 366, 310, 2190, 2910, 387, 193, 7788, 7344, 401, 2497)

# ==========================================================================================
# HPV/Infection State Progression -----------------------------------------
# ==========================================================================================
# Below are the equations for relationships between all states from HPV/infection:
# Ages 15-24:
HPV_Normal_15 <- (1 - exp(-0.7 * 1.5))
# Ages 25-29:
HPV_Normal_25 <- (1 - exp(-.5 * 1.5))
# Ages ≥30:
HPV_Normal_30plus <- (1 - exp(-.15 * 1.5))

# HPV to Normal -----------------------------------------------------------
# All state progressions are conditional on not regressing to Normal. Therefore, can model
# all other parameters using one input distribution.

# Parameters for this regression:
alpha.HPVtoNormal_15to24 <- beta_params(mean = HPV_Normal_15, sigma = 0.05)$alpha
beta.HPVtoNormal_15to24 <- beta_params(mean = HPV_Normal_15, sigma = 0.05)$beta

alpha.HPVtoNormal_25to29 <- beta_params(mean = HPV_Normal_25, sigma = 0.05)$alpha
beta.HPVtoNormal_25to29 <- beta_params(mean = HPV_Normal_25, sigma = 0.05)$beta

alpha.HPVtoNormal_30toPlus <- beta_params(mean = HPV_Normal_30plus, sigma = 0.05)$alpha
beta.HPVtoNormal_30toPlus <- beta_params(mean = HPV_Normal_30plus, sigma = 0.05)$beta

# ====================================================================================
# Sub-Model for progressions from LSIL --------------------
# ====================================================================================
# To Infection or Normal ----------------------------------------
# Simulation distribution parameter values Ages 15-34:

# It is important to place an upper bound on this distribution in order to make the probabilities
# sensible due to the South African mortality context.
LSIL_15to34 <- (1 - exp(-0.65 * 6))

alpha.LSIL_15to34 <- beta_params(mean = LSIL_15to34, sigma = 0.05)$alpha
beta.LSIL_15to34 <- beta_params(mean = LSIL_15to34, sigma = 0.05)$beta

# Simulation distribution parameter values Ages 25-50:
LSIL_35to85 <- (1 - exp(-0.4 * 6))

alpha.LSIL_35to85 <- beta_params(mean = LSIL_35to85, sigma = 0.05)$alpha
beta.LSIL_35to85 <- beta_params(mean = LSIL_35to85, sigma = 0.05)$beta

# ====================================================================================
# Sub-Model for progressions from HSIL --------------------
# ====================================================================================
HSIL <- (1 - exp(-.35 * 6))

alpha.HSIL <- beta_params(mean = HSIL, sigma = 0.05)$alpha
beta.HSIL <- beta_params(mean = HSIL, sigma = 0.05)$beta

# ==========================================================================================
# All Cancer States Progression ------------------------------------------------
# ==========================================================================================
# Stage I -----------------------------------------------------------------
# Stage I to Stage II probability:
StageI.mu <- (1 - exp(-0.9 * 4))

# Distribution parameters
alpha.StageI <- beta_params(mean = StageI.mu, sigma = 0.05)$alpha
beta.StageI <- beta_params(mean = StageI.mu, sigma = 0.05)$beta

# Annual probability of symptoms from Stage I to Treatment:
StageI.Detect.mu <- 0.15

# Stage II ----------------------------------------------------------------
# Stage II to Stage III probability:
StageII.mu <- (1 - exp(-0.9 * 3))

# Distribution parameters
alpha.StageII <- beta_params(mean = StageII.mu, sigma = 0.05)$alpha
beta.StageII <- beta_params(mean = StageII.mu, sigma = 0.05)$beta

# Annual probability of symptoms from Stage II to Treatment:
StageII.Detect.mu <- 0.225

# Stage III ---------------------------------------------------------------
# Stage III to Stage IV probability:
StageIII.mu <- 1 - exp(-0.9 * 2)

# Distribution parameters
alpha.StageIII <- beta_params(mean = StageIII.mu, sigma = 0.05)$alpha
beta.StageIII <- beta_params(mean = StageIII.mu, sigma = 0.05)$beta

# Annual probability of symptoms from Stage III to Treatment:
StageIII.Detect.mu <- 0.6

# Stage IV ----------------------------------------------------------------
# Annual probability of symptoms from Stage III to Treatment: 
StageIV.Detect.mu <- 0.9

# Distribution parameters
alpha.StageIV <- beta_params(mean = StageIV.Detect.mu, sigma = 0.05)$alpha
beta.StageIV <- beta_params(mean = StageIV.Detect.mu, sigma = 0.05)$beta

# ===========================================================================================
# Stage I 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.9688
# Die = 1 - (0.9688 + mortality_lessHPV).. etc.
alpha.StageI_YearI <- beta_params(mean = 0.9688, sigma = 0.05)$alpha
beta.StageI_YearI <- beta_params(mean = 0.9688, sigma = 0.05)$beta

# Year 2
# Survive = 0.9525
# Die = 1 - (0.9525+ mortality_lessHPV)
alpha.StageI_YearII <- beta_params(mean = 0.9525, sigma = 0.05)$alpha
beta.StageI_YearII <- beta_params(mean = 0.9525, sigma = 0.05)$beta

# Year 3
# Survive = 0.9544
# Die = 1 - (0.9544 + mortality_lessHPV)
alpha.StageI_YearIII <- beta_params(mean = 0.9544, sigma = 0.05)$alpha
beta.StageI_YearIII <- beta_params(mean = 0.9544, sigma = 0.05)$beta

# Year 4
# Survive = 0.9760
# Die = 1 - 0.9760
alpha.StageI_YearIV <- beta_params(mean = 0.9760, sigma = 0.05)$alpha
beta.StageI_YearIV <- beta_params(mean = 0.9760, sigma = 0.05)$beta

# Year 5
# Survive = 0.9761
# Die = 1 - 0.9761
alpha.StageI_YearV <- beta_params(mean = 0.9761, sigma = 0.05)$alpha
beta.StageI_YearV <- beta_params(mean = 0.9761, sigma = 0.05)$beta

# 5 year survival:
((0.9688) * (0.9525) * (0.9544) * (0.9760) * (0.9761))

# ===========================================================================================
# Stage II 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.9066
# Die = 1 - 0.9066
alpha.StageII_YearI <- beta_params(mean = 0.9066, sigma = 0.05)$alpha
beta.StageII_YearI <- beta_params(mean = 0.9066, sigma = 0.05)$beta

# Year 2
# Survive = 0.8760
# Die = 1 - 0.8760
alpha.StageII_YearII <- beta_params(mean = 0.8760, sigma = 0.05)$alpha
beta.StageII_YearII <- beta_params(mean = 0.8760, sigma = 0.05)$beta

# Year 3
# Survive = 0.9225
# Die = 1 - 0.9225
alpha.StageII_YearIII <- beta_params(mean = 0.9225, sigma = 0.05)$alpha
beta.StageII_YearIII <- beta_params(mean = 0.9225, sigma = 0.05)$beta

# Year 4
# Survive = 0.9332
# Die = 1 - 0.9332
alpha.StageII_YearIV <- beta_params(mean = 0.9332, sigma = 0.05)$alpha
beta.StageII_YearIV <- beta_params(mean = 0.9332, sigma = 0.05)$beta

# Year 5
# Survive = 0.9604
# Die = 1 - 0.9604
alpha.StageII_YearV <- beta_params(mean = 0.9604, sigma = 0.05)$alpha
beta.StageII_YearV <- beta_params(mean = 0.9604, sigma = 0.05)$beta

# 5 year survival:
((0.9066) * (0.8760) * (0.9225) * (0.9332) * (0.9604))

# ===========================================================================================
# Stage III 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.7064
# Die = 1 - 0.7064
alpha.StageIII_YearI <- beta_params(mean = 0.7064, sigma = 0.05)$alpha
beta.StageIII_YearI <- beta_params(mean = 0.7064, sigma = 0.05)$beta

# Year 2
# Survive = 0.7378
# Die = 1 - 0.7378
alpha.StageIII_YearII <- beta_params(mean = 0.7378, sigma = 0.05)$alpha
beta.StageIII_YearII <- beta_params(mean = 0.7378, sigma = 0.05)$beta

# Year 3
# Survive = 0.8610
# Die = 1 - 0.8610
alpha.StageIII_YearIII <- beta_params(mean = 0.8610, sigma = 0.05)$alpha
beta.StageIII_YearIII <- beta_params(mean = 0.8610, sigma = 0.05)$beta

# Year 4
# Survive = 0.9231
# Die = 1 - 0.9231
alpha.StageIII_YearIV <- beta_params(mean = 0.9231, sigma = 0.05)$alpha
beta.StageIII_YearIV <- beta_params(mean = 0.9231, sigma = 0.05)$beta

# Year 5
# Survive = 0.9142
# Die = 1 - 0.9142
alpha.StageIII_YearV <- beta_params(mean = 0.9142, sigma = 0.05)$alpha
beta.StageIII_YearV <- beta_params(mean = 0.9142, sigma = 0.05)$beta

# 5 year survival:
((0.7064) * (0.7378) * (0.8610) * (0.9231) * (0.9142))

# ===========================================================================================
# Stage IV 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.3986
# Die = 1 - 0.3986
alpha.StageIV_YearI <- beta_params(mean = 0.3986, sigma = 0.05)$alpha
beta.StageIV_YearI <- beta_params(mean = 0.3986, sigma = 0.05)$beta

# Year 2
# Survive = 0.4982
# Die = 1 - 0.4982
alpha.StageIV_YearII <- beta_params(mean = 0.4982, sigma = 0.05)$alpha
beta.StageIV_YearII <- beta_params(mean = 0.4982, sigma = 0.05)$beta

# Year 3
# Survive = 0.7638
# Die = 1 - 0.7638
alpha.StageIV_YearIII <- beta_params(mean = 0.7638, sigma = 0.05)$alpha
beta.StageIV_YearIII <- beta_params(mean = 0.7638, sigma = 0.05)$beta

# Year 4
# Survive = 0.8652
# Die = 1 - 0.8652
alpha.StageIV_YearIV <- beta_params(mean = 0.8652, sigma = 0.05)$alpha
beta.StageIV_YearIV <- beta_params(mean = 0.8652, sigma = 0.05)$beta

# Year 5
# Survive = 0.8592
# Die = 1 - 0.8592
alpha.StageIV_YearV <- beta_params(mean = 0.8592, sigma = 0.05)$alpha
beta.StageIV_YearV <- beta_params(mean = 0.8592, sigma = 0.05)$beta

# 5 year survival:
((0.3986) * (0.4982) * (0.7638) * (0.8652) * (0.8592))

# End file ----------------------------------------------------------------