# Note: the following script is the source code for all model input parameter values.


# Informative prior -------------------------------------------------------
# Most data is derived from:
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# See Collated Studies.R file. for vaccine RCT data.

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
length(Prevalence)
mu.a.log <- log(Prevalence)

# Optional 'emprical bayes' method:
# mu.a.log <- lnorm_params(m = Prevalence, v = .01)$mu
# sigma.a.log <- lnorm_params(m = Prevalence, v = .01)$sigma
# prec.age <- 1 / (sigma.a.log * sigma.a.log)
# prec.age

# Normal to Death ---------------------------------------------------------
# Import ASSA mortality table:
mort_data <- read_excel(
"/Users/joshuamusson/Desktop/Analytics/R/Dissertation/Evidence_Synthesis/mortality tables.xls", 
                        sheet = "ASSA data", range = "A3:H94")
mort_data

# Total Pop. divided by total deaths:
v_p_mort_All <- mort_data[1:86, 3] / mort_data[1:86, 2]
v_p_mort_All

# Probability of mortality Cervical Cancer from HPV:
v_p_mort_HPV <- mort_data[1:86, 8]
colnames(v_p_mort_HPV) <- "Death_HPV"
v_p_mort_HPV

# Probability of mortality less HPV:
v_p_mort_lessHPV <- as.matrix(v_p_mort_All[1:86, ] - v_p_mort_HPV[1:86, ])
colnames(v_p_mort_lessHPV) <- c("Death_less.HPV")
v_p_mort_lessHPV

# annual.mortality_1to14 <- v_p_mort_lessHPV[1:15, 2]
# Ages 15:24
# annual.mortality_15to24 <- v_p_mort_lessHPV[16:25, 2]
# Ages 25:29
# annual.mortality_25to29 <- v_p_mort_lessHPV[26:30, 2]
# Ages ≥30:
# annual.mortality_30plus <- v_p_mort_lessHPV[31:86, 2]

# State equations ---------------------------------------------------------
# Well to death:
# v_p_mort_lessHPV[1:86]
# Well to Infection:
# ((1 - v_p_mort_lessHPV[1:86]) * omega.age[, 1:86])
# Well to Well:
# (1 - (v_p_mort_lessHPV[1:86] + ((1 - v_p_mort_lessHPV[1:86]) * omega.age[, 1:86])))

# LOTP check:
# (1 - (v_p_mort_lessHPV[1:86] + ((1 - v_p_mort_lessHPV[1:86]) * omega.age[, 1:86]))) + (
#  (1 - v_p_mort_lessHPV[1:86]) * omega.age[, 1:86]) + v_p_mort_lessHPV[1:86]


# ==========================================================================================
# Vaccine efficacy data ----------------------------------------------------
# ==========================================================================================
# PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 11
rA.vac <- c(435, 41, 28, 38, 219, 15, 10, 21, 56, 20, 17)
nA.vac <- c(5375, 355, 277, 2239, 2924, 392, 175, 7838, 7312, 372, 2502)
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
alpha.HPVtoNormal_15to24 <- beta_params(mean = HPV_Normal_15, sigma = 0.025)$alpha
beta.HPVtoNormal_15to24 <- beta_params(mean = HPV_Normal_15, sigma = 0.025)$beta

alpha.HPVtoNormal_25to29 <- beta_params(mean = HPV_Normal_25, sigma = 0.025)$alpha
beta.HPVtoNormal_25to29 <- beta_params(mean = HPV_Normal_25, sigma = 0.025)$beta

alpha.HPVtoNormal_30toPlus <- beta_params(mean = HPV_Normal_30plus, sigma = 0.025)$alpha
beta.HPVtoNormal_30toPlus <- beta_params(mean = HPV_Normal_30plus, sigma = 0.025)$beta


# State Equations ---------------------------------------------------------
# HPV to Death ------------------------------------------------------------
# Death independent on regression to normal...
# Ages 15-24:
HPV_Death_15 <- v_p_mort_lessHPV[15:24]
HPV_Death_15

# Ages 25-29:
HPV_Death_25 <- v_p_mort_lessHPV[25:29]
HPV_Death_25

# Ages ≥30:
HPV_Death_30plus <- v_p_mort_lessHPV[30:86]
HPV_Death_30plus

# The complement set of Normal by age group:
# Ages 15-24:
(1 - (HPV_Normal_15 + HPV_Death_15))

# Ages 25-29:
1 - HPV_Normal_25

# Ages ≥30:
1 - HPV_Normal_30plus


# HPV to LSIL -------------------------------------------------------------
# Ages 15-24:
HPV_LSIL_15 <- (1 - (HPV_Normal_15 + HPV_Death_15)) * (
 (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))
HPV_LSIL_15

# Ages 25-29:
HPV_LSIL_25 <- (1 - (HPV_Normal_25 + HPV_Death_25)) * (
 (1 - exp(-.2 * 3)) - ((1 - exp(-.2 * 3)) * .1))
HPV_LSIL_25

# Ages ≥30:
HPV_LSIL_30plus <- (1 - (HPV_Normal_30plus + HPV_Death_30plus))  * ((1 - exp(-.2 * 3)) - 
                                               ((1 - exp(-.2 * 3)) * .1))
HPV_LSIL_30plus

# HPV to HSIL -------------------------------------------------------------
# Ages 15-24:
HPV_HSIL_15 <- ((1 - (HPV_Normal_15 + HPV_Death_15)) * ((1 - exp(-.2 * 3)) * .1))
HPV_HSIL_15

# Ages 25-29:
HPV_HSIL_25 <- ((1 - (HPV_Normal_25 + HPV_Death_25)) * ((1 - exp(-.2 * 3)) * .1))
HPV_HSIL_25

# Ages ≥30:
HPV_HSIL_30plus <- ((1 - (HPV_Normal_30plus + HPV_Death_30plus)) * ((1 - exp(-.2 * 3)) * .1))
HPV_HSIL_30plus

# HPV to HPV --------------------------------------------------------------
# Ages 15-24:
HPV_HPV_15 <- (1 - (HPV_Normal_15 + HPV_Death_15)) - (HPV_LSIL_15 + 
                                HPV_HSIL_15)

# Ages 25-29:
HPV_HPV_25 <- (1 - (HPV_Normal_25 + HPV_Death_25)) - (HPV_LSIL_25 + 
                                HPV_HSIL_25)

# Ages ≥30:
HPV_HPV_30up <- (1 - (HPV_Normal_30plus + HPV_Death_30plus)) - (HPV_LSIL_30plus + 
                                HPV_HSIL_30plus)

# State Calibration -------------------------------------------------------------
# Calibration - Ages 15-34. Must equal 1 across all ages:
HPV_HPV_15 + HPV_HSIL_15 + HPV_LSIL_15 + HPV_Death_15 + HPV_Normal_15

# Calibration - Ages 25-29. Must equal 1 across all ages:
HPV_HPV_25 + HPV_HSIL_25 + HPV_LSIL_25 + HPV_Death_25 + HPV_Normal_25

# Calibration - Ages ≥30. Must equal 1 across all ages:
HPV_HPV_30up + HPV_HSIL_30plus + HPV_LSIL_30plus + HPV_Death_30plus + HPV_Normal_30plus

# ====================================================================================
# Sub-Model for progressions from LSIL --------------------
# ====================================================================================
# To Infection or Normal ----------------------------------------
# Simulation distribution parameter values Ages 15-34:

# It is important to place an upper bound on this distribution in order to make the probabilities
# sensible. Truncation is placed on either distribution as T(0, 0.75-0.90).
LSIL_15to34 <- (1 - exp(-0.65 * 6))

alpha.LSIL_15to34 <- beta_params(mean = LSIL_15to34, sigma = 0.025)$alpha
beta.LSIL_15to34 <- beta_params(mean = LSIL_15to34, sigma = 0.025)$beta

# Simulation distribution parameter values Ages 25-50:
LSIL_35to85 <- (1 - exp(-0.4 * 6))

alpha.LSIL_35to85 <- beta_params(mean = LSIL_35to85, sigma = 0.025)$alpha
beta.LSIL_35to85 <- beta_params(mean = LSIL_35to85, sigma = 0.025)$beta

# State Equations ---------------------------------------------------------
# Transition to Death:
LSILtoDeath <- (v_p_mort_lessHPV[16:35])
# To Normal:
LSILNORM <-  (LSIL_15to24 * 0.9)
LSILNORM
# To Infection:
LSILINF <- (LSIL_15to24) -  (LSIL_15to24 * 0.9)
LSILINF

# Transition to HSIL
LSILtoHSIL <- (1 - LSIL_15to24) * (1 - exp(-0.1 * 6))
LSILtoHSIL

# Transition to LSIL
LSILtoLSIL <- 1 - (LSILtoDeath + LSILNORM + LSILINF + LSILtoHSIL)
LSILtoLSIL
# LOTP check:
LSILtoDeath + LSILNORM + LSILINF + LSILtoHSIL + LSILtoLSIL

# ====================================================================================
# Sub-Model for progressions from HSIL --------------------
# ====================================================================================
HSIL <- (1 - exp(-.35 * 6))
HSIL
alpha.HSIL <- beta_params(mean = HSIL, sigma = 0.025)$alpha
beta.HSIL <- beta_params(mean = HSIL, sigma = 0.025)$beta


# State Equations ---------------------------------------------------------
# Transition to Normal:
HSILtoNormal <- ((1 - exp(-.35 * 6)) * .5)
# Transition to LSIL:
HSILtoLSIL <- (1 - exp(-.35 * 6)) - ((1 - exp(-.35 * 6)) * .5)
# Transition to Stage I Cancer:
HSILtoStageI <- ((1 - HSIL) * (1 - exp(-0.4 * 10)))


# Transition to HSIL
HSILtoHSIL <- 1 - (HSILtoNormal + HSILtoLSIL + HSILtoStageI)
HSILtoHSIL

# LOTP check:
HSILtoHSIL + HSILtoLSIL + HSILtoNormal + HSILtoStageI

# ==========================================================================================
# All Cancer States Progression ------------------------------------------------
# ==========================================================================================
# Stage I -----------------------------------------------------------------
# Stage I to Stage II probability:
StageI.mu <- (1 - exp(-0.9 * 4))

# Distribution parameters
alpha.StageI <- beta_params(mean = StageI.mu, sigma = 0.025)$alpha
beta.StageI <- beta_params(mean = StageI.mu, sigma = 0.025)$beta

# Annual probability of symptoms from Stage I to Treatment:
StageI.Detect.mu <- 0.15

# Stage II ----------------------------------------------------------------
# Stage II to Stage III probability:
StageII.mu <- (1 - exp(-0.9 * 3))

# Distribution parameters
alpha.StageII <- beta_params(mean = StageII.mu, sigma = 0.025)$alpha
beta.StageII <- beta_params(mean = StageII.mu, sigma = 0.025)$beta

# Annual probability of symptoms from Stage II to Treatment:
StageII.Detect.mu <- 0.225

# Stage III ---------------------------------------------------------------
# Stage III to Stage IV probability:
StageIII.mu <- 1 - exp(-0.9 * 2)

# Distribution parameters
alpha.StageIII <- beta_params(mean = StageIII.mu, sigma = 0.025)$alpha
beta.StageIII <- beta_params(mean = StageIII.mu, sigma = 0.025)$beta

# Annual probability of symptoms from Stage III to Treatment:
StageIII.Detect.mu <- 0.6

# Stage IV ----------------------------------------------------------------
# Annual probability of symptoms from Stage III to Treatment: 
StageIV.Detect.mu <- 0.9

# Distribution parameters
alpha.StageIV <- beta_params(mean = StageIV.Detect.mu, sigma = 0.025)$alpha
beta.StageIV <- beta_params(mean = StageIV.Detect.mu, sigma = 0.025)$beta

# ===========================================================================================
# Stage I 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.9688
# Die = 1 - (0.9688 + mortality_lessHPV).. etc.
alpha.StageI_YearI <- beta_params(mean = 0.9688, sigma = 0.025)$alpha
beta.StageI_YearI <- beta_params(mean = 0.9688, sigma = 0.025)$beta

# Year 2
# Survive = 0.9525
# Die = 1 - (0.9525+ mortality_lessHPV)
alpha.StageI_YearII <- beta_params(mean = 0.9525, sigma = 0.025)$alpha
beta.StageI_YearII <- beta_params(mean = 0.9525, sigma = 0.025)$beta

# Year 3
# Survive = 0.9544
# Die = 1 - (0.9544 + mortality_lessHPV)
alpha.StageI_YearIII <- beta_params(mean = 0.9544, sigma = 0.025)$alpha
beta.StageI_YearIII <- beta_params(mean = 0.9544, sigma = 0.025)$beta

# Year 4
# Survive = 0.9760
# Die = 1 - 0.9760
alpha.StageI_YearIV <- beta_params(mean = 0.9760, sigma = 0.025)$alpha
beta.StageI_YearIV <- beta_params(mean = 0.9760, sigma = 0.025)$beta

# Year 5
# Survive = 0.9761
# Die = 1 - 0.9761
alpha.StageI_YearV <- beta_params(mean = 0.9761, sigma = 0.025)$alpha
beta.StageI_YearV <- beta_params(mean = 0.9761, sigma = 0.025)$beta

# 5 year survival:
((0.9688) * (0.9525) * (0.9544) * (0.9760) * (0.9761))

# ===========================================================================================
# Stage II 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.9066
# Die = 1 - 0.9066
alpha.StageII_YearI <- beta_params(mean = 0.9066, sigma = 0.025)$alpha
beta.StageII_YearI <- beta_params(mean = 0.9066, sigma = 0.025)$beta

# Year 2
# Survive = 0.8760
# Die = 1 - 0.8760
alpha.StageII_YearII <- beta_params(mean = 0.8760, sigma = 0.025)$alpha
beta.StageII_YearII <- beta_params(mean = 0.8760, sigma = 0.025)$beta

# Year 3
# Survive = 0.9225
# Die = 1 - 0.9225
alpha.StageII_YearIII <- beta_params(mean = 0.9225, sigma = 0.025)$alpha
beta.StageII_YearIII <- beta_params(mean = 0.9225, sigma = 0.025)$beta

# Year 4
# Survive = 0.9332
# Die = 1 - 0.9332
alpha.StageII_YearIV <- beta_params(mean = 0.9332, sigma = 0.025)$alpha
beta.StageII_YearIV <- beta_params(mean = 0.9332, sigma = 0.025)$beta

# Year 5
# Survive = 0.9604
# Die = 1 - 0.9604
alpha.StageII_YearV <- beta_params(mean = 0.9604, sigma = 0.025)$alpha
beta.StageII_YearV <- beta_params(mean = 0.9604, sigma = 0.025)$beta

# 5 year survival:
((0.9066) * (0.8760) * (0.9225) * (0.9332) * (0.9604))

# ===========================================================================================
# Stage III 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.7064
# Die = 1 - 0.7064
alpha.StageIII_YearI <- beta_params(mean = 0.7064, sigma = 0.025)$alpha
beta.StageIII_YearI <- beta_params(mean = 0.7064, sigma = 0.025)$beta

# Year 2
# Survive = 0.7378
# Die = 1 - 0.7378
alpha.StageIII_YearII <- beta_params(mean = 0.7378, sigma = 0.025)$alpha
beta.StageIII_YearII <- beta_params(mean = 0.7378, sigma = 0.025)$beta

# Year 3
# Survive = 0.8610
# Die = 1 - 0.8610
alpha.StageIII_YearIII <- beta_params(mean = 0.8610, sigma = 0.025)$alpha
beta.StageIII_YearIII <- beta_params(mean = 0.8610, sigma = 0.025)$beta

# Year 4
# Survive = 0.9231
# Die = 1 - 0.9231
alpha.StageIII_YearIV <- beta_params(mean = 0.9231, sigma = 0.025)$alpha
beta.StageIII_YearIV <- beta_params(mean = 0.9231, sigma = 0.025)$beta

# Year 5
# Survive = 0.9142
# Die = 1 - 0.9142
alpha.StageIII_YearV <- beta_params(mean = 0.9142, sigma = 0.025)$alpha
beta.StageIII_YearV <- beta_params(mean = 0.9142, sigma = 0.025)$beta

# 5 year survival:
((0.7064) * (0.7378) * (0.8610) * (0.9231) * (0.9142))

# ===========================================================================================
# Stage IV 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.3986
# Die = 1 - 0.3986
alpha.StageIV_YearI <- beta_params(mean = 0.3986, sigma = 0.025)$alpha
beta.StageIV_YearI <- beta_params(mean = 0.3986, sigma = 0.025)$beta

# Year 2
# Survive = 0.4982
# Die = 1 - 0.4982
alpha.StageIV_YearII <- beta_params(mean = 0.4982, sigma = 0.025)$alpha
beta.StageIV_YearII <- beta_params(mean = 0.4982, sigma = 0.025)$beta

# Year 3
# Survive = 0.7638
# Die = 1 - 0.7638
alpha.StageIV_YearIII <- beta_params(mean = 0.7638, sigma = 0.025)$alpha
beta.StageIV_YearIII <- beta_params(mean = 0.7638, sigma = 0.025)$beta

# Year 4
# Survive = 0.8652
# Die = 1 - 0.8652
alpha.StageIV_YearIV <- beta_params(mean = 0.8652, sigma = 0.025)$alpha
beta.StageIV_YearIV <- beta_params(mean = 0.8652, sigma = 0.025)$beta

# Year 5
# Survive = 0.8592
# Die = 1 - 0.8592
alpha.StageIV_YearV <- beta_params(mean = 0.8592, sigma = 0.025)$alpha
beta.StageIV_YearV <- beta_params(mean = 0.8592, sigma = 0.025)$beta

# 5 year survival:
((0.3986) * (0.4982) * (0.7638) * (0.8652) * (0.8592))

# State equations ---------------------------------------------------------
# Death from cancer, for example:
# 1 - (surv.StageI_year1 - v_p_mort_lessHPV[1, 2])

# Death from cancer, for example:
# 1 - (surv.StageI_year1 - v_p_mort_lessHPV[1, 2])

# Stage I-IV equations for example DELETE WHEN FINISHED:
# StageItoStageII <- (Stage.I.canc - (Stage.I.canc * 0.15))
# StageItoStageII

# StageItoTreat <- (Stage.I.canc * 0.15)
# StageItoTreat

# StageItoDeath <- v_p_mort_lessHPV[22]

# StageItoStageI <- 1 - (StageItoDeath + StageItoTreat + StageItoStageII)

# LOTP Check:
# (StageItoDeath + StageItoTreat + StageItoStageII + StageItoStageI)

# End file ----------------------------------------------------------------