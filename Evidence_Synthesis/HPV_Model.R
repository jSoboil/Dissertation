library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(parallel)

options(mc.cores = detectCores())

# ==========================================================================================
# Vaccine efficacy data ----------------------------------------------------
# ==========================================================================================

# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 7
rA.vac <- c(41, 46, 38, 53, 24, 60, 42)
nA.vac <- c(1607, 233, 1583, 1245, 422, 2279, 5260)
rB.vac <- c(4, 2, 1, 5, 3, 0, 1)
nB.vac <- c(1615, 235, 1578, 1276, 419, 2261, 5305)

# ===========================================================================================
# Vaccine coverage --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:
mu.vac_alpha <- beta_params(mean = .8, sigma = .1)$alpha
mu.vac_beta <- beta_params(mean = .8, sigma = .1)$beta
# Assumes roughly 80% coverage.

# ===========================================================================================
# Vaccine compliance --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# Assume roughly 100% compliance:
iota.comp_alpha <- beta_params(mean = .9999, sigma = .009)$alpha
iota.comp_beta <- beta_params(mean = .9999, sigma = .009)$beta

# ===========================================================================================
# Annual Screening Coverage  --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. HPV and Related Diseases Report
# Proportion of pop. screening coverage every three years according to age group, converted 
# to an annual rate, and then annual probability:

sigma.screening_alpha <- c(rep(beta_params(mean = .05, sigma = .05)$alpha, length(15:29)),
                           rep(beta_params(mean = .08, sigma = .05)$alpha, length(30:39)),
                           rep(beta_params(mean = .04, sigma = .05)$alpha, length(40:49)),
                           rep(beta_params(mean = .03,sigma = .05)$alpha, length(50:59)),
                           rep(beta_params(mean = .02, sigma = .05)$alpha, length(60:90))
                           )
sigma.screening_alpha



sigma.screening_beta <- c(rep(beta_params(mean = .05, sigma = .05)$beta, length(15:29)),
                           rep(beta_params(mean = .08, sigma = .05)$beta, length(30:39)),
                           rep(beta_params(mean = .04, sigma = .05)$beta, length(40:49)),
                           rep(beta_params(mean = .03,sigma = .05)$beta, length(50:59)),
                           rep(beta_params(mean = .02, sigma = .05)$beta, length(60:90))
                           )
sigma.screening_beta

# ==========================================================================================
# Age-specific infection ----------------------------------------------------
# ==========================================================================================
# Note: will use a relative risk parameter for proportion of HIV+ population who have 
# increased risk of being reinfected rather than clearing HPV. See Italian study.

# Informative prior -------------------------------------
# Mix of the following two sources:
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# 2. Data collected from HPV and Related Diseases Report. Used age groups 25+.
1 - exp(-.005)
# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤25", "25-34", 
               "35-44", "45-54", "55-64")

# Estimated Prevalence:
Prevalence <- c(rep(0.09516258, length(15:16)), rep(0.1130796, length(17:17)), 
                rep(0.139292, length(18:18)), rep(0.1563352, length(19:19)), 
                rep(0.139292, length(20:20)), rep(.435, length(21:24)), 
                rep(.3643, length(25:34)), rep(0.2051, length(35:44)), 
                rep(0.1852, length(45:54)), rep(0.1654, length(55:64)), 
                rep(0.004987521, length(65:90)))

omega.mu.log <- lnorm_params(m = Prevalence, v = .05)$mu
omega.sigma.log <- 1 / (lnorm_params(m = Prevalence, v = .05)$sigma) ^ 2
omega.mu.log
omega.sigma.log
cbind(Age = 15:90, Prevalence)

# Prior sampling model:
# omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])

# ==========================================================================================
# Age-specific Regression of Infection to Exposure -------------------------------
# ==========================================================================================
# Informative prior:

# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

delta.regr_alpha <- c(
  # Ages 15-24:
  rep(beta_params(mean = .719, sigma = .05)$alpha, length(15:24)), 
  # Ages 25-29:
  rep(beta_params(mean = .699, sigma = .05)$alpha, length(25:29)),
  # Ages 30-39:
  rep(beta_params(mean = .35, sigma = .05)$alpha, length(30:39)),
  # Ages 40-49:
  rep(beta_params(mean = .201, sigma = .05)$alpha, length(40:49)),
  # Ages 50:90:
  rep(beta_params(mean = .099, sigma = .05)$alpha, length(50:90))
  )

delta.regr_beta <- c(
  # Ages 15-24:
  rep(beta_params(mean = .719, sigma = .05)$beta, length(15:24)), 
  # Ages 25-29:
  rep(beta_params(mean = .699, sigma = .05)$beta, length(25:29)),
  # Ages 30-39:
  rep(beta_params(mean = .35, sigma = .05)$beta, length(30:39)),
  # Ages 40-49:
  rep(beta_params(mean = .201, sigma = .05)$beta, length(40:49)),
  # Ages 50:90:
  rep(beta_params(mean = .099, sigma = .05)$beta, length(50:90))
  )


# ==========================================================================================
# LSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From Infection progression to Low-Grade Squamous Intraepithelial Lesions:
delta.LSIL_alpha <- beta_params(mean = .11, sigma = .1)$alpha
delta.LSIL_beta <- beta_params(mean = .11, sigma = .1)$beta

# ==========================================================================================
# HSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From Infection progression to High-Grade Squamous Intraepithelial Lesions:
delta.HSIL_alpha <- beta_params(mean = .022, sigma = .01)$alpha
delta.HSIL_beta <- beta_params(mean = .022, sigma = .01)$beta

# ===========================================================================================
# Subject to No-Treatment  --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1.Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From LSIL to Clearance:
pi.I.LSIL_Clear_alpha <- beta_params(mean = .71, sigma = .05)$alpha
pi.I.LSIL_Clear_beta <- beta_params(mean = .71, sigma = .05)$beta
# From LSIL to HSIL:
pi.I.LSIL_HSIL_alpha <- beta_params(mean = .224, sigma = .05)$alpha
pi.I.LSIL_HSIL_beta <- beta_params(mean = .224, sigma = .05)$beta
# From HSIL to Clearance:
pi.I.HSIL_Clear_alpha <- beta_params(mean = .355, sigma = .05)$alpha
pi.I.HSIL_Clear_beta <- beta_params(mean = .355, sigma = .05)$beta
# From HSIL to LSIL:
pi.I.HSIL_LSIL_alpha <- beta_params(mean = .25, sigma = .05)$alpha
pi.I.HSIL_LSIL_beta <- beta_params(mean = .25, sigma = .05)$beta
# From HSIL to Cancer:
pi.I.HSIL_Cancer_alpha <- beta_params(mean = .05, sigma = .05)$alpha
pi.I.HSIL_Cancer_beta <- beta_params(mean = .05, sigma = .05)$beta

# ===========================================================================================
# Subject to Treatment (Conization)  --------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From LSIL to Clearance:
beta_params(mean = .8990, sigma = .05)
# From LSIL to HSIL:
beta_params(mean = .12, sigma = .05)
# From HSIL to Clearance:
beta_params(mean = .861, sigma = .05)
# From HSIL to cancer:
beta_params(mean = .015, sigma = .05)

# Probability of Conization
# Immediate:
beta_params(mean = .302, sigma = .05)
# Delayed:
beta_params(mean = .17, sigma = .05)

# ==========================================================================================
# Cervical Adenocarcinoma -------------------------------------------------
# ==========================================================================================

# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# Proportion of inviduals in the four stages who progress to cancer:
# FIGO I 
# dDirc(.5500)
# FIGO II
# dDirc(.1500)
# FIGO III
# dDirc(.1200)
# FIGO IV
# dDirc(.1800)
cancer_alpha <- list(.55, .15, .12, .18)

# Survival probabilities for each FIGO stage, over 5 years.

# Year 1 survival by stage:
# FIGO I
alphaI.year_1 <- beta_params(mean = .9770, sigma = .05)$alpha
betaI.year_1 <- beta_params(mean = .9770, sigma = .05)$beta
# FIGO II
alphaII.year_1 <- beta_params(mean = .8290, sigma = .05)$alpha
betaII.year_1 <- beta_params(mean = .8290, sigma = .05)$beta
# FIGO III
alphaIII.year_1 <- beta_params(mean = .59, sigma = .05)$alpha
betaIII.year_1 <- beta_params(mean = .59, sigma = .05)$beta
# FIGO IV
alphaIV.year_1 <- beta_params(mean = .5020, sigma = .05)$alpha
betaIV.year_1 <- beta_params(mean = .5020, sigma = .05)$beta

# Year 2 survival by stage:
# FIGO I
alphaI.year_2 <- beta_params(mean = .9790, sigma = .05)$alpha
betaI.year_2 <- beta_params(mean = .9790, sigma = .05)$beta
# FIGO II
alphaII.year_2 <- beta_params(mean = .8330, sigma = .05)$alpha
betaII.year_2 <- beta_params(mean = .8330, sigma = .05)$beta
# FIGO III
alphaIII.year_2 <- beta_params(mean = .6930, sigma = .05)$alpha
betaIII.year_2 <- beta_params(mean = .6930, sigma = .05)$beta
# FIGO IV
alphaIV.year_2 <- beta_params(mean = .7820, sigma = .05)$alpha
betaIV.year_2 <- beta_params(mean = .7820, sigma = .05)$beta

# Year 3 survival by stage:
# FIGO I
alphaI.year_3 <- beta_params(mean = .9630, sigma = .05)$alpha
betaI.year_3 <- beta_params(mean = .9630, sigma = .05)$beta
# FIGO II
alphaII.year_3 <- beta_params(mean = .7550, sigma = .05)$alpha
betaII.year_3 <- beta_params(mean = .7550, sigma = .05)$beta
# FIGO III
alphaIII.year_3 <- beta_params(mean = .7780, sigma = .05)$alpha
betaIII.year_3 <- beta_params(mean = .7780, sigma = .05)$beta
# FIGO IV
alphaIV.year_3 <- beta_params(mean = .7220, sigma = .05)$alpha
betaIV.year_3 <- beta_params(mean = .7220, sigma = .05)$beta

# Year 4 survival by stage:
# FIGO I
alphaI.year_4 <- beta_params(mean = .9890, sigma = .05)$alpha
betaI.year_4 <- beta_params(mean = .9890, sigma = .05)$beta
# FIGO II
alphaII.year_4 <- beta_params(mean = .8690, sigma = .05)$alpha
betaII.year_4 <- beta_params(mean = .8690, sigma = .05)$beta
# FIGO III
alphaIII.year_4 <- beta_params(mean = .9290, sigma = .05)$alpha
betaIII.year_4 <- beta_params(mean = .9290, sigma = .05)$beta
# FIGO IV
alphaIV.year_4 <- beta_params(mean = .9250, sigma = .05)$alpha
betaIV.year_4 <- beta_params(mean = .9250, sigma = .05)$beta

# ===========================================================================================
# Relative Risk Increase (HIV+) ----------------------------------------------------
# ===========================================================================================

# Study: Mbulawa Z., et al. 2015 ------------------------------------------
# Human papillomavirus prevalence in South African women and men according to age and human 
# immunodeficiency virus status.

# Then use proportion of HIV+ pop. for how many at risk at time j. See Italian study.

# Participants were 208 HIV-negative women, 278 HIV-positive women, 325 HIV-negative men and
# 161 HIV-positive men between the ages of 18–66 years. HPV types were determined in cervical
# and penile cells by Roche Linear Array HPV genotyping assay.

# HIV+ women positive to HPV: 205/277
# HIV- women positive to HPV: 76/207

rHIV_pos <- 205
nHIV_pos <- 277
rHIV_neg <- 76
nHIV_neg <- 207

# Proportion of pop. at increased risk due to HIV:
zeta_alpha <- beta_params(mean = .131, sigma = .05)$alpha
zeta_beta <- beta_params(mean = .131, sigma = .05)$beta

# ==========================================================================================
# Model ------------------------------------------------------
# ==========================================================================================
model_String <- "
model {
  
# SUB-MODEL 1: VACCINE-EFFICACY.
  # model parameters abbreviated by .vac
  for (i in 1:Nstud.vac) {
    # Likelihood
    rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
    rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])
    # Logistic function:
    logit(pA.vac[i]) <- mu.vac[i]
    logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
    
    # Average effect:
    mu.vac[i] ~ dnorm(0, 1.0e-2)
    # Random population effect:
    delta.vac[i] ~ dnorm(psi.vac, prec.vac)
 }
  # Priors for sub-model 1:
   psi.vac ~ dnorm(0, 1.0e-1)
   tau.vac ~ dunif(0, 1)
   prec.vac <- 1 / pow(tau.vac, 2)
   # Odds ratio:
   OR.vac <- exp(psi.vac)
   # Probability of no protection given
   # vaccine:
   pEfficacy.vac <- (OR.vac / 1 + OR.vac)

# SUB-MODEL 2: AVERAGE VACCINE COVERAGE.
  # model parameters abbreviated by mu.vac
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
  mu.vac_coverage ~ dbeta(mu.vac_alpha, mu.vac_beta)

# SUB-MODEL 3: VACCINE COMPLIANCE.
  # model parameters abbreviated by iota.
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
  iota.compli ~ dbeta(iota.comp_alpha, iota.comp_beta)
 
# SUB-MODEL 4: CROSS-PROTECTION EFFECT.
  # model parameters abbreviated by chi.
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
  chi.xProtect ~ dlnorm(chi.crossPro_mu, chi.crossPro_sigma)

# SUB-MODEL 5: SCREENING PROBABILITY BY AGE
# GROUP.
  # model parameters abbreviated by .screen
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
 for (i in 1:76) {
  sigma.screenCov[i] ~ dbeta(sigma.screening_alpha[i], sigma.screening_beta[i])
 }

# SUB-MODEL 6: AGE-SPECIFIC INFECTION.
  # model parameters abbreviated by .age
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
 for (i in 1:76) {
   omega.age[i] ~ dlnorm(omega.mu.log[i], omega.sigma.log[i])
 }
  
# SUB-MODEL 7: AGE-SPECIFIC REGRESSION OF
# INFECTION.
  # model parameters abbreviated by .regr
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
 for (i in 1:76) {
  delta.regr[i] ~ dbeta(delta.regr_alpha[i], delta.regr_beta[i])  
 }
 
# SUB-MODEL 8: PROGRESSION TO LSIL OR
# HSIL FROM INFECTION.
  # model parameters abbreviated by .regr
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
 delta.LSIL ~ dbeta(delta.LSIL_alpha, delta.LSIL_beta)
 delta.HSIL ~ dbeta(delta.HSIL_alpha, delta.HSIL_beta)

# SUB-MODEL 9: PROGRESSION AND REGRESSION FROM
# LSIL AND HSIL GIVEN NO TREATMENT.
  # model parameters abbreviated by pi.I.
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior:
    
 pi.I.LSIL_Clear ~ dbeta(pi.I.LSIL_Clear_alpha, pi.I.LSIL_Clear_beta)
 pi.I.LSIL_HSIL ~ dbeta(pi.I.LSIL_HSIL_alpha, pi.I.LSIL_HSIL_beta)
 pi.I.HSIL_Clear ~ dbeta(pi.I.HSIL_Clear_alpha, pi.I.HSIL_Clear_beta)
 pi.I.HSIL_LSIL ~ dbeta(pi.I.HSIL_LSIL_alpha, pi.I.HSIL_LSIL_beta)
 pi.I.HSIL_Cancer ~ dbeta(pi.I.HSIL_Cancer_alpha, pi.I.HSIL_Cancer_beta)

# SUB-MODEL 10: CERVICAL CANCER.
# the distribution of invidiuals who have 
# cancer according to FIGO stages I-IV, and 
# the associated survival probabilities according 
# to each stage, over 4 years.

  # Note: equivalent to standard PSA, as it is 
  # sampling direclty from the prior.
    
# Distribution of Cervical Cancer stages: 
 FIGO ~ ddirch(cancer_alpha)
 
# Surivival probabilities:
 # Year 1:
 FIGO_I.year_1 ~ dbeta(alphaI.year_1, betaI.year_1)
 FIGO_II.year_1 ~ dbeta(alphaII.year_1, betaII.year_1)
 FIGO_III.year_1 ~ dbeta(alphaIII.year_1, betaIII.year_1)
 FIGO_IV.year_1 ~ dbeta(alphaIV.year_1, betaIV.year_1)
 
 # Year 2:
 FIGO_I.year_2 ~ dbeta(alphaI.year_2, betaI.year_2)
 FIGO_II.year_2 ~ dbeta(alphaII.year_2, betaII.year_2)
 FIGO_III.year_2 ~ dbeta(alphaIII.year_2, betaIII.year_2)
 FIGO_IV.year_2 ~ dbeta(alphaIV.year_2, betaIV.year_2)
 
 # Year 3:
 FIGO_I.year_3 ~ dbeta(alphaI.year_3, betaI.year_3)
 FIGO_II.year_3 ~ dbeta(alphaII.year_3, betaII.year_3)
 FIGO_III.year_3 ~ dbeta(alphaIII.year_3, betaIII.year_3)
 FIGO_IV.year_3 ~ dbeta(alphaIV.year_3, betaIV.year_3)
 
 # Year 4:
 FIGO_I.year_4 ~ dbeta(alphaI.year_4, betaI.year_4)
 FIGO_II.year_4 ~ dbeta(alphaII.year_4, betaII.year_4)
 FIGO_III.year_4 ~ dbeta(alphaIII.year_4, betaIII.year_4)
 FIGO_IV.year_4 ~ dbeta(alphaIV.year_4, betaIV.year_4)
 
# SUB-MODEL 11: POPULATION AT INCREASED RISK 
# DUE TO HIV+.
  # model parameters abbreviated by HIV.
    
 # Binomial Likelihood:
  rHIV_pos ~ dbin(pHIV_pos, nHIV_pos)
  rHIV_neg ~ dbin(pHIV_neg, nHIV_neg)
  # Logistic function:
  logit(pHIV_pos) <- mu.hiv + delta.hiv
  logit(pHIV_neg) <- mu.hiv
 
   # Average effect:
   mu.hiv ~ dnorm(0, 1 / 10)
   # Random effect:
   delta.hiv ~ dnorm(nu.hiv, prec.hiv)
 
  # Hyperpriors on delta.hiv:
  nu.hiv ~ dnorm(0, 1 / 10)
  prec.hiv <- 1 / psi.hiv * psi.hiv
  psi.hiv ~ dunif(0, 1)
 
  HIV_odds <- exp(nu.hiv)
  HIV_RR <- psi.hiv / (1 + psi.hiv)
 
  # Proportion of population at
  # increased risk:
  zeta.HIVpop ~ dbeta(zeta_alpha, zeta_beta)
 
}
"
writeLines(text = model_String, con = "HPV_model.txt")

# Compiled data for JAGS input --------------------------------------------
data_JAGS <- list(
 # SUB-MODEL 1 DATA:
 Nstud.vac = Nstud.vac, rA.vac = rA.vac, rB.vac = rB.vac, 
                  nA.vac = nA.vac, nB.vac = nB.vac,
 # SUB-MODEL 2 DATA:
 mu.vac_alpha = mu.vac_alpha, mu.vac_beta = mu.vac_beta,
 # SUB-MODEL 3 DATA:
 iota.comp_alpha = iota.comp_alpha, iota.comp_beta = iota.comp_beta,
 # SUB-MODEL 4 DATA:
 chi.crossPro_mu = chi.crossPro_mu, chi.crossPro_sigma = chi.crossPro_sigma,
 # SUB-MODEL 5DATA:
 sigma.screening_alpha = sigma.screening_alpha, 
 sigma.screening_beta = sigma.screening_beta,
 # SUB-MODEL 6 DATA:
 omega.mu.log = omega.mu.log, omega.sigma.log = omega.sigma.log,
 # SUB-MODEL 7 DATA:
 delta.regr_alpha = delta.regr_alpha, delta.regr_beta = delta.regr_beta,
 # SUB-MODEL 8 DATA:
 delta.LSIL_alpha = delta.LSIL_alpha, delta.LSIL_beta = delta.LSIL_beta,
 delta.HSIL_alpha = delta.HSIL_alpha, delta.HSIL_beta = delta.HSIL_beta,
 # SUB-MODEL 9 DATA:
 pi.I.LSIL_Clear_alpha = pi.I.LSIL_Clear_alpha, pi.I.LSIL_Clear_beta = pi.I.LSIL_Clear_beta,
 pi.I.LSIL_HSIL_alpha = pi.I.LSIL_HSIL_alpha, pi.I.LSIL_HSIL_beta = pi.I.LSIL_HSIL_beta,
 pi.I.HSIL_Clear_alpha = pi.I.HSIL_Clear_alpha, pi.I.HSIL_Clear_beta = pi.I.HSIL_Clear_beta,
 pi.I.HSIL_LSIL_alpha = pi.I.HSIL_LSIL_alpha, pi.I.HSIL_LSIL_beta = pi.I.HSIL_LSIL_beta,
 pi.I.HSIL_Cancer_alpha = pi.I.HSIL_Cancer_alpha,
 pi.I.HSIL_Cancer_beta = pi.I.HSIL_Cancer_beta,
 
 # SUB-MODEL 10 DATA:
 cancer_alpha = cancer_alpha, 
 alphaI.year_1 = alphaI.year_1, betaI.year_1 = betaI.year_1,
 alphaI.year_2 = alphaI.year_2, betaI.year_2 = betaI.year_2, 
 alphaI.year_3 = alphaI.year_3, betaI.year_3 = betaI.year_3, 
 alphaI.year_4 = alphaI.year_4, betaI.year_4 = betaI.year_4,
 alphaII.year_1 = alphaII.year_1, betaII.year_1 = betaII.year_1, 
 alphaII.year_2 = alphaII.year_2, betaII.year_2 = betaII.year_2, 
 alphaII.year_3 = alphaII.year_3, betaII.year_3 = betaII.year_3, 
 alphaII.year_4 = alphaII.year_4, betaII.year_4 = betaII.year_4,
 alphaIII.year_1 = alphaIII.year_1, betaIII.year_1 = betaIII.year_1,
 alphaIII.year_2 = alphaIII.year_2, betaIII.year_2 = betaIII.year_2,
 alphaIII.year_3 = alphaII.year_3, betaIII.year_3 = betaIII.year_3,
 alphaIII.year_4 = alphaIII.year_4, betaIII.year_4 = betaIII.year_4,
 alphaIV.year_1 = alphaIV.year_1, betaIV.year_1 = betaIV.year_1,
 alphaIV.year_2 = alphaIV.year_2, betaIV.year_2 = betaIV.year_2,
 alphaIV.year_3 = alphaIV.year_3, betaIV.year_3 = betaIV.year_3,
 alphaIV.year_4 = alphaIV.year_4, betaIV.year_4 = betaIV.year_4,
 # SUB-MODEL 11 DATA:
 rHIV_pos = rHIV_pos, rHIV_neg = rHIV_neg, 
 nHIV_pos = nHIV_pos, nHIV_neg = nHIV_neg, 
 zeta_alpha = zeta_alpha, zeta_beta = zeta_beta
 
 )

# Set parameters to be monitored:
params <- c(
  # Vaccine efficacy probability and Odds:
  "OR.vac", "pEfficacy.vac",
  # Average vaccine coverage probability:
  "mu.vac_coverage",
  # Average vaccine compliance:
  "iota.compli",
  # Cross protection effect probability:
  "chi.xProtect",
  # Screening coverage probability:
  "sigma.screenCov",
  # Age-specific infection probability:
  "omega.age",
  # Age-specific regression from infection
  # to exposure:
  "delta.regr",
  # Progression from Infection to LSIL or HSIL:
  "delta.LSIL", "delta.HSIL",
  # Progression and regression from LSIL and HSIL:
  "pi.I.LSIL_Clear", "pi.I.LSIL_HSIL", "pi.I.HSIL_Clear", 
  "pi.I.HSIL_LSIL", "pi.I.HSIL_Cancer",
  # Cervical cancer parameters:
  "FIGO", 
  "FIGO_I.year_1", "FIGO_II.year_1", "FIGO_III.year_1", "FIGO_IV.year_1", 
  "FIGO_I.year_2", "FIGO_II.year_2", "FIGO_III.year_2", "FIGO_IV.year_2", 
  "FIGO_I.year_3", "FIGO_II.year_3", "FIGO_III.year_3", "FIGO_IV.year_3", 
  "FIGO_I.year_4", "FIGO_II.year_4", "FIGO_III.year_4", "FIGO_IV.year_4",
  # Relative Risk parameters:
  "HIV_odds", "HIV_RR", "zeta.HIVpop"
            )
# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 40000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)
# Run model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "HPV_model.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
options(max.print = 4000)
mod_JAGS
attach.jags(mod_JAGS)
