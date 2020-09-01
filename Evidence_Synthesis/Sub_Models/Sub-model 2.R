library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(parallel)

options(mc.cores = detectCores())

# The following sub-model simulates age-specific prevalence and average vaccine efficiacy 
# using data from multiple RCTs. Assumed any end-point - infection or disease - signalling 
# efficacy.

# ===========================================================================================
# Collated studies on vaccine efficacy -------------------------------------------------
# ===========================================================================================
# Study 1: Apter et al. (2015)   ------------------------------------------
# Efficacy of Human Papillomavirus 16 and 18 (HPV-16/18) AS04- Adjuvanted Vaccine against 
# Cervical Infection and Precancer in Young Women: Final Event-Driven Analysis of the 
# Randomized, Double- Blind PATRICIA Trial.

# Data collected from 6-month HPV 16/18 ATP cohort.

# Control (A):
# n = 5375
# cases = 435

# Vaccine (B):
# n = 5406
# cases = 32

# Study 2: Harper et al. (2004) -------------------------------------------
# Efficacy of a bivalent L1 virus-like particle vaccine in prevention of infection with human 
# papillomavirus types 16 and 18 in young women: a randomised controlled trial.

# Data collected from 27 month HPV 16/18 ATP cohort.

# Control (A):
# n = 355
# cases = 41

# Vaccine (B):
# n = 366
# cases = 12

# Study 3: Harper et al. (2006) -------------------------------------------
# Sustained efficacy up to 4·5 years of a bivalent L1 virus-like particle vaccine against 
# human papillomavirus types 16 and 18: follow-up from a randomised control trial.

# Data collected from 54 month HPV 16/18 ATP cohort.

# Control (A):
# n = 277
# cases = 28

# Vaccine (B):
# n = 310
# cases = 1

# Study 4: Herrero et al. (2011) ------------------------------------------
# Prevention of Persistent Human Papillomavirus Infection by an HPV16/18 Vaccine: A 
# Community-Based Randomized Clinical Trial in Guanacaste, Costa Rica.

# Data collected from 22-34 month follow-up HPV 16/18 ATP cohort.

# Control (A):
# n = 2239
# cases = 38

# Vaccine (B):
# n = 2190
# cases = 3

# Study 5: Herrero et al. (2013) ------------------------------------------
# Reduced Prevalence of Oral Human Papillomavirus (HPV) 4 Years after Bivalent HPV Vaccination
# in a Randomized Clinical Trial in Costa Rica.

# Data collected from 48 month 16/18 group. Note!: likely need to discount this study study has
# combined ATP and other cohorts.

# Control (A):
# n = 2924
# cases = 219

# Vaccine (B):
# n = 2910
# cases = 61

# Study 6: Konno et al. (2010) --------------------------------------------
# Efficacy of Human Papillomavirus Type 16/18 AS04YAdjuvanted Vaccine in Japanese Women Aged 
# 20 to 25 Years.

# Data collected from 6 month HPV 16/18 ATP cohort.

# Control (A):
# n = 392
# cases = 15

# Vaccine (B):
# n = 387
# cases = 0

# Study 7: Naud et al. (2014) ---------------------------------------------
# Sustained efficacy, immunogenicity, and safety of the HPV-16/18 AS04-adjuvanted vaccine: 
# Final analysis of a long-term follow-up study up to 9.4 years post-vaccination.

# Data collected from 12 month persistent infection HPV 16/18 ATP cohort, using the combined
# analysis from all follow-up studies, HPV-001/007/023.

# Control (A):
# n = 175
# cases = 10

# Vaccine (B):
# n = 193
# cases = 0

# Study 8: Paavonen et al. (2007) -----------------------------------------
# Efficacy of a prophylactic adjuvanted bivalent L1 virus-like-particle vaccine against 
# infection with human papillomavirus types 16 and 18 in young women: an interim analysis of a
# phase III double-blind, randomised controlled trial.

# Data collected from primary endpoint CIN2+ HPV 16/18.

# Control (A):
# n = 7838
# cases = 21

# Vaccine (B):
# n = 7788
# cases = 2

# Study 9: Paavonen et al. (2009) -----------------------------------------
# Efficacy of human papillomavirus (HPV)-16/18 AS04-adjuvanted vaccine against cervical 
# infection and precancer caused by oncogenic HPV types (PATRICIA): final analysis of a 
# double-blind, randomised study in young women.

# Data collected from 34·9 month HPB 16/18 ATP cohort. Indicator used CIN2+.

# Control (A):
# n = 7312
# cases = 56

# Vaccine (B):
# n = 7344
# cases = 4

# Study 10: Romanowski et al. (2009) --------------------------------------
# Sustained efficacy and immunogenicity of the human papillomavirus (HPV)-16/18 AS04-adjuvanted
# vaccine: analysis of a randomised placebo-controlled trial up to 6·4 years.

# Data collected from 12-month persistent infection HPV 16/18 ATP cohort.

# Control (A):
# n = 372
# cases = 20

# Vaccine (B):
# n = 401
# cases = 0

# Study 11: Zhu et al. (2014) ---------------------------------------------
# Efficacy, immunogenicity and safety of the HPV-16/18 AS04- adjuvanted vaccine in healthy 
# Chinese women aged 18–25 years: Results from a randomized controlled trial.

# Data collected from 6-Month PI and/or CIN1+ (primary endpoint) HPV 16/18 ATP cohort.

# Control (A):
# n = 2502
# cases = 17

# Vaccine (B):
# n = 2497
# cases = 1

# Vaccine efficacy data ----------------------------------------------------
# PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 11
rA.vac <- c(435, 41, 28, 38, 219, 15, 10, 21, 56, 20, 17)
nA.vac <- c(5375, 355, 277, 2239, 2924, 392, 175, 7838, 7312, 372, 2502)
rB.vac <- c(32, 12, 1, 3, 61, 0, 0, 2, 4, 0, 1)
nB.vac <- c(5406, 366, 310, 2190, 2910, 387, 193, 7788, 7344, 401, 2497)
# Check:
length(rA.vac) == length(nA.vac)
length(rB.vac) == length(nB.vac)

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

# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤25", "25-34", 
               "35-44", "45-54", "55-64")

# Estimated Prevalence:
Prevalence <- c(0.09516258, 0.1130796, 0.139292, 0.1563352, 0.139292, 
                .435, .3643, 0.2051, 0.1852, 0.1654)
mu.a.log <- lnorm_params(m = Prevalence, v = .1)$mu
sigma.a.log <- 1 / (lnorm_params(m = Prevalence, v = .1)$sigma) ^ 2
mu.a.log
sigma.a.log
cbind(age_group, Prevalence)
# Informative prior sampling model:
# omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])

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
alpha.c <- list(.55, .15, .12, .18)

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

# ==========================================================================================
# Sub-models 3 ------------------------------------------------------
# ==========================================================================================
model_String <- "model {
# SUB-MODEL 1: AGE-SPECIFIC PREVALENCE,
  # model parameters abbreviated by .age
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior.
  for (i in 1:10) {
    omega.age[i] ~ dlnorm(mu.a.log[i], sigma.a.log[i])
    }
  
  
# SUB-MODEL 2: VACCINE-EFFICACY,
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
  
  # Priors for sub-model 2:
  psi.vac ~ dnorm(0, 1.0e-1)
  tau.vac ~ dunif(0, 1)
  prec.vac <- 1 / pow(tau.vac, 2)
  # Odds ratio:
  OR.vac <- exp(psi.vac)
  # Probability of no protection given
  # vaccine:
  pEfficacy.vac <- (OR.vac / 1 + OR.vac)
  

# SUB-MODEL 3: CERVICAL CANCER,
# the distribution of invidiuals who have 
# cancer according to FIGO stages I-IV, and 
# the associated survival probabilities according 
# to each stage, over 4 years.

  # Note: equivalent to standard PSA, as it is 
  # sampling direclty from the prior.
    
# Distribution of Cervical Cancer stages: 
 FIGO ~ ddirch(alpha.c)
 
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
 
}
"
writeLines(text = model_String, con = "Age_and_Efficacy.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
 Nstud.vac = Nstud.vac, rA.vac = rA.vac, rB.vac = rB.vac, 
                  nA.vac = nA.vac, nB.vac = nB.vac,
 # Population prevalence data:
 mu.a.log = mu.a.log, sigma.a.log = sigma.a.log,
 # Cancer data:
 alpha.c = alpha.c, 
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
 alphaIV.year_4 = alphaIV.year_4, betaIV.year_4 = betaIV.year_4
 )

# Initial JAGS sampler values:
# Parameters to monitor:
params <- c(
  # Vaccine efficacy parameters:
  "OR.vac", "pEfficacy.vac", "omega.age", 
  # Cervical cancer parameters:
  "FIGO", 
  "FIGO_I.year_1", "FIGO_II.year_1", "FIGO_III.year_1", "FIGO_IV.year_1", 
  "FIGO_I.year_2", "FIGO_II.year_2", "FIGO_III.year_2", "FIGO_IV.year_2", 
  "FIGO_I.year_3", "FIGO_II.year_3", "FIGO_III.year_3", "FIGO_IV.year_3", 
  "FIGO_I.year_4", "FIGO_II.year_4", "FIGO_III.year_4", "FIGO_IV.year_4"
            )
# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 40000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "Age_and_Efficacy.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS
# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("pEfficacy.vac", "omega.age[1]", "omega.age[2]", 
                               "omega.age[3]", "omega.age[4]", "omega.age[5]", 
                               "FIGO_IV.year_3", "FIGO[1]", "FIGO[3]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , ], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("pEfficacy.vac", "omega.age[1]", "omega.age[2]", 
                               "omega.age[3]", "omega.age[4]", "omega.age[5]", 
                               "FIGO_III.year_4", "FIGO[1]", "FIGO[3]"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("pEfficacy.vac", "omega.age[1]", "omega.age[2]", 
                               "omega.age[3]", "omega.age[4]", "omega.age[5]",
                             "FIGO_II.year_3", "FIGO[1]", "FIGO[3]"),
         lags = 150)

# Calibration: ------------------------------------------------------------
# Age-specific HPV prevalence:
plot(apply(omega.age, 2, mean), type = "l", ylab = "Prevalence", xlab = "By Age",
     col = "red", lwd = 2)
# Distribution of cancer stages:
plot(apply(FIGO, 2, mean), type = "l")
# Vaccine efficacy probability:
1 - (apply(pEfficacy.vac, 2, mean))
# Example: say 21% of some age pop. group has HPV:
.21 * apply(pEfficacy.vac, 2, mean)

# End file ----------------------------------------------------------------