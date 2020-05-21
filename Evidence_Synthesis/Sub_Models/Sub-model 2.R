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
# Vaccine efficacy --------------------------------------------------------
# ===========================================================================================
# Vaccine efficacy is expressed as the reduction in the occurrence of HPV due to vaccination.
# Individuals who are vaccinated fully experience a rate of occurrence of HPV that is 
# (1 – alpha) times that of those who are not vaccinated (Favato G., et al. 2012).

# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).

# Study 1: Munoz N., et al. 2009 ------------------------------------------
# Safety, immunogenicity, and eðcacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) recombinant vaccine in women aged 24–45 years: a randomised, double-blind trial.

# Women aged 24–45 years with no history of genital warts or cervical disease were enrolled 
# from community health centres, academic health centres, and primary health-care providers 
# into an ongoing multicentre, parallel, randomised, placebo-controlled, double-blind study.

# Efficacy against the incidence of infection .

# Events by group:
# tA:41
# nA:1607
# tB:4
# nB:1615

 # Study 2: Villa L, et al. 2006 -------------------------------------------
# High sustained efficacy of a prophylactic quadrivalent human papillomavirus types 
# 6/11/16/18 L1 virus-like particle vaccine through 5 years of follow-up.

# The study enrolled nonpregnant, healthy women who had no prior abnormal Pap smears, and 
# reported a lifetime history of four or fewer male sex partners. By end point persistent 
# infection.

# Events by group:
# tA:46
# nA:233
# tB:2
# nB:235

# Study 3: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years. By end point persistent infection.

# Events by group:
# tA:38
# nA:1583
# tB:1
# nB:1578

# Study 4: Wei L., et al. 2018 --------------------------------------------
# Efficacy of quadrivalent human papillomavirus vaccine against persistent infection and 
# genital disease in Chinese women: A randomized, placebo-controlled trial with 78-month 
# follow-up.

# Pregnant women and those with a history of genital warts or significant cervical disease, 
# active cervical disease, or prior HPV vaccine recipients were excluded. 12-Month cervical
# Persisitent Infection.

# Events by group:
# tA:53
# nA:1245
# tB:5
# nB:1276

# Study 5: Yoshikawa H., et al. 2013 --------------------------------------
# Efficacy of quadrivalent human papillomavirus (types 6, 11, 16 and 18) vaccine (GARDASIL) 
# in Japanese women aged 18–26 years.

# A randomized double-blind placebo-controlled phase II trial was conducted to evaluate the 
# efficacy of a prophylactic quadrivalent vaccine targeting the human papillomavirus (HPV) 
# types most fre- quently associated with cervical cancer (types 16 ⁄ 18) and genital warts 
# (types 6 ⁄ 11) in Japanese women aged 18–26 years.

# Events by group:
# tA:24
# nA:422
# tB:3
# nB:419

# Study 6: Garland S., et al. Future I. 2007 --------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# A phase 3 trial was conducted to evaluate the efficacy of a prophylactic quadrivalent 
# vaccine in preventing anogenital diseases associated with human papillomavirus (HPV) types
# 6, 11, 16, and 18.

# Events by group:
# tA:60
# nA:2279
# tB:0
# nB:2261

# Study 7: Koutsky L., et al. Future II. 2007 -----------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent High-Grade Cervical Lesions.

# Randomized, double-blind trial, we assigned 12,167 women between the ages of 15 and 26 
# years to receive three doses of either HPV-6/11/16/18 vaccine or placebo, administered at 
# day 1, month 2, and month 6. The primary analysis was performed for a per-protocol 
# susceptible population that included 5305 women in the vaccine group and 5260 in the 
# placebo group who had no virologic evidence of infection with HPV-16 or HPV-18 through 1 
# month after the third dose (month 7).

# Events by group:
# tA:42
# nA:5260
# tB:1
# nB:5305

# Vaccine efficacy data ----------------------------------------------------
# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 7
rA.vac <- c(41, 46, 38, 53, 24, 60, 42)
nA.vac <- c(1607, 233, 1583, 1245, 422, 2279, 5260)
rB.vac <- c(4, 2, 1, 5, 3, 0, 1)
nB.vac <- c(1615, 235, 1578, 1276, 419, 2261, 5305)

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
# SUB-MODEL 2: AGE-SPECIFIC PREVALENCE,
  # model parameters abbreviated by .age
    # Note: equivalent to standard PSA, as it is 
    # sampling direclty from the prior.
  for (i in 1:10) {
    omega.age[i] ~ dlnorm(mu.a.log[i], sigma.a.log[i])
    }
  
  
# SUB-MODEL 3: VACCINE-EFFICACY,
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
  
  # Priors for sub-model 3:
  psi.vac ~ dnorm(0, 1.0e-1)
  tau.vac ~ dunif(0, 1)
  prec.vac <- 1 / pow(tau.vac, 2)
  # Odds ratio:
  OR.vac <- exp(psi.vac)
  # Probability of no protection given
  # vaccine:
  pEfficacy.vac <- (OR.vac / 1 + OR.vac)
  

# SUB-MODEL 4: CERVICAL CANCER,
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